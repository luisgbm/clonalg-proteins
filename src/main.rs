use std::process::ExitCode;
use std::time::{Duration, Instant};

use clap::Parser;
use clonalg_proteins::{Clonalg, Config, Protein};

const TORTILLA_20: &str = "HPHPPHHPHPPHPHHPPHPH";

#[derive(Debug, Parser)]
#[command(version, about)]
struct Args {
    /// Protein sequence using H (hydrophobic) and P (polar).
    #[arg(short, long, default_value = TORTILLA_20)]
    sequence: String,

    /// Number of antibodies retained in each generation.
    #[arg(long, default_value_t = 10)]
    population: usize,

    /// Number of mutated clones produced per antibody.
    #[arg(long, default_value_t = 2)]
    clone_factor: usize,

    /// Maximum mutation fraction, applied most strongly to low-affinity antibodies.
    #[arg(long, default_value_t = 0.40)]
    mutation_rate: f64,

    /// Fraction of each generation replaced by random antibodies.
    #[arg(long, default_value_t = 0.10)]
    random_injection_rate: f64,

    /// Maximum number of generations.
    #[arg(short, long, default_value_t = 200)]
    generations: usize,

    /// Known optimum contact count, used by affinity-dependent hypermutation.
    #[arg(long, default_value_t = 11)]
    optimal_contacts: usize,

    /// Stop once this many hydrophobic contacts are found.
    #[arg(long)]
    target_contacts: Option<usize>,

    /// Seed for reproducible pseudo-random evolution.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Number of independent runs; subsequent runs increment the seed.
    #[arg(long, default_value_t = 1)]
    runs: usize,
}

fn main() -> ExitCode {
    let args = Args::parse();
    match execute(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn execute(args: Args) -> Result<(), String> {
    let protein = Protein::parse(&args.sequence).map_err(|error| error.to_string())?;
    if protein.len() < 2 {
        return Err("protein sequence must contain at least two residues".into());
    }
    if args.runs == 0 {
        return Err("runs must be greater than zero".into());
    }

    let config = Config {
        population_size: args.population,
        clone_factor: args.clone_factor,
        mutation_rate: args.mutation_rate,
        random_injection_rate: args.random_injection_rate,
        generations: args.generations,
        optimal_contacts: Some(args.optimal_contacts),
        target_contacts: args.target_contacts,
        seed: args.seed,
    };

    if args.runs > 1 {
        return run_experiment(&args, &protein, config);
    }

    let result = Clonalg::new(config)?.run(&protein);

    println!("Sequence:    {}", args.sequence.to_ascii_uppercase());
    println!("Length:      {}", protein.len());
    println!("Energy:      {}", result.energy);
    println!("H-H contacts: {}", result.contacts);
    println!("Generations: {}", result.generations);
    println!("Evaluations: {}", result.evaluations);
    println!(
        "Moves:       {}",
        result
            .best
            .moves()
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!("Coordinates:");
    for (index, (residue, point)) in protein
        .residues()
        .iter()
        .zip(result.best.coordinates())
        .enumerate()
    {
        let residue = match residue {
            clonalg_proteins::Residue::Hydrophobic => 'H',
            clonalg_proteins::Residue::Polar => 'P',
        };
        println!(
            "  {:>3} {} ({:>3}, {:>3}, {:>3})",
            index + 1,
            residue,
            point.x,
            point.y,
            point.z
        );
    }

    Ok(())
}

fn run_experiment(args: &Args, protein: &Protein, mut config: Config) -> Result<(), String> {
    let mut contacts = Vec::with_capacity(args.runs);
    let mut elapsed = Duration::ZERO;

    for offset in 0..args.runs {
        config.seed = args
            .seed
            .checked_add(offset as u64)
            .ok_or("seed overflow across experiment runs")?;
        let solver = Clonalg::new(config.clone())?;
        let started = Instant::now();
        let result = solver.run(protein);
        elapsed += started.elapsed();
        contacts.push(result.contacts);
    }

    let mean = contacts.iter().sum::<usize>() as f64 / contacts.len() as f64;
    let sample_variance = if contacts.len() > 1 {
        contacts
            .iter()
            .map(|value| (*value as f64 - mean).powi(2))
            .sum::<f64>()
            / (contacts.len() - 1) as f64
    } else {
        0.0
    };
    let minimum = contacts.iter().min().copied().unwrap_or_default();
    let maximum = contacts.iter().max().copied().unwrap_or_default();
    let optimum_count = contacts
        .iter()
        .filter(|contacts| **contacts >= args.optimal_contacts)
        .count();

    println!(
        "Sequence:              {}",
        args.sequence.to_ascii_uppercase()
    );
    println!("Runs:                  {}", contacts.len());
    println!("Seed range:            {}..={}", args.seed, config.seed);
    println!("Mean H-H contacts:     {mean:.3}");
    println!("Sample standard dev.:  {:.3}", sample_variance.sqrt());
    println!("Minimum contacts:      {minimum}");
    println!("Maximum contacts:      {maximum}");
    println!(
        "Known optimum reached: {} ({:.2}%)",
        optimum_count,
        optimum_count as f64 * 100.0 / contacts.len() as f64
    );
    println!(
        "Mean runtime:          {:.3} ms",
        elapsed.as_secs_f64() * 1_000.0 / contacts.len() as f64
    );

    println!("Distribution:");
    for score in minimum..=maximum {
        let count = contacts
            .iter()
            .filter(|contacts| **contacts == score)
            .count();
        if count > 0 {
            println!("  {score:>2} contacts: {count:>3}");
        }
    }

    Ok(())
}
