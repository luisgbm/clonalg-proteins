use std::collections::HashSet;

use rand::SeedableRng;
use rand::rngs::StdRng;

use crate::hp::{Candidate, Protein};

#[derive(Clone, Debug)]
pub struct Config {
    pub population_size: usize,
    pub clone_factor: usize,
    pub mutation_rate: f64,
    pub random_injection_rate: f64,
    pub generations: usize,
    pub optimal_contacts: Option<usize>,
    pub target_contacts: Option<usize>,
    pub seed: u64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            population_size: 10,
            clone_factor: 2,
            mutation_rate: 0.40,
            random_injection_rate: 0.10,
            generations: 200,
            optimal_contacts: Some(11),
            target_contacts: None,
            seed: 42,
        }
    }
}

impl Config {
    pub fn validate(&self) -> Result<(), String> {
        if self.population_size == 0 {
            return Err("population size must be greater than zero".into());
        }
        if self.clone_factor == 0 {
            return Err("clone factor must be greater than zero".into());
        }
        if !(0.0..=1.0).contains(&self.mutation_rate) {
            return Err("mutation rate must be between 0 and 1".into());
        }
        if !(0.0..1.0).contains(&self.random_injection_rate) {
            return Err("random injection rate must be at least 0 and less than 1".into());
        }
        if self.generations == 0 {
            return Err("generations must be greater than zero".into());
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct RunResult {
    pub best: Candidate,
    pub contacts: usize,
    pub energy: i32,
    pub generations: usize,
    pub evaluations: usize,
}

pub struct Clonalg {
    config: Config,
}

impl Clonalg {
    pub fn new(config: Config) -> Result<Self, String> {
        config.validate()?;
        Ok(Self { config })
    }

    pub fn run(&self, protein: &Protein) -> RunResult {
        let mut rng = StdRng::seed_from_u64(self.config.seed);
        let mut population: Vec<Candidate> = (0..self.config.population_size)
            .map(|_| Candidate::random(protein.len(), &mut rng))
            .collect();
        let mut evaluations = population.len();
        rank(&mut population, protein);

        let mut completed_generations = 0;
        for generation in 1..=self.config.generations {
            if reached_target(&population[0], protein, self.config.target_contacts) {
                break;
            }

            let offspring_capacity = population.len() * self.config.clone_factor;
            let mut hypermutated = Vec::with_capacity(offspring_capacity);
            let mut hypermacromutated = Vec::with_capacity(offspring_capacity);
            for antibody in &population {
                let max_mutations = self.max_mutations(antibody, protein);
                for _ in 0..self.config.clone_factor {
                    hypermutated.push(antibody.hypermutate(max_mutations, &mut rng));
                    hypermacromutated.push(antibody.hypermacromutate(&mut rng));
                }
            }
            evaluations += hypermutated.len() + hypermacromutated.len();

            population.append(&mut hypermutated);
            population.append(&mut hypermacromutated);
            rank(&mut population, protein);
            population = select_unique(population, protein, self.survivor_count());

            while population.len() < self.config.population_size {
                population.push(Candidate::random(protein.len(), &mut rng));
                evaluations += 1;
            }
            rank(&mut population, protein);
            completed_generations = generation;
        }

        let best = population.remove(0);
        let contacts = best.contacts(protein);
        RunResult {
            energy: -(contacts as i32),
            best,
            contacts,
            generations: completed_generations,
            evaluations,
        }
    }

    fn survivor_count(&self) -> usize {
        let injected = ((self.config.population_size as f64) * self.config.random_injection_rate)
            .ceil() as usize;
        self.config.population_size.saturating_sub(injected).max(1)
    }

    fn max_mutations(&self, antibody: &Candidate, protein: &Protein) -> usize {
        let alpha = self.config.mutation_rate * protein.len() as f64;
        let affinity = antibody.contacts(protein);
        match (self.config.optimal_contacts, affinity) {
            (Some(optimum), affinity) if affinity > 0 => {
                ((1.0 + optimum as f64 / affinity as f64) * alpha).ceil() as usize
            }
            (Some(optimum), _) => ((2 + optimum) as f64 * alpha).ceil() as usize,
            (None, _) => alpha.ceil() as usize,
        }
        .max(1)
    }
}

fn reached_target(candidate: &Candidate, protein: &Protein, target: Option<usize>) -> bool {
    target.is_some_and(|target| candidate.contacts(protein) >= target)
}

fn rank(population: &mut [Candidate], protein: &Protein) {
    population.sort_unstable_by_key(|candidate| std::cmp::Reverse(candidate.contacts(protein)));
}

fn select_unique(population: Vec<Candidate>, protein: &Protein, count: usize) -> Vec<Candidate> {
    let mut seen = HashSet::with_capacity(count);
    let mut selected = Vec::with_capacity(count);

    for candidate in population {
        if seen.insert(candidate.clone()) {
            selected.push(candidate);
            if selected.len() == count {
                break;
            }
        }
    }
    rank(&mut selected, protein);
    selected
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_configuration() {
        let config = Config {
            population_size: 0,
            ..Config::default()
        };
        assert!(Clonalg::new(config).is_err());
    }

    #[test]
    fn run_is_reproducible_for_a_seed() {
        let protein = Protein::parse("HHPPHHPPHH").unwrap();
        let solver = Clonalg::new(Config {
            generations: 20,
            seed: 123,
            ..Config::default()
        })
        .unwrap();

        let first = solver.run(&protein);
        let second = solver.run(&protein);

        assert_eq!(first.best, second.best);
        assert_eq!(first.energy, second.energy);
        assert_eq!(first.evaluations, second.evaluations);
    }

    #[test]
    fn stops_when_target_is_reached() {
        let protein = Protein::parse("HP").unwrap();
        let solver = Clonalg::new(Config {
            target_contacts: Some(0),
            ..Config::default()
        })
        .unwrap();

        let result = solver.run(&protein);
        assert_eq!(result.generations, 0);
        assert_eq!(result.evaluations, 10);
    }
}
