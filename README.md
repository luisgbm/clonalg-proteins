# clonalg-proteins

A sequential Rust implementation of the Clonal Selection Algorithm (Clonalg)
for protein folding in the three-dimensional hydrophobic-polar lattice model
(3D HP).

**[Explore the interactive 3D protein-folding visualization](https://luisgbm.github.io/clonalg-proteins/protein-folding-3d.html)**

The project models a protein as a sequence of hydrophobic (`H`) and polar (`P`)
residues and searches for a low-energy, self-avoiding conformation on a cubic
lattice. Each non-consecutive pair of adjacent hydrophobic residues contributes
`-1` to the conformation's energy, so maximizing H-H contacts minimizes energy.

## Protein folding problem

Proteins are chains of amino-acid residues whose biological behavior depends
strongly on their three-dimensional conformation. Although the amino-acid
sequence describes the protein's primary structure, the chain can assume a very
large number of spatial arrangements. Protein structure prediction asks which
conformation is energetically most favorable for a given sequence.

Experimental structure determination can be expensive, while an exhaustive
computational search quickly becomes impractical as the sequence grows. A
simulation therefore needs both a simplified physical model and a search method
capable of exploring promising conformations without enumerating the complete
space.

### The 3D HP model

The hydrophobic-polar model reduces the twenty amino acids to two classes:

- **Hydrophobic (`H`)** residues tend to avoid water and cluster inside the
  folded protein.
- **Polar (`P`)** residues interact more readily with the surrounding aqueous
  environment.

A conformation is represented as a walk on a three-dimensional cubic lattice.
Every residue occupies one lattice coordinate and each consecutive pair must be
one unit apart. Two residues cannot occupy the same coordinate, making valid
conformations self-avoiding.

The model assigns one favorable contact for each pair of hydrophobic residues
that are adjacent in space but not consecutive in the protein chain:

```text
energy = -(number of non-local H-H contacts)
```

Consequently, a more negative energy denotes a better fold. This is an
abstraction rather than a complete molecular-energy model, but it preserves the
hydrophobic effect as the primary force guiding the simulated folding process.

## Clonal selection approach

Clonalg is inspired by the adaptive immune system's clonal selection principle:
an antibody with high affinity for an antigen is selected, cloned, and mutated
to produce related antibodies that may have still greater affinity.

The protein-folding simulation maps these concepts as follows:

| Immune-system concept | Protein-folding representation |
|---|---|
| Antigen | The fixed H/P protein sequence |
| Antibody | A candidate lattice conformation |
| Affinity | Number of non-local H-H contacts |
| Cloning | Copies of candidate conformations |
| Somatic hypermutation | Position changes controlled by current affinity |
| Selection | Survival of conformations with more H-H contacts |
| Repertoire diversity | Unique candidates and random replacements |

Lower-affinity candidates receive a larger mutation allowance, encouraging
exploration, while high-affinity candidates are changed more conservatively.
Hypermacromutation complements this process by perturbing a randomly selected
interval of the conformation. Both operators use the
**first-constructive-mutation rule**: mutations are attempted until the first
valid self-avoiding conformation is found; if no constructive change is found
within the operator's limit, the original candidate is retained.

### Simulation cycle

Each independent simulation performs the following steps:

1. Generate a random population of feasible, self-avoiding conformations.
2. Evaluate affinity by counting every non-local H-H contact.
3. Clone each population member according to the configured clone factor.
4. Apply affinity-dependent hypermutation to one clone population.
5. Apply interval hypermacromutation to a second clone population.
6. Combine the parents and both offspring populations.
7. Rank the candidates, remove duplicate conformations, and retain the best.
8. Introduce random feasible candidates to restore population size and
   diversity.
9. Repeat until the generation limit or an optional contact target is reached.

Because initialization and mutation are stochastic, a single run is not enough
to characterize the algorithm. Reproduction experiments use independent seeds
and report the mean, sample standard deviation, minimum, maximum, optimum
frequency, score distribution, and mean runtime across all runs.

## Scientific background

This implementation follows the research developed by **Carolina Paula de
Almeida** in her 2007 master's dissertation, *Aplicação de Sistemas
Imunológicos Artificiais para a Predição da Estrutura de Proteínas*. Her work
applied Artificial Immune Systems to protein structure prediction in the 3D HP
model and studied Clonalg, immune networks, feasible and penalized infeasible
search spaces, aging operators, affinity maturation, fuzzy inference, and Tabu
Search. The dissertation found Clonalg-based models more effective than the
evaluated immune-network models, with its strongest results combining
penalized search, fuzzy aging, and intensive maturation through Tabu Search.

The later PIBIC work by **Luís Guilherme Bergamini Mendes**, supervised by
**Prof. Dr. Myriam Regattieri B. S. Delgado**, was presented as *Computação
Paralela: Suporte para Pesquisas de Computação Natural com Aplicação em
Bioinformática*. It adapted the Java Clonalg implementation for sequential and
MPI-based execution, then evaluated parallel speedup and migration between
independently evolving populations. For the 20-residue benchmark, the
sequential experiment performed 96 independent runs of 200 generations and
reported a mean contact score of `9.79`, a standard deviation of `1.21`, and a
known optimum of `11` contacts (energy `-11`).

This Rust project focuses exclusively on the **sequential Clonalg search**. It
does not implement MPI, clusters, island populations, or migration.

## 3D visualization

Open [`protein-folding-3d.html`](protein-folding-3d.html) directly in a modern
browser to watch the benchmark protein change conformation until it reaches a
valid energy `-11` fold.

The viewer is a single offline HTML file with Three.js and OrbitControls
embedded. It requires no web server, package installation, build step, or
network connection. It provides:

- real-time 3D residue and bond animation;
- hydrophobic, polar, and non-local contact highlighting;
- live generation, energy, contact, and search-state metrics;
- orbit, zoom, play/pause, reset, progress, and speed controls;
- responsive desktop and mobile layouts.

The publication reports aggregate scores and the final optimum, but does not
publish the coordinate history of an individual run. The animation is therefore
a deterministic Clonalg-inspired reconstruction, not a claim that these were
the article's exact intermediate conformations. Its terminal coordinates form a
self-avoiding fold produced by this Rust implementation with 11 non-local H-H
contacts.

## Implemented model

- Six-direction cubic lattice (`±X`, `±Y`, and `±Z`)
- Self-avoiding protein conformations
- Non-local hydrophobic-contact energy function
- Random feasible population initialization
- Clonal expansion
- Affinity-dependent hypermutation
- Interval hypermacromutation
- First-constructive-mutation feasibility rule
- Selection across parents and both offspring populations
- Duplicate-conformation removal
- Random replacement to preserve population diversity
- Deterministic seeded runs
- Early stopping at a requested contact target
- Multi-run statistical experiments

The default protein is the 20-residue Tortilla benchmark:

```text
HPHPPHHPHPPHPHHPPHPH
```

Its known optimum is `11` H-H contacts, corresponding to energy `-11`.

## PIBIC reproduction

The Rust experiment uses the same benchmark shape as the PIBIC sequential
evaluation: 96 independent runs, 200 generations, population size 10, clone
factor 2, and hypermutation rate 40%.

Run it with:

```console
cargo run --release -- --runs 96 --seed 0 --generations 200
```

Results measured on the current implementation:

| Metric | Rust reproduction | PIBIC article |
|---|---:|---:|
| Independent runs | 96 | 96 |
| Generations per run | 200 | 200 |
| Mean H-H contacts | **10.094** | **9.79** |
| Sample standard deviation | 0.834 | 1.21 |
| Minimum contacts | 8 | Not reported |
| Maximum contacts | **11** | **11** |
| Best energy | **-11** | **-11** |
| Runs reaching the optimum | 32 (33.33%) | Not reported |
| Mean runtime per run | 67.552 ms | 26.49 s |

The Rust solver reproduces the known optimum and slightly exceeds the article's
reported mean solution quality. Runtime values are not directly comparable:
the PIBIC measurements used a Java implementation on 2009-era 2.4 GHz
hardware, whereas the Rust measurement was made on current hardware.
Furthermore, the publication does not provide its complete seed list, exact
replacement count, or all low-level implementation choices.

For seeds `0..=95`, the Rust score distribution was:

| H-H contacts | Runs |
|---:|---:|
| 8 | 6 |
| 9 | 11 |
| 10 | 47 |
| 11 | 32 |

## Usage

Run one fold using the default benchmark and parameters:

```console
cargo run --release
```

Use a custom HP sequence and reproducible seed:

```console
cargo run --release -- --sequence HHHHPPHH --seed 7 --generations 500
```

Stop when a known target is reached:

```console
cargo run --release -- --target-contacts 11 --generations 5000
```

List every command-line option:

```console
cargo run -- --help
```

A single run prints the best energy, number of H-H contacts, generation and
evaluation counts, lattice moves, and coordinates for every residue. A
multi-run experiment prints aggregate statistics and the score distribution.

## Development

Format, test, and lint the project with:

```console
cargo fmt --check
cargo test
cargo clippy --all-targets -- -D warnings
```
