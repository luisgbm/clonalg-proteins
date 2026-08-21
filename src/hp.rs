use std::collections::HashSet;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use rand::Rng;
use rand::seq::SliceRandom;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Residue {
    Hydrophobic,
    Polar,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Protein {
    residues: Vec<Residue>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProteinParseError {
    position: usize,
    character: char,
}

impl Protein {
    pub fn parse(sequence: &str) -> Result<Self, ProteinParseError> {
        let mut residues = Vec::new();

        for (position, character) in sequence.chars().enumerate() {
            let residue = match character {
                'H' | 'h' => Residue::Hydrophobic,
                'P' | 'p' => Residue::Polar,
                character if character.is_whitespace() => continue,
                character => {
                    return Err(ProteinParseError {
                        position,
                        character,
                    });
                }
            };
            residues.push(residue);
        }

        Ok(Self { residues })
    }

    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }
}

impl Display for ProteinParseError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "invalid residue '{}' at character {}; expected H or P",
            self.character,
            self.position + 1
        )
    }
}

impl Error for ProteinParseError {}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Point {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl Point {
    pub const ORIGIN: Self = Self { x: 0, y: 0, z: 0 };

    fn step(self, direction: Direction) -> Self {
        let (x, y, z) = direction.delta();
        Self {
            x: self.x + x,
            y: self.y + y,
            z: self.z + z,
        }
    }

    fn is_adjacent(self, other: Self) -> bool {
        (self.x - other.x).abs() + (self.y - other.y).abs() + (self.z - other.z).abs() == 1
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Direction {
    PosX,
    NegX,
    PosY,
    NegY,
    PosZ,
    NegZ,
}

impl Direction {
    pub const ALL: [Self; 6] = [
        Self::PosX,
        Self::NegX,
        Self::PosY,
        Self::NegY,
        Self::PosZ,
        Self::NegZ,
    ];

    fn delta(self) -> (i32, i32, i32) {
        match self {
            Self::PosX => (1, 0, 0),
            Self::NegX => (-1, 0, 0),
            Self::PosY => (0, 1, 0),
            Self::NegY => (0, -1, 0),
            Self::PosZ => (0, 0, 1),
            Self::NegZ => (0, 0, -1),
        }
    }
}

impl Display for Direction {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::PosX => "+X",
            Self::NegX => "-X",
            Self::PosY => "+Y",
            Self::NegY => "-Y",
            Self::PosZ => "+Z",
            Self::NegZ => "-Z",
        })
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Candidate {
    moves: Vec<Direction>,
}

impl Candidate {
    pub fn from_moves(moves: Vec<Direction>) -> Option<Self> {
        let candidate = Self { moves };
        candidate.is_self_avoiding().then_some(candidate)
    }

    pub fn random<R: Rng + ?Sized>(protein_len: usize, rng: &mut R) -> Self {
        if protein_len <= 1 {
            return Self { moves: Vec::new() };
        }

        loop {
            if let Some(moves) = grow_chain(protein_len - 1, &[Point::ORIGIN], rng) {
                return Self { moves };
            }
        }
    }

    pub fn moves(&self) -> &[Direction] {
        &self.moves
    }

    pub fn coordinates(&self) -> Vec<Point> {
        let mut coordinates = Vec::with_capacity(self.moves.len() + 1);
        let mut current = Point::ORIGIN;
        coordinates.push(current);

        for direction in &self.moves {
            current = current.step(*direction);
            coordinates.push(current);
        }

        coordinates
    }

    pub fn contacts(&self, protein: &Protein) -> usize {
        debug_assert_eq!(self.moves.len() + 1, protein.len());
        let coordinates = self.coordinates();
        let mut contacts = 0;

        for left in 0..protein.len() {
            if protein.residues[left] != Residue::Hydrophobic {
                continue;
            }
            for right in (left + 2)..protein.len() {
                if protein.residues[right] == Residue::Hydrophobic
                    && coordinates[left].is_adjacent(coordinates[right])
                {
                    contacts += 1;
                }
            }
        }

        contacts
    }

    pub fn energy(&self, protein: &Protein) -> i32 {
        -(self.contacts(protein) as i32)
    }

    pub fn hypermutate<R: Rng + ?Sized>(&self, max_mutations: usize, rng: &mut R) -> Self {
        if self.moves.is_empty() {
            return self.clone();
        }

        let mut moves = self.moves.clone();
        for _ in 0..max_mutations.max(1) {
            let position = rng.random_range(0..moves.len());
            moves[position] = Self::random_other_direction(moves[position], rng);
            if let Some(candidate) = Self::from_moves(moves.clone()) {
                return candidate;
            }
        }

        self.clone()
    }

    pub fn hypermacromutate<R: Rng + ?Sized>(&self, rng: &mut R) -> Self {
        if self.moves.is_empty() {
            return self.clone();
        }

        let start = rng.random_range(0..self.moves.len());
        let end = rng.random_range(start..self.moves.len());
        let forward = rng.random_bool(0.5);
        let positions: Box<dyn Iterator<Item = usize>> = if forward {
            Box::new(start..=end)
        } else {
            Box::new((start..=end).rev())
        };
        let mut moves = self.moves.clone();

        for position in positions {
            moves[position] = Self::random_other_direction(moves[position], rng);
            if let Some(candidate) = Self::from_moves(moves.clone()) {
                return candidate;
            }
        }

        self.clone()
    }

    fn is_self_avoiding(&self) -> bool {
        let coordinates = self.coordinates();
        coordinates.iter().collect::<HashSet<_>>().len() == coordinates.len()
    }

    fn random_other_direction<R: Rng + ?Sized>(current: Direction, rng: &mut R) -> Direction {
        loop {
            let direction = Direction::ALL[rng.random_range(0..Direction::ALL.len())];
            if direction != current {
                return direction;
            }
        }
    }
}

fn grow_chain<R: Rng + ?Sized>(
    remaining: usize,
    occupied_prefix: &[Point],
    rng: &mut R,
) -> Option<Vec<Direction>> {
    let mut occupied: HashSet<Point> = occupied_prefix.iter().copied().collect();
    let mut current = *occupied_prefix.last()?;
    let mut moves = Vec::with_capacity(remaining);

    for _ in 0..remaining {
        let mut directions = Direction::ALL;
        directions.shuffle(rng);

        let (direction, next) = directions.into_iter().find_map(|direction| {
            let next = current.step(direction);
            (!occupied.contains(&next)).then_some((direction, next))
        })?;

        moves.push(direction);
        occupied.insert(next);
        current = next;
    }

    Some(moves)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn parses_hp_sequences_case_insensitively() {
        let protein = Protein::parse("Hp pH").unwrap();
        assert_eq!(
            protein.residues(),
            &[
                Residue::Hydrophobic,
                Residue::Polar,
                Residue::Polar,
                Residue::Hydrophobic
            ]
        );
    }

    #[test]
    fn rejects_unknown_residues() {
        let error = Protein::parse("HPX").unwrap_err();
        assert_eq!(
            error.to_string(),
            "invalid residue 'X' at character 3; expected H or P"
        );
    }

    #[test]
    fn counts_only_non_local_hydrophobic_contacts() {
        let protein = Protein::parse("HHHH").unwrap();
        let square =
            Candidate::from_moves(vec![Direction::PosX, Direction::PosY, Direction::NegX]).unwrap();

        assert_eq!(square.contacts(&protein), 1);
        assert_eq!(square.energy(&protein), -1);
    }

    #[test]
    fn rejects_colliding_chains() {
        assert!(Candidate::from_moves(vec![Direction::PosX, Direction::NegX]).is_none());
    }

    #[test]
    fn random_and_mutated_candidates_remain_self_avoiding() {
        let mut rng = StdRng::seed_from_u64(7);
        let candidate = Candidate::random(60, &mut rng);
        let hypermutated = candidate.hypermutate(20, &mut rng);
        let macromutated = candidate.hypermacromutate(&mut rng);

        assert_eq!(candidate.coordinates().len(), 60);
        assert_eq!(hypermutated.coordinates().len(), 60);
        assert_eq!(macromutated.coordinates().len(), 60);
        assert!(candidate.is_self_avoiding());
        assert!(hypermutated.is_self_avoiding());
        assert!(macromutated.is_self_avoiding());
    }
}
