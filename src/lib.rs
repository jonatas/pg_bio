use pgrx::prelude::*;
use serde::{Deserialize, Serialize};

pgrx::pg_module_magic!();

// =====================================================================
// 1. 3D SPATIAL INDEXING (Z-ORDER)
// =====================================================================

#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct ResidueCoord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub name: String, 
}

/// SQL constructor to safely create a ResidueCoord without relying on internal CBOR string parsing
#[pg_extern(immutable, parallel_safe)]
pub fn create_residue_coord(x: f64, y: f64, z: f64, name: &str) -> ResidueCoord {
    ResidueCoord { x, y, z, name: name.to_string() }
}

#[inline]
fn split_by_3(x: u64) -> u64 {
    let mut res = 0;
    for i in 0..21 {
        res |= ((x >> i) & 1) << (3 * i);
    }
    res
}

#[pg_extern(immutable, parallel_safe)]
pub fn z_order_encode(x: f64, y: f64, z: f64) -> i64 {
    let map = |v: f64| -> u64 {
        let v = v.clamp(-500.0, 500.0) + 500.0;
        let scaled = (v * 1000.0) as u64; 
        scaled & 0x1FFFFF 
    };
    let xx = split_by_3(map(x));
    let yy = split_by_3(map(y));
    let zz = split_by_3(map(z));
    (xx | (yy << 1) | (zz << 2)) as i64
}

#[pg_extern(immutable, parallel_safe)]
pub fn residue_z_index(coord: ResidueCoord) -> i64 {
    z_order_encode(coord.x, coord.y, coord.z)
}

#[pg_extern(immutable, parallel_safe)]
pub fn distance_angstroms(a: ResidueCoord, b: ResidueCoord) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// =====================================================================
// 2. HIGH-DIMENSIONAL VECTOR EMBEDDINGS (ESM Computational models)
// =====================================================================

#[pg_extern(immutable, parallel_safe)]
pub fn embedding_cosine_distance(a: Vec<f32>, b: Vec<f32>) -> f64 {
    if a.len() != b.len() || a.is_empty() { return 2.0; }
    
    let mut dot_product: f64 = 0.0;
    let mut norm_a: f64 = 0.0;
    let mut norm_b: f64 = 0.0;

    for i in 0..a.len() {
        let val_a = a[i] as f64;
        let val_b = b[i] as f64;
        dot_product += val_a * val_b;
        norm_a += val_a * val_a;
        norm_b += val_b * val_b;
    }

    if norm_a == 0.0 || norm_b == 0.0 { return 2.0; }
    let similarity = dot_product / (norm_a.sqrt() * norm_b.sqrt());
    1.0 - similarity
}

#[pg_extern(immutable, parallel_safe)]
pub fn get_esm_embedding(sequence: &str) -> Vec<f32> {
    let mut vec = vec![0.0f32; 1280];
    if sequence.is_empty() { return vec; }
    
    let bytes = sequence.as_bytes();
    let k = 3; // Use tri-peptides
    if bytes.len() < k {
        for &b in bytes {
            let idx = (b as usize) % 1280;
            vec[idx] += 1.0;
        }
    } else {
        for window in bytes.windows(k) {
            let hash = (window[0] as usize).wrapping_mul(73)
                .wrapping_add((window[1] as usize).wrapping_mul(179))
                .wrapping_add((window[2] as usize).wrapping_mul(283));
            let idx = hash % 1280;
            vec[idx] += 1.0;
        }
    }
    
    // L2 Normalization
    let mut sum_sq = 0.0;
    for &v in vec.iter() {
        sum_sq += v * v;
    }
    if sum_sq > 0.0 {
        let norm = sum_sq.sqrt();
        for v in vec.iter_mut() {
            *v /= norm;
        }
    }
    vec
}

// =====================================================================
// 3. SPARSE ATTENTION MATRICES (Custom Memory-Optimized Storage)
// =====================================================================

/// Represents a single non-zero interaction in an Attention Map.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AttentionEntry {
    pub source_residue: i32, // Amino acid index initiating attention
    pub target_residue: i32, // Amino acid index being attended to
    pub weight: f32,         // Strength of the interaction (0.0 to 1.0)
}

/// A Postgres type that compresses an N x N attention matrix by only storing
/// the highly correlated (non-zero) weights, similar to a Compressed Sparse Row (CSR) format.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct SparseAttentionMap {
    pub sequence_length: i32,
    pub entries: Vec<AttentionEntry>,
}

/// SQL constructor to build a SparseAttentionMap from arrays
#[pg_extern(immutable, parallel_safe)]
pub fn create_sparse_map(length: i32, sources: Vec<i32>, targets: Vec<i32>, weights: Vec<f32>) -> SparseAttentionMap {
    let mut entries = Vec::new();
    let count = sources.len().min(targets.len()).min(weights.len());
    for i in 0..count {
        entries.push(AttentionEntry {
            source_residue: sources[i],
            target_residue: targets[i],
            weight: weights[i],
        });
    }
    SparseAttentionMap { sequence_length: length, entries }
}

/// A highly optimized native function that traverses the compressed attention map
/// to find which amino acids are physically/evolutionarily interacting with a specific target.
#[pg_extern(immutable, parallel_safe)]
pub fn get_top_interacting_residues(map: SparseAttentionMap, target: i32, limit: i32) -> Vec<i32> {
    let mut interactions: Vec<&AttentionEntry> = map.entries
        .iter()
        .filter(|e| e.target_residue == target || e.source_residue == target)
        .collect();

    // Sort by weight descending
    interactions.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));

    interactions.into_iter()
        .take(limit as usize)
        .map(|e| if e.target_residue == target { e.source_residue } else { e.target_residue })
        .collect()
}

// =====================================================================
// 4. BASIC CHEMISTRY & SEQUENCE MATH
// =====================================================================

/// Computes the molecular weight of an amino acid sequence (in Daltons).
#[pg_extern(immutable, parallel_safe)]
pub fn molecular_weight(sequence: &str) -> f64 {
    let mut weight = 18.015; // H2O (N-term H and C-term OH)
    for c in sequence.chars() {
        weight += match c.to_ascii_uppercase() {
            'A' => 71.079,
            'R' => 156.188,
            'N' => 114.104,
            'D' => 115.089,
            'C' => 103.145,
            'E' => 129.116,
            'Q' => 128.131,
            'G' => 57.052,
            'H' => 137.141,
            'I' => 113.160,
            'L' => 113.160,
            'K' => 128.174,
            'M' => 131.199,
            'F' => 147.177,
            'P' => 97.117,
            'S' => 87.078,
            'T' => 101.105,
            'W' => 186.213,
            'Y' => 163.176,
            'V' => 99.133,
            _ => 0.0, // Ignore unknown characters
        };
    }
    weight
}

/// Calculates the center of mass for a collection of atoms
#[pg_extern(immutable, parallel_safe)]
pub fn center_of_mass(atoms: Vec<ResidueCoord>) -> ResidueCoord {
    if atoms.is_empty() {
        return ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "COM".to_string() };
    }
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    let count = atoms.len() as f64;

    for atom in &atoms {
        sum_x += atom.x;
        sum_y += atom.y;
        sum_z += atom.z;
    }

    ResidueCoord {
        x: sum_x / count,
        y: sum_y / count,
        z: sum_z / count,
        name: "COM".to_string(),
    }
}

/// Computes the dihedral angle (in degrees) formed by four atoms (e.g. for Ramachandran Phi/Psi angles).
#[pg_extern(immutable, parallel_safe)]
pub fn dihedral_angle(a: ResidueCoord, b: ResidueCoord, c: ResidueCoord, d: ResidueCoord) -> f64 {
    let b1 = (b.x - a.x, b.y - a.y, b.z - a.z);
    let b2 = (c.x - b.x, c.y - b.y, c.z - b.z);
    let b3 = (d.x - c.x, d.y - c.y, d.z - c.z);

    let cross = |v1: (f64, f64, f64), v2: (f64, f64, f64)| -> (f64, f64, f64) {
        (
            v1.1 * v2.2 - v1.2 * v2.1,
            v1.2 * v2.0 - v1.0 * v2.2,
            v1.0 * v2.1 - v1.1 * v2.0,
        )
    };

    let dot = |v1: (f64, f64, f64), v2: (f64, f64, f64)| -> f64 {
        v1.0 * v2.0 + v1.1 * v2.1 + v1.2 * v2.2
    };

    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);

    let b2_mag = (b2.0 * b2.0 + b2.1 * b2.1 + b2.2 * b2.2).sqrt();
    let x = dot(n1, n2);
    let y = dot(b1, n2) * b2_mag;
    
    y.atan2(x).to_degrees()
}

/// Checks if two atoms are clashing sterically (closer than sum of Van der Waals radii minus some tolerance)
#[pg_extern(immutable, parallel_safe)]
pub fn is_steric_clash(a: ResidueCoord, b: ResidueCoord, vdw_radius_sum: f64) -> bool {
    // Avoid checking an atom against itself
    if a.x == b.x && a.y == b.y && a.z == b.z {
        return false;
    }
    crate::distance_angstroms(a, b) < vdw_radius_sum
}

/// A highly simplified native regex builder for biological PROSITE patterns.
/// E.g. [ST]-x(2)-[RK] -> matches Sequences natively in the database.
#[pg_extern(immutable, parallel_safe)]
pub fn prosite_match(sequence: &str, pattern: &str) -> bool {
    let mut regex_pattern = String::new();
    let mut chars = pattern.chars().peekable();
    
    while let Some(c) = chars.next() {
        match c {
            '-' => continue,
            'x' | 'X' => {
                regex_pattern.push('.');
            },
            '(' => {
                regex_pattern.push('{');
                while let Some(&next_c) = chars.peek() {
                    if next_c == ')' {
                        chars.next();
                        regex_pattern.push('}');
                        break;
                    }
                    regex_pattern.push(chars.next().unwrap());
                }
            },
            '{' => {
                regex_pattern.push_str("[^");
                while let Some(&next_c) = chars.peek() {
                    if next_c == '}' {
                        chars.next();
                        regex_pattern.push(']');
                        break;
                    }
                    regex_pattern.push(chars.next().unwrap());
                }
            },
            _ => regex_pattern.push(c)
        }
    }
    
    if let Ok(re) = regex::Regex::new(&regex_pattern) {
        re.is_match(sequence)
    } else {
        false
    }
}

// =====================================================================
// 5. SMALL MOLECULE CHEMISTRY (Cheminformatics: Hückel 4n+2, WL/ECFP4, VF2, Salt/Tautomer Normalization)
// =====================================================================

pub const BOND_SINGLE: u8 = 1;
pub const BOND_DOUBLE: u8 = 2;
pub const BOND_TRIPLE: u8 = 3;
pub const BOND_AROMATIC: u8 = 4;

/// Compact atom node inside a normalized `Molecule` graph.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct PgBioAtom {
    pub atomic_num: u8,
    pub formal_charge: i8,
    pub h_count: u8,
    pub is_aromatic: bool,
    pub in_ring: bool,
}

/// Compact undirected bond edge inside a normalized `Molecule` graph.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct PgBioBond {
    pub src_idx: u16,
    pub dst_idx: u16,
    pub bond_type: u8,
    pub in_ring: bool,
}

/// Represents a small molecule parsed from a SMILES string and normalized via:
/// 1. Connected-component salt stripping & prototropic tautomer normalization
/// 2. Hückel 4n+2 cycle aromaticity perception
/// 3. Weisfeiler-Lehman (WL) / Morgan ECFP4 64-bit canonical hashing & 1024-bit bitset indexing
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct Molecule {
    pub smiles: String,
    pub is_valid: bool,
    #[serde(default)]
    pub wl_hash: i64,
    #[serde(default)]
    pub ecfp_bits: Vec<u64>,
    #[serde(default)]
    pub atoms: Vec<PgBioAtom>,
    #[serde(default)]
    pub bonds: Vec<PgBioBond>,
}

#[inline]
fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}

fn symbol_to_atomic_num(sym: &str) -> Option<(u8, bool)> {
    match sym {
        "*" => Some((0, false)),
        "H" => Some((1, false)),
        "He" => Some((2, false)),
        "Li" => Some((3, false)),
        "Be" => Some((4, false)),
        "B" => Some((5, false)),
        "b" => Some((5, true)),
        "C" => Some((6, false)),
        "c" => Some((6, true)),
        "N" => Some((7, false)),
        "n" => Some((7, true)),
        "O" => Some((8, false)),
        "o" => Some((8, true)),
        "F" => Some((9, false)),
        "Na" => Some((11, false)),
        "Mg" => Some((12, false)),
        "Al" => Some((13, false)),
        "Si" => Some((14, false)),
        "P" => Some((15, false)),
        "p" => Some((15, true)),
        "S" => Some((16, false)),
        "s" => Some((16, true)),
        "Cl" => Some((17, false)),
        "K" => Some((19, false)),
        "Ca" => Some((20, false)),
        "Ti" => Some((22, false)),
        "V" => Some((23, false)),
        "Cr" => Some((24, false)),
        "Mn" => Some((25, false)),
        "Fe" => Some((26, false)),
        "Co" => Some((27, false)),
        "Ni" => Some((28, false)),
        "Cu" => Some((29, false)),
        "Zn" => Some((30, false)),
        "Ga" => Some((31, false)),
        "Ge" => Some((32, false)),
        "As" => Some((33, false)),
        "as" => Some((33, true)),
        "Se" => Some((34, false)),
        "se" => Some((34, true)),
        "Br" => Some((35, false)),
        "Ag" => Some((47, false)),
        "Cd" => Some((48, false)),
        "Sn" => Some((50, false)),
        "Sb" => Some((51, false)),
        "Te" => Some((52, false)),
        "te" => Some((52, true)),
        "I" => Some((53, false)),
        "Ba" => Some((56, false)),
        "Pt" => Some((78, false)),
        "Au" => Some((79, false)),
        "Hg" => Some((80, false)),
        "Pb" => Some((82, false)),
        "Bi" => Some((83, false)),
        _ => None,
    }
}

/// Parses a bracketed SMILES atom specification `[...]` (excluding `[` and `]`).
fn parse_bracket_atom(inner: &[u8]) -> Option<(u8, i8, u8, bool)> {
    let n = inner.len();
    if n == 0 {
        return None;
    }
    let mut i = 0;

    // 1. Optional isotope digits
    while i < n && inner[i].is_ascii_digit() {
        i += 1;
    }
    if i >= n {
        return None;
    }

    // 2. Element symbol (try 2-letter first, then 1-letter)
    let (atomic_num, is_aromatic) = if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
        let two = std::str::from_utf8(&inner[i..i + 2]).ok()?;
        if let Some(res) = symbol_to_atomic_num(two) {
            i += 2;
            res
        } else {
            let one = std::str::from_utf8(&inner[i..i + 1]).ok()?;
            let res = symbol_to_atomic_num(one)?;
            i += 1;
            res
        }
    } else {
        let one = std::str::from_utf8(&inner[i..i + 1]).ok()?;
        let res = symbol_to_atomic_num(one)?;
        i += 1;
        res
    };

    // 3. Optional chirality (@, @@, @TH1, etc.)
    if i < n && inner[i] == b'@' {
        i += 1;
        if i < n && inner[i] == b'@' {
            i += 1;
        } else {
            while i < n && (inner[i].is_ascii_uppercase() || inner[i].is_ascii_digit()) && inner[i] != b'H' {
                i += 1;
            }
        }
    }

    // 4. Optional explicit hydrogen count (H, H1, H2, ...)
    let mut h_count: u8 = 0;
    if i < n && inner[i] == b'H' {
        i += 1;
        if i < n && inner[i].is_ascii_digit() {
            let mut hc: u8 = 0;
            while i < n && inner[i].is_ascii_digit() {
                hc = hc.saturating_mul(10).saturating_add(inner[i] - b'0');
                i += 1;
            }
            h_count = hc;
        } else {
            h_count = 1;
        }
    }

    // 5. Optional formal charge (+, ++, +2, -, --, -2)
    let mut formal_charge: i8 = 0;
    if i < n && (inner[i] == b'+' || inner[i] == b'-') {
        let sign: i8 = if inner[i] == b'+' { 1 } else { -1 };
        let ch = inner[i];
        i += 1;
        if i < n && inner[i].is_ascii_digit() {
            let mut mag: i8 = 0;
            while i < n && inner[i].is_ascii_digit() {
                mag = mag.saturating_mul(10).saturating_add((inner[i] - b'0') as i8);
                i += 1;
            }
            formal_charge = sign * mag;
        } else {
            let mut count: i8 = 1;
            while i < n && inner[i] == ch {
                count += 1;
                i += 1;
            }
            formal_charge = sign * count;
        }
    }

    // 6. Optional atom class (:1, :2, ...)
    if i < n && inner[i] == b':' {
        i += 1;
        while i < n && inner[i].is_ascii_digit() {
            i += 1;
        }
    }

    if i != n {
        return None;
    }

    Some((atomic_num, formal_charge, h_count, is_aromatic))
}

/// Lexes and parses a raw SMILES string into `(atoms, bonds)` before salt/tautomer/aromaticity normalization.
fn pgbio_lex_smiles(smiles: &str) -> Option<(Vec<PgBioAtom>, Vec<PgBioBond>)> {
    let s = smiles.trim();
    if s.is_empty() {
        return None;
    }
    let bytes = s.as_bytes();
    let n = bytes.len();

    let mut atoms: Vec<PgBioAtom> = Vec::new();
    let mut is_bracket_atom: Vec<bool> = Vec::new();
    let mut bonds: Vec<PgBioBond> = Vec::new();

    let mut stack: Vec<usize> = Vec::new();
    let mut current_atom: Option<usize> = None;
    let mut pending_bond: Option<u8> = None;
    let mut just_saw_dot = false;
    // Ring closures 0..100: stores (opening_atom_idx, optional_explicit_bond_type)
    let mut rings: [Option<(usize, Option<u8>)>; 100] = [None; 100];

    let mut i = 0;
    while i < n {
        let c = bytes[i];
        match c {
            b'(' => {
                if pending_bond.is_some() || just_saw_dot {
                    return None;
                }
                let curr = current_atom?;
                stack.push(curr);
                i += 1;
            }
            b')' => {
                if pending_bond.is_some() || just_saw_dot {
                    return None;
                }
                current_atom = Some(stack.pop()?);
                i += 1;
            }
            b'.' => {
                if pending_bond.is_some() || current_atom.is_none() || just_saw_dot {
                    return None;
                }
                current_atom = None;
                just_saw_dot = true;
                i += 1;
            }
            b'-' | b'=' | b'#' | b':' | b'/' | b'\\' => {
                if pending_bond.is_some() || current_atom.is_none() {
                    return None;
                }
                pending_bond = Some(match c {
                    b'-' | b'/' | b'\\' => BOND_SINGLE,
                    b'=' => BOND_DOUBLE,
                    b'#' => BOND_TRIPLE,
                    b':' => BOND_AROMATIC,
                    _ => BOND_SINGLE,
                });
                i += 1;
            }
            b'0'..=b'9' | b'%' => {
                let curr = current_atom?;
                let ring_num: usize = if c == b'%' {
                    if i + 2 >= n || !bytes[i + 1].is_ascii_digit() || !bytes[i + 2].is_ascii_digit() {
                        return None;
                    }
                    let num = ((bytes[i + 1] - b'0') as usize) * 10 + ((bytes[i + 2] - b'0') as usize);
                    i += 3;
                    num
                } else {
                    i += 1;
                    (c - b'0') as usize
                };

                let explicit_b = pending_bond.take();
                if let Some((open_atom, open_bond)) = rings[ring_num].take() {
                    if open_atom == curr {
                        return None; // Self-loop invalid
                    }
                    // Determine ring closure bond type
                    let btype = match (open_bond, explicit_b) {
                        (Some(b1), Some(b2)) => {
                            if b1 != b2 {
                                return None;
                            }
                            b1
                        }
                        (Some(b1), None) => b1,
                        (None, Some(b2)) => b2,
                        (None, None) => {
                            if atoms[open_atom].is_aromatic && atoms[curr].is_aromatic {
                                BOND_AROMATIC
                            } else {
                                BOND_SINGLE
                            }
                        }
                    };
                    // Avoid duplicate edge
                    if bonds.iter().any(|b| {
                        (b.src_idx == open_atom as u16 && b.dst_idx == curr as u16)
                            || (b.src_idx == curr as u16 && b.dst_idx == open_atom as u16)
                    }) {
                        return None;
                    }
                    bonds.push(PgBioBond {
                        src_idx: open_atom as u16,
                        dst_idx: curr as u16,
                        bond_type: btype,
                        in_ring: true,
                    });
                } else {
                    rings[ring_num] = Some((curr, explicit_b));
                }
            }
            b'[' => {
                let start = i + 1;
                let mut end = start;
                while end < n && bytes[end] != b']' {
                    end += 1;
                }
                if end >= n {
                    return None;
                }
                let (atomic_num, formal_charge, h_count, is_aromatic) =
                    parse_bracket_atom(&bytes[start..end])?;
                let new_idx = atoms.len();
                atoms.push(PgBioAtom {
                    atomic_num,
                    formal_charge,
                    h_count,
                    is_aromatic,
                    in_ring: false,
                });
                is_bracket_atom.push(true);

                if let Some(prev) = current_atom {
                    let btype = pending_bond.take().unwrap_or_else(|| {
                        if atoms[prev].is_aromatic && is_aromatic {
                            BOND_AROMATIC
                        } else {
                            BOND_SINGLE
                        }
                    });
                    bonds.push(PgBioBond {
                        src_idx: prev as u16,
                        dst_idx: new_idx as u16,
                        bond_type: btype,
                        in_ring: false,
                    });
                } else if pending_bond.is_some() {
                    return None;
                }
                current_atom = Some(new_idx);
                just_saw_dot = false;
                i = end + 1;
            }
            _ => {
                // Unbracketed organic subset atom: Cl, Br, B, C, N, O, P, S, F, I, *, b, c, n, o, p, s
                let (atomic_num, is_aromatic, step) = if i + 1 < n
                    && ((c == b'C' && bytes[i + 1] == b'l') || (c == b'B' && bytes[i + 1] == b'r'))
                {
                    if c == b'C' {
                        (17u8, false, 2usize)
                    } else {
                        (35u8, false, 2usize)
                    }
                } else {
                    match c {
                        b'*' => (0u8, false, 1usize),
                        b'B' => (5u8, false, 1usize),
                        b'b' => (5u8, true, 1usize),
                        b'C' => (6u8, false, 1usize),
                        b'c' => (6u8, true, 1usize),
                        b'N' => (7u8, false, 1usize),
                        b'n' => (7u8, true, 1usize),
                        b'O' => (8u8, false, 1usize),
                        b'o' => (8u8, true, 1usize),
                        b'F' => (9u8, false, 1usize),
                        b'P' => (15u8, false, 1usize),
                        b'p' => (15u8, true, 1usize),
                        b'S' => (16u8, false, 1usize),
                        b's' => (16u8, true, 1usize),
                        b'I' => (53u8, false, 1usize),
                        _ => return None,
                    }
                };

                let new_idx = atoms.len();
                atoms.push(PgBioAtom {
                    atomic_num,
                    formal_charge: 0,
                    h_count: 0,
                    is_aromatic,
                    in_ring: false,
                });
                is_bracket_atom.push(false);

                if let Some(prev) = current_atom {
                    let btype = pending_bond.take().unwrap_or_else(|| {
                        if atoms[prev].is_aromatic && is_aromatic {
                            BOND_AROMATIC
                        } else {
                            BOND_SINGLE
                        }
                    });
                    bonds.push(PgBioBond {
                        src_idx: prev as u16,
                        dst_idx: new_idx as u16,
                        bond_type: btype,
                        in_ring: false,
                    });
                } else if pending_bond.is_some() {
                    return None;
                }
                current_atom = Some(new_idx);
                just_saw_dot = false;
                i += step;
            }
        }
    }

    if !stack.is_empty() || pending_bond.is_some() || just_saw_dot || atoms.is_empty() {
        return None;
    }
    if rings.iter().any(|r| r.is_some()) {
        return None;
    }

    // Fold explicit [H] atoms connected to heavy atoms into their neighbor's h_count
    let mut remove_atom = vec![false; atoms.len()];
    for b in &bonds {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        if atoms[u].atomic_num == 1 && atoms[u].formal_charge == 0 && atoms[v].atomic_num > 1 {
            atoms[v].h_count = atoms[v].h_count.saturating_add(1);
            remove_atom[u] = true;
        } else if atoms[v].atomic_num == 1 && atoms[v].formal_charge == 0 && atoms[u].atomic_num > 1 {
            atoms[u].h_count = atoms[u].h_count.saturating_add(1);
            remove_atom[v] = true;
        }
    }
    if remove_atom.iter().any(|&r| r) && remove_atom.iter().any(|&r| !r) {
        let mut remap = vec![0u16; atoms.len()];
        let mut new_atoms = Vec::new();
        let mut new_bracket = Vec::new();
        for (idx, a) in atoms.into_iter().enumerate() {
            if !remove_atom[idx] {
                remap[idx] = new_atoms.len() as u16;
                new_atoms.push(a);
                new_bracket.push(is_bracket_atom[idx]);
            }
        }
        let mut new_bonds = Vec::new();
        for b in bonds {
            if !remove_atom[b.src_idx as usize] && !remove_atom[b.dst_idx as usize] {
                new_bonds.push(PgBioBond {
                    src_idx: remap[b.src_idx as usize],
                    dst_idx: remap[b.dst_idx as usize],
                    bond_type: b.bond_type,
                    in_ring: b.in_ring,
                });
            }
        }
        atoms = new_atoms;
        is_bracket_atom = new_bracket;
        bonds = new_bonds;
    }

    // Compute implicit hydrogens for all unbracketed organic atoms
    for idx in 0..atoms.len() {
        if is_bracket_atom[idx] {
            continue;
        }
        let mut degree: u8 = 0;
        let mut explicit_mult: u8 = 0;
        for b in &bonds {
            if b.src_idx as usize == idx || b.dst_idx as usize == idx {
                degree = degree.saturating_add(1);
                if b.bond_type == BOND_DOUBLE {
                    explicit_mult = explicit_mult.saturating_add(1);
                } else if b.bond_type == BOND_TRIPLE {
                    explicit_mult = explicit_mult.saturating_add(2);
                }
            }
        }
        let a = &mut atoms[idx];
        if a.is_aromatic {
            // OpenSMILES aromatic implicit valence rule
            a.h_count = match a.atomic_num {
                6 => {
                    if explicit_mult > 0 {
                        4u8.saturating_sub(degree.saturating_add(explicit_mult))
                    } else {
                        3u8.saturating_sub(degree)
                    }
                }
                5 => 2u8.saturating_sub(degree),
                // Unbracketed aromatic n, o, s, p in SMILES have 0 implicit hydrogens (pyrrole NH is [nH])
                _ => 0,
            };
        } else {
            let bo_sum = degree.saturating_add(explicit_mult);
            let target_val = match a.atomic_num {
                5 => 3,
                6 => 4,
                7 | 15 => {
                    if bo_sum <= 3 {
                        3
                    } else if bo_sum <= 5 {
                        5
                    } else {
                        bo_sum
                    }
                }
                8 => 2,
                16 => {
                    if bo_sum <= 2 {
                        2
                    } else if bo_sum <= 4 {
                        4
                    } else if bo_sum <= 6 {
                        6
                    } else {
                        bo_sum
                    }
                }
                9 | 17 | 35 | 53 => 1,
                _ => bo_sum,
            };
            a.h_count = target_val.saturating_sub(bo_sum);
        }
    }

    Some((atoms, bonds))
}

/// Feature #3: Salt Stripping & Prototropic Tautomer Normalization.
/// Splits disconnected `.` components, keeps the largest heavy-atom parent component,
/// neutralizes counter-ion salt charges (e.g. `CC(=O)[O-].[Na+]` -> `CC(=O)O`),
/// and normalizes aliphatic 1,3-prototropic tautomers (`C=C-OH` -> `CH-C=O`, `N=C-OH` -> `NH-C=O`).
fn pgbio_strip_salts_and_normalize_tautomers(
    atoms: &mut Vec<PgBioAtom>,
    bonds: &mut Vec<PgBioBond>,
) {
    let n = atoms.len();
    if n > 1 {
        let mut comp_id = vec![usize::MAX; n];
        let mut components: Vec<Vec<usize>> = Vec::new();

        for start in 0..n {
            if comp_id[start] != usize::MAX {
                continue;
            }
            let cid = components.len();
            let mut comp_nodes = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            comp_id[start] = cid;
            queue.push_back(start);

            while let Some(u) = queue.pop_front() {
                comp_nodes.push(u);
                for b in bonds.iter() {
                    let v = if b.src_idx as usize == u {
                        b.dst_idx as usize
                    } else if b.dst_idx as usize == u {
                        b.src_idx as usize
                    } else {
                        continue;
                    };
                    if comp_id[v] == usize::MAX {
                        comp_id[v] = cid;
                        queue.push_back(v);
                    }
                }
            }
            components.push(comp_nodes);
        }

        if components.len() > 1 {
            // Score each component by (heavy_atom_count, sum_of_atomic_numbers)
            let mut best_cid = 0;
            let mut best_score = (0usize, 0usize);
            for (cid, nodes) in components.iter().enumerate() {
                let heavy = nodes.iter().filter(|&&u| atoms[u].atomic_num > 1).count();
                let z_sum: usize = nodes.iter().map(|&u| atoms[u].atomic_num as usize).sum();
                let score = (heavy, z_sum);
                if score > best_score {
                    best_score = score;
                    best_cid = cid;
                }
            }

            // Sum net charge of stripped counter-ions
            let mut stripped_charge: i32 = 0;
            for (cid, nodes) in components.iter().enumerate() {
                if cid != best_cid {
                    for &u in nodes {
                        stripped_charge += atoms[u].formal_charge as i32;
                    }
                }
            }

            let keep_nodes = &components[best_cid];
            let mut remap = vec![u16::MAX; n];
            let mut new_atoms = Vec::with_capacity(keep_nodes.len());
            for &old_u in keep_nodes {
                remap[old_u] = new_atoms.len() as u16;
                new_atoms.push(atoms[old_u].clone());
            }
            let mut new_bonds = Vec::new();
            for b in bonds.iter() {
                let ru = remap[b.src_idx as usize];
                let rv = remap[b.dst_idx as usize];
                if ru != u16::MAX && rv != u16::MAX {
                    new_bonds.push(PgBioBond {
                        src_idx: ru,
                        dst_idx: rv,
                        bond_type: b.bond_type,
                        in_ring: b.in_ring,
                    });
                }
            }
            *atoms = new_atoms;
            *bonds = new_bonds;

            // If positively charged counter-ion(s) were stripped (e.g. [Na+], [K+], [Ca+2]),
            // protonate matching deprotonated oxygens/sulfurs ([O-] -> OH) in the parent molecule.
            if stripped_charge > 0 {
                let mut pos_rem = stripped_charge;
                for a in atoms.iter_mut() {
                    if pos_rem > 0 && (a.atomic_num == 8 || a.atomic_num == 16) && a.formal_charge == -1 {
                        a.formal_charge = 0;
                        a.h_count = a.h_count.saturating_add(1);
                        pos_rem -= 1;
                    }
                }
            } else if stripped_charge < 0 {
                // If halide/anion counter-ion was stripped (e.g. .[Cl-]), neutralize protonated amine [NH3+] -> NH2
                let mut neg_rem = -stripped_charge;
                for a in atoms.iter_mut() {
                    if neg_rem > 0 && a.atomic_num == 7 && a.formal_charge == 1 && a.h_count > 0 {
                        a.formal_charge = 0;
                        a.h_count -= 1;
                        neg_rem -= 1;
                    }
                }
            }
        }
    }

    // Enumerate simple cycles first so we know which bonds are in rings before tautomer shifts
    let cycles = pgbio_find_cycles(atoms.len(), bonds);
    for (_, cyc_bonds) in &cycles {
        for &b_idx in cyc_bonds {
            bonds[b_idx as usize].in_ring = true;
        }
    }

    // 1,3-Prototropic Tautomer Normalization:
    // Normalize aliphatic/non-ring enols and imidic acids: X=C-OH -> XH-C=O (where X is C or N and X=C is not in a cycle)
    for b_oh_idx in 0..bonds.len() {
        if bonds[b_oh_idx].bond_type != BOND_SINGLE || bonds[b_oh_idx].in_ring {
            continue;
        }
        let u = bonds[b_oh_idx].src_idx as usize;
        let v = bonds[b_oh_idx].dst_idx as usize;
        let (c_idx, o_idx) = if atoms[u].atomic_num == 6 && atoms[v].atomic_num == 8 {
            (u, v)
        } else if atoms[v].atomic_num == 6 && atoms[u].atomic_num == 8 {
            (v, u)
        } else {
            continue;
        };

        // O must be neutral -OH with degree 1
        if atoms[o_idx].formal_charge != 0 || atoms[o_idx].h_count != 1 {
            continue;
        }
        let o_deg = bonds
            .iter()
            .filter(|b| b.src_idx as usize == o_idx || b.dst_idx as usize == o_idx)
            .count();
        if o_deg != 1 {
            continue;
        }

        // Look for an non-ring double bond C=X (where X is neutral C or N) incident to c_idx
        let mut double_bond_idx: Option<(usize, usize)> = None;
        for (b_idx, b) in bonds.iter().enumerate() {
            if b.bond_type == BOND_DOUBLE && !b.in_ring {
                let x_idx = if b.src_idx as usize == c_idx {
                    b.dst_idx as usize
                } else if b.dst_idx as usize == c_idx {
                    b.src_idx as usize
                } else {
                    continue;
                };
                if (atoms[x_idx].atomic_num == 6 || atoms[x_idx].atomic_num == 7)
                    && atoms[x_idx].formal_charge == 0
                    && !atoms[x_idx].is_aromatic
                {
                    double_bond_idx = Some((b_idx, x_idx));
                    break;
                }
            }
        }

        if let Some((b_cx_idx, x_idx)) = double_bond_idx {
            atoms[o_idx].h_count = 0;
            atoms[x_idx].h_count = atoms[x_idx].h_count.saturating_add(1);
            bonds[b_oh_idx].bond_type = BOND_DOUBLE;
            bonds[b_cx_idx].bond_type = BOND_SINGLE;
        }
    }
}

/// Finds all simple cycles of length 3..=10 in the molecular graph.
/// Returns pairs of `(cycle_atom_indices, cycle_bond_indices)`.
fn pgbio_find_cycles(num_atoms: usize, bonds: &[PgBioBond]) -> Vec<(Vec<u16>, Vec<u16>)> {
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); num_atoms];
    for (b_idx, b) in bonds.iter().enumerate() {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        if u < num_atoms && v < num_atoms {
            adj[u].push((v, b_idx));
            adj[v].push((u, b_idx));
        }
    }

    let mut cycles: Vec<(Vec<u16>, Vec<u16>)> = Vec::new();
    let mut seen_edge_masks: Vec<u128> = Vec::new();

    // DFS from each start vertex `s` restricted to vertices >= `s` to find all simple cycles up to length 10
    for s in 0..num_atoms {
        let mut path_atoms: Vec<u16> = vec![s as u16];
        let mut path_bonds: Vec<u16> = Vec::new();
        let mut visited = vec![false; num_atoms];
        visited[s] = true;

        fn dfs_cycle(
            s: usize,
            u: usize,
            adj: &[Vec<(usize, usize)>],
            visited: &mut [bool],
            path_atoms: &mut Vec<u16>,
            path_bonds: &mut Vec<u16>,
            cycles: &mut Vec<(Vec<u16>, Vec<u16>)>,
            seen_edge_masks: &mut Vec<u128>,
        ) {
            for &(v, b_idx) in &adj[u] {
                if !path_bonds.is_empty() && path_bonds[path_bonds.len() - 1] == b_idx as u16 {
                    continue; // Do not immediately traverse back along the same edge
                }
                if v == s && path_atoms.len() >= 3 {
                    let mut cyc_bonds = path_bonds.clone();
                    cyc_bonds.push(b_idx as u16);
                    let mut mask: u128 = 0;
                    for &eb in &cyc_bonds {
                        if (eb as usize) < 128 {
                            mask |= 1u128 << (eb as usize);
                        }
                    }
                    if !seen_edge_masks.contains(&mask) {
                        seen_edge_masks.push(mask);
                        cycles.push((path_atoms.clone(), cyc_bonds));
                    }
                } else if v > s && !visited[v] && path_atoms.len() < 10 {
                    visited[v] = true;
                    path_atoms.push(v as u16);
                    path_bonds.push(b_idx as u16);
                    dfs_cycle(s, v, adj, visited, path_atoms, path_bonds, cycles, seen_edge_masks);
                    path_bonds.pop();
                    path_atoms.pop();
                    visited[v] = false;
                }
            }
        }

        dfs_cycle(
            s,
            s,
            &adj,
            &mut visited,
            &mut path_atoms,
            &mut path_bonds,
            &mut cycles,
            &mut seen_edge_masks,
        );
    }

    // Sort cycles by length ascending so smaller fundamental rings (5..7) are normalized first
    cycles.sort_by_key(|(ca, _)| ca.len());
    cycles
}

/// Return pi-electron contribution for atom `atom_idx` inside a candidate ring,
/// or -1 if the atom breaks sp2 conjugation (cannot be aromatic).
#[inline]
pub fn pgbio_atom_pi_electrons(
    atoms: &[PgBioAtom],
    bonds: &[PgBioBond],
    atom_idx: usize,
) -> i32 {
    let a = &atoms[atom_idx];
    let mut has_ring_double_or_arom = false;
    let mut has_exocyclic_electroneg_double = false;
    let mut has_any_double = false;
    let mut explicit_bonds: i32 = 0;

    for b in bonds {
        if b.src_idx as usize != atom_idx && b.dst_idx as usize != atom_idx {
            continue;
        }
        let other = if b.src_idx as usize == atom_idx {
            b.dst_idx as usize
        } else {
            b.src_idx as usize
        };

        if b.bond_type == BOND_TRIPLE {
            return -1; // sp-hybridized triple bonds break standard planar aromaticity
        }

        if b.bond_type == BOND_DOUBLE || b.bond_type == BOND_AROMATIC {
            if b.in_ring {
                has_ring_double_or_arom = true;
            } else if b.bond_type == BOND_DOUBLE
                && (atoms[other].atomic_num == 7
                    || atoms[other].atomic_num == 8
                    || atoms[other].atomic_num == 16)
            {
                has_exocyclic_electroneg_double = true; // e.g., C(=O) inside Caffeine/Uracil/Pyridone
            }
        }
        if b.bond_type == BOND_DOUBLE {
            has_any_double = true;
        }
        explicit_bonds += 1;
    }

    let sigma_degree = explicit_bonds + (a.h_count as i32);
    if sigma_degree >= 4 {
        return -1; // Saturated sp3 atom (like -CH2- in Cyclohexane/Cyclopentadiene) -> breaks conjugation!
    }

    match a.atomic_num {
        6 => {
            // Carbon with exocyclic =O / =N / =S or carbocation [C+] contributes 0 pi electrons but preserves cycle
            if (has_exocyclic_electroneg_double && a.formal_charge == 0)
                || (a.formal_charge == 1 && sigma_degree == 3 && !has_any_double)
            {
                return 0;
            }
            // Carbanion [C-] (e.g., cyclopentadienyl anion [cH-]1cccc1 / [CH-]1C=CC=C1) -> donates 2 pi electrons
            if a.formal_charge == -1 && sigma_degree == 3 && !has_any_double {
                return 2;
            }
            // Standard sp2 ring carbon with in-ring = or : bond
            if (a.is_aromatic || has_ring_double_or_arom) && sigma_degree == 3 && a.formal_charge == 0 {
                return 1;
            }
            -1
        }
        7 | 15 => {
            if a.formal_charge == 0 {
                // Pyrrole-like N (with 1 H or exocyclic substituent, 3 single/aromatic sigma bonds, no explicit = bond) -> donates lone pair (2 e-)
                if !has_any_double && sigma_degree == 3 {
                    return 2;
                }
                // Pyridine-like N (=N- in Kekule or lowercase n with 2 sigma bonds and 0 H) -> contributes 1 pi electron
                if (a.is_aromatic || has_ring_double_or_arom) && sigma_degree == 2 {
                    return 1;
                }
            } else if a.formal_charge == 1 && sigma_degree == 3 && (a.is_aromatic || has_ring_double_or_arom) {
                // Pyridinium [nH+] -> contributes 1 pi electron
                return 1;
            } else if a.formal_charge == -1 && sigma_degree == 2 && !has_any_double {
                // Pyrrolide anion [n-] -> donates 2 pi electrons
                return 2;
            }
            -1
        }
        8 | 16 | 34 => {
            // Neutral O, S, Se with 2 single/aromatic ring bonds -> donates lone pair (2 e-)
            if a.formal_charge == 0 && sigma_degree == 2 && !has_any_double {
                return 2;
            }
            // Pyrylium / thiopyrylium [o+] / [s+] -> contributes 1 pi electron
            if a.formal_charge == 1 && sigma_degree == 3 && (a.is_aromatic || has_ring_double_or_arom) {
                return 1;
            }
            -1
        }
        5 => {
            // Trigonal planar Boron with empty p-orbital -> 0 pi electrons
            if a.formal_charge == 0 && sigma_degree == 3 && !has_any_double {
                return 0;
            }
            -1
        }
        _ => -1,
    }
}

/// Normalize a detected cycle (array of vertex indices and edge indices) using Hückel's 4n+2 rule:
/// `(pi_sum >= 6) && ((pi_sum - 2) & 3) == 0`.
pub fn pgbio_normalize_cycle(
    atoms: &mut [PgBioAtom],
    bonds: &mut [PgBioBond],
    cycle_atoms: &[u16],
    cycle_bonds: &[u16],
) -> bool {
    let mut pi_sum: i32 = 0;

    for &atom_idx in cycle_atoms {
        let pi = pgbio_atom_pi_electrons(atoms, bonds, atom_idx as usize);
        if pi < 0 {
            return false; // Not conjugated -> leave bonds as explicit single/double
        }
        pi_sum += pi;
    }

    // Hückel 4n + 2 check: 6, 10, 14, 18...
    if pi_sum >= 6 && ((pi_sum - 2) & 3) == 0 {
        for &atom_idx in cycle_atoms {
            atoms[atom_idx as usize].is_aromatic = true;
        }
        for &bond_idx in cycle_bonds {
            bonds[bond_idx as usize].bond_type = BOND_AROMATIC;
        }
        return true;
    }
    false
}

/// Runs Hückel 4n+2 aromaticity perception across all cycles in the molecular graph.
fn pgbio_perceive_aromaticity(atoms: &mut [PgBioAtom], bonds: &mut [PgBioBond]) {
    let cycles = pgbio_find_cycles(atoms.len(), bonds);

    // Mark in_ring flags on all atoms and bonds belonging to any cycle
    for b in bonds.iter_mut() {
        b.in_ring = false;
    }
    for a in atoms.iter_mut() {
        a.in_ring = false;
    }
    for (cyc_atoms, cyc_bonds) in &cycles {
        for &u in cyc_atoms {
            atoms[u as usize].in_ring = true;
        }
        for &b_idx in cyc_bonds {
            bonds[b_idx as usize].in_ring = true;
        }
    }

    // Non-ring bonds between two aromatic atoms (e.g. biphenyl `c1ccccc1c2ccccc2`) are single bonds in SMILES,
    // whereas acyclic query fragments (where neither atom is in a ring, e.g. `cc`) retain BOND_AROMATIC.
    for b in bonds.iter_mut() {
        if !b.in_ring
            && b.bond_type == BOND_AROMATIC
            && (atoms[b.src_idx as usize].in_ring || atoms[b.dst_idx as usize].in_ring)
        {
            b.bond_type = BOND_SINGLE;
        }
    }

    if cycles.is_empty() {
        return;
    }

    let mut aromatic_atom_mask = vec![false; atoms.len()];
    let mut aromatic_bond_mask = vec![false; bonds.len()];

    // Pass 1 & 2: Normalize fundamental rings (e.g. 5..7) and fused perimeters (up to 10)
    for _pass in 0..2 {
        for (cyc_atoms, cyc_bonds) in &cycles {
            if pgbio_normalize_cycle(atoms, bonds, cyc_atoms, cyc_bonds) {
                for &u in cyc_atoms {
                    aromatic_atom_mask[u as usize] = true;
                }
                for &b_idx in cyc_bonds {
                    aromatic_bond_mask[b_idx as usize] = true;
                }
            }
        }
    }

    // Demote any ring atoms/bonds that did not belong to any valid 4n+2 Hückel cycle
    // (e.g. anti-aromatic cyclobutadiene `c1ccc1` or cyclooctatetraene `c1ccccccc1`)
    for (idx, a) in atoms.iter_mut().enumerate() {
        if a.in_ring && !aromatic_atom_mask[idx] {
            a.is_aromatic = false;
        }
    }
    for (idx, b) in bonds.iter_mut().enumerate() {
        if b.in_ring && !aromatic_bond_mask[idx] && b.bond_type == BOND_AROMATIC {
            b.bond_type = BOND_SINGLE;
        }
    }
}

/// Computes the canonical 64-bit Weisfeiler-Lehman (WL) graph hash and 1024-bit ECFP4 fingerprint
/// (`16 x u64` bitset) for a normalized molecular graph.
pub fn pgbio_compute_wl_ecfp(atoms: &[PgBioAtom], bonds: &[PgBioBond]) -> (i64, Vec<u64>) {
    let mut ecfp_bits = vec![0u64; 16];
    let n = atoms.len();
    if n == 0 {
        return (0, ecfp_bits);
    }

    // Precompute adjacency list: (neighbor_idx, bond_type)
    let mut adj: Vec<Vec<(usize, u8)>> = vec![Vec::new(); n];
    for b in bonds {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        if u < n && v < n {
            adj[u].push((v, b.bond_type));
            adj[v].push((u, b.bond_type));
        }
    }

    // For ring-automorphism invariance on heterocycles (e.g., Imidazole tautomer `C1=CN=CN1` vs `c1c[nH]cn1`),
    // aromatic ring nitrogens sharing an aromatic ring with both =N- and -NH- normalize tautomer-mobile H in the atom invariant.
    let aromatic_n_count = atoms
        .iter()
        .filter(|a| a.is_aromatic && a.in_ring && a.atomic_num == 7 && a.formal_charge == 0)
        .count();
    let aromatic_nh_count = atoms
        .iter()
        .filter(|a| {
            a.is_aromatic && a.in_ring && a.atomic_num == 7 && a.formal_charge == 0 && a.h_count == 1
        })
        .count();
    let has_tautomeric_arom_n = aromatic_n_count >= 2 && aromatic_nh_count >= 1 && aromatic_nh_count < aromatic_n_count;

    // Radius 0: Vertex-permutation-invariant initial atom invariants
    let mut hashes: Vec<u64> = (0..n)
        .map(|i| {
            let a = &atoms[i];
            let deg = adj[i].len() as u64;
            let effective_h = if has_tautomeric_arom_n
                && a.is_aromatic
                && a.in_ring
                && a.atomic_num == 7
                && a.formal_charge == 0
                && deg == 2
            {
                1u64
            } else {
                a.h_count as u64
            };
            let packed: u64 = ((a.atomic_num as u64) & 0xFF)
                | (((a.is_aromatic as u64) & 1) << 8)
                | ((effective_h & 0xF) << 9)
                | ((((a.formal_charge as i16 + 8) as u64) & 0xF) << 13)
                | ((deg & 0xF) << 17)
                | (((a.in_ring as u64) & 1) << 21);
            let h0 = fmix64(packed ^ 0x9e3779b97f4a7c15);
            let bit = (h0 as usize) & 1023;
            ecfp_bits[bit >> 6] |= 1u64 << (bit & 63);
            h0
        })
        .collect();

    let mut all_wl_features: Vec<u64> = hashes.iter().map(|&h| fmix64(h ^ 0x100)).collect();

    // Weisfeiler-Lehman / Morgan ECFP4 message passing (Radius 1..=2)
    for radius in 1u64..=2u64 {
        let mut next_hashes = vec![0u64; n];
        for i in 0..n {
            let mut nbr_pairs: Vec<(u8, u64)> = adj[i]
                .iter()
                .map(|&(nbr, btype)| (btype, hashes[nbr]))
                .collect();
            nbr_pairs.sort_unstable();

            let mut acc = fmix64(hashes[i] ^ (radius.wrapping_mul(0x517cc1b727220a95)));
            for (btype, nh) in nbr_pairs {
                let edge_token = fmix64(nh.wrapping_add((btype as u64).wrapping_mul(0x9e3779b97f4a7c15)));
                acc = fmix64(acc.wrapping_mul(1099511628211).wrapping_add(edge_token));
            }
            next_hashes[i] = acc;
            let bit = (acc as usize) & 1023;
            ecfp_bits[bit >> 6] |= 1u64 << (bit & 63);
            all_wl_features.push(fmix64(acc ^ (radius << 32)));
        }
        hashes = next_hashes;
    }

    // Sort all multi-radius WL features to produce a canonical, permutation-invariant 64-bit graph hash
    all_wl_features.sort_unstable();
    let mut graph_hash: u64 = fmix64((n as u64) ^ ((bonds.len() as u64) << 16));
    for feat in all_wl_features {
        graph_hash = fmix64(graph_hash.wrapping_mul(0x100000001b3).wrapping_add(feat));
    }

    (graph_hash as i64, ecfp_bits)
}

/// Feature #2: VF2 Subgraph Isomorphism state machine.
/// Determines whether `query` is isomorphic to a subgraph of `target`.
pub fn pgbio_vf2_is_subgraph(target: &Molecule, query: &Molecule) -> bool {
    if !target.is_valid || !query.is_valid {
        return false;
    }
    let n_q = query.atoms.len();
    let n_t = target.atoms.len();
    if n_q == 0 || n_q > n_t || query.bonds.len() > target.bonds.len() {
        return false;
    }

    // Fast O(1) exact match shortcut if atom & bond counts match and canonical WL hashes match
    if n_q == n_t && query.bonds.len() == target.bonds.len() && query.wl_hash == target.wl_hash {
        return true;
    }

    // Build adjacency tables of bond types (0 = no bond)
    let mut q_deg = vec![0usize; n_q];
    let mut q_adj: Vec<Vec<(usize, u8)>> = vec![Vec::new(); n_q];
    for b in &query.bonds {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        q_deg[u] += 1;
        q_deg[v] += 1;
        q_adj[u].push((v, b.bond_type));
        q_adj[v].push((u, b.bond_type));
    }

    let mut t_deg = vec![0usize; n_t];
    let mut t_adj: Vec<Vec<(usize, u8)>> = vec![Vec::new(); n_t];
    for b in &target.bonds {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        t_deg[u] += 1;
        t_deg[v] += 1;
        t_adj[u].push((v, b.bond_type));
        t_adj[v].push((u, b.bond_type));
    }

    // Order query vertices to maximize pruning: connected DFS order starting from highest-degree/rarest atom
    let mut q_order: Vec<usize> = Vec::with_capacity(n_q);
    let mut in_order = vec![false; n_q];
    while q_order.len() < n_q {
        // Pick next vertex connected to already-ordered set (or highest degree if starting a component)
        let mut best_u: Option<usize> = None;
        let mut best_priority = (0usize, 0usize, 0usize);
        for u in 0..n_q {
            if in_order[u] {
                continue;
            }
            let conn_to_mapped = q_adj[u].iter().filter(|&&(nbr, _)| in_order[nbr]).count();
            let priority = (conn_to_mapped, q_deg[u], query.atoms[u].atomic_num as usize);
            if best_u.is_none() || priority > best_priority {
                best_u = Some(u);
                best_priority = priority;
            }
        }
        let u = best_u.unwrap();
        in_order[u] = true;
        q_order.push(u);
    }

    let mut q_to_t: Vec<Option<usize>> = vec![None; n_q];
    let mut t_used: Vec<bool> = vec![false; n_t];

    fn vf2_backtrack(
        step: usize,
        q_order: &[usize],
        query: &Molecule,
        target: &Molecule,
        q_deg: &[usize],
        t_deg: &[usize],
        q_adj: &[Vec<(usize, u8)>],
        t_adj: &[Vec<(usize, u8)>],
        q_to_t: &mut [Option<usize>],
        t_used: &mut [bool],
    ) -> bool {
        if step == q_order.len() {
            return true;
        }
        let q_u = q_order[step];
        let qa = &query.atoms[q_u];

        for t_v in 0..target.atoms.len() {
            if t_used[t_v] {
                continue;
            }
            let ta = &target.atoms[t_v];

            // Node compatibility check: atomic number (0 = wildcard `*`), aromaticity, and degree bound
            if qa.atomic_num != 0 && qa.atomic_num != ta.atomic_num {
                continue;
            }
            if qa.is_aromatic != ta.is_aromatic {
                continue;
            }
            if t_deg[t_v] < q_deg[q_u] {
                continue;
            }
            if qa.in_ring && !ta.in_ring {
                continue;
            }

            // Edge consistency check against all already-mapped neighbors of q_u
            let mut edges_ok = true;
            for &(q_nbr, q_btype) in &q_adj[q_u] {
                if let Some(t_nbr) = q_to_t[q_nbr] {
                    let found_bond = t_adj[t_v]
                        .iter()
                        .any(|&(cand_nbr, t_btype)| cand_nbr == t_nbr && t_btype == q_btype);
                    if !found_bond {
                        edges_ok = false;
                        break;
                    }
                }
            }
            if !edges_ok {
                continue;
            }

            q_to_t[q_u] = Some(t_v);
            t_used[t_v] = true;
            if vf2_backtrack(
                step + 1,
                q_order,
                query,
                target,
                q_deg,
                t_deg,
                q_adj,
                t_adj,
                q_to_t,
                t_used,
            ) {
                return true;
            }
            q_to_t[q_u] = None;
            t_used[t_v] = false;
        }

        false
    }

    vf2_backtrack(
        0,
        &q_order,
        query,
        target,
        &q_deg,
        &t_deg,
        &q_adj,
        &t_adj,
        &mut q_to_t,
        &mut t_used,
    )
}

/// Parses, validates, strips salts, normalizes tautomers, perceives Hückel 4n+2 aromaticity,
/// and computes the Weisfeiler-Lehman hash & 1024-bit ECFP4 fingerprint for a SMILES string.
#[pg_extern(immutable, parallel_safe)]
pub fn parse_smiles(smiles_str: &str) -> Molecule {
    if let Some((mut atoms, mut bonds)) = pgbio_lex_smiles(smiles_str) {
        pgbio_strip_salts_and_normalize_tautomers(&mut atoms, &mut bonds);
        pgbio_perceive_aromaticity(&mut atoms, &mut bonds);
        let (wl_hash, ecfp_bits) = pgbio_compute_wl_ecfp(&atoms, &bonds);
        Molecule {
            smiles: smiles_str.to_string(),
            is_valid: true,
            wl_hash,
            ecfp_bits,
            atoms,
            bonds,
        }
    } else {
        Molecule {
            smiles: smiles_str.to_string(),
            is_valid: false,
            wl_hash: 0,
            ecfp_bits: Vec::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }
}

/// Returns the canonical 64-bit Weisfeiler-Lehman graph hash after Hückel 4n+2 aromaticity perception
/// and salt/tautomer normalization (e.g., `molecule_wl_hash('C1=CC=CC=C1') == molecule_wl_hash('c1ccccc1')`).
#[pg_extern(immutable, parallel_safe)]
pub fn molecule_wl_hash(smiles: &str) -> i64 {
    parse_smiles(smiles).wl_hash
}

/// Checks whether two SMILES strings represent the exact same normalized molecular graph
/// (e.g. Kekulé vs. aromatic `C1=CC=CC=C1` vs. `c1ccccc1`, or salted `CC(=O)O.[Na+]` vs. `CC(=O)O`).
#[pg_extern(immutable, parallel_safe)]
pub fn molecules_are_identical(smiles_a: &str, smiles_b: &str) -> bool {
    let mol_a = parse_smiles(smiles_a);
    let mol_b = parse_smiles(smiles_b);
    if !mol_a.is_valid || !mol_b.is_valid {
        return false;
    }
    mol_a.wl_hash == mol_b.wl_hash
        && mol_a.atoms.len() == mol_b.atoms.len()
        && mol_a.bonds.len() == mol_b.bonds.len()
        && pgbio_vf2_is_subgraph(&mol_a, &mol_b)
}

/// Performs an exact VF2 Subgraph Isomorphism check (with Hückel 4n+2 aromaticity normalization)
/// to test whether `target` contains `substructure` as a molecular subgraph.
#[pg_extern(immutable, parallel_safe)]
pub fn smiles_contains(target: Molecule, substructure: &str) -> bool {
    if !target.is_valid {
        return false;
    }
    let target_norm = if target.atoms.is_empty() && !target.smiles.is_empty() {
        parse_smiles(&target.smiles)
    } else {
        target
    };
    let query_mol = parse_smiles(substructure);
    if !query_mol.is_valid {
        return false;
    }
    pgbio_vf2_is_subgraph(&target_norm, &query_mol)
}

/// Subgraph containment operator procedure (`target @> query`) combining bitset screening with VF2 subgraph isomorphism recheck.
#[pg_extern(immutable, parallel_safe)]
pub fn molecule_contains_molecule(target: Molecule, query: Molecule) -> bool {
    if !target.is_valid || !query.is_valid {
        return false;
    }
    let target_norm = if target.atoms.is_empty() && !target.smiles.is_empty() {
        parse_smiles(&target.smiles)
    } else {
        target
    };
    let query_norm = if query.atoms.is_empty() && !query.smiles.is_empty() {
        parse_smiles(&query.smiles)
    } else {
        query
    };
    pgbio_vf2_is_subgraph(&target_norm, &query_norm)
}

/// Hardware popcount (`count_ones`) Tanimoto similarity directly on two `Molecule` 1024-bit ECFP4 bitsets.
#[pg_extern(immutable, parallel_safe)]
pub fn molecule_tanimoto(a: Molecule, b: Molecule) -> f64 {
    if !a.is_valid || !b.is_valid || a.ecfp_bits.len() != 16 || b.ecfp_bits.len() != 16 {
        return 0.0;
    }
    let mut inter: u32 = 0;
    let mut union: u32 = 0;
    for i in 0..16 {
        inter += (a.ecfp_bits[i] & b.ecfp_bits[i]).count_ones();
        union += (a.ecfp_bits[i] | b.ecfp_bits[i]).count_ones();
    }
    if union == 0 {
        1.0
    } else {
        (inter as f64) / (union as f64)
    }
}

/// Native Morgan (ECFP4) 1024-bit Fingerprint generation in Rust with Hückel 4n+2 aromaticity normalization.
#[pg_extern(immutable, parallel_safe)]
pub fn smiles_to_fingerprint(smiles: &str) -> Vec<f32> {
    let mol = parse_smiles(smiles);
    let mut fp = vec![0.0f32; 1024];
    if !mol.is_valid || mol.ecfp_bits.len() != 16 {
        return fp;
    }
    for word_idx in 0..16 {
        let word = mol.ecfp_bits[word_idx];
        if word == 0 {
            continue;
        }
        for bit in 0..64 {
            if (word & (1u64 << bit)) != 0 {
                fp[(word_idx << 6) | bit] = 1.0;
            }
        }
    }
    fp
}

/// Tanimoto similarity function for chemical fingerprints
#[pg_extern(immutable, parallel_safe)]
pub fn tanimoto_similarity(a: Vec<f32>, b: Vec<f32>) -> f64 {
    if a.len() != b.len() || a.is_empty() { return 0.0; }
    
    let mut dot_product: f64 = 0.0;
    let mut norm_a: f64 = 0.0;
    let mut norm_b: f64 = 0.0;

    for i in 0..a.len() {
        let val_a = a[i] as f64;
        let val_b = b[i] as f64;
        dot_product += val_a * val_b;
        norm_a += val_a * val_a;
        norm_b += val_b * val_b;
    }

    if norm_a == 0.0 && norm_b == 0.0 { return 1.0; }
    let denominator = norm_a + norm_b - dot_product;
    if denominator == 0.0 { return 0.0; }
    
    dot_product / denominator
}

pgrx::extension_sql!(
    r#"
-- 5. Create Tanimoto similarity and VF2 substructure operators
CREATE OPERATOR % (
    LEFTARG = real[],
    RIGHTARG = real[],
    PROCEDURE = tanimoto_similarity,
    COMMUTATOR = %
);

CREATE OPERATOR % (
    LEFTARG = Molecule,
    RIGHTARG = Molecule,
    PROCEDURE = molecule_tanimoto,
    COMMUTATOR = %
);

CREATE OPERATOR @> (
    LEFTARG = Molecule,
    RIGHTARG = Molecule,
    PROCEDURE = molecule_contains_molecule
);
"#,
    name = "pg_bio_tanimoto_operator",
    requires = [tanimoto_similarity, molecule_tanimoto, molecule_contains_molecule]
);

// =====================================================================
// 6. SEQUENCE COMPLEXITY
// =====================================================================

/// Computes the Shannon entropy of a sequence (0.0 to ~4.32 for 20 amino acids)
/// Low complexity sequences (like poly-G) will have very low entropy.
#[pg_extern(immutable, parallel_safe)]
pub fn sequence_entropy(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    
    let mut counts = [0; 256];
    for b in sequence.bytes() {
        counts[b as usize] += 1;
    }
    
    let len = sequence.len() as f64;
    let mut entropy = 0.0;
    
    for &count in counts.iter() {
        if count > 0 {
            let p = count as f64 / len;
            entropy -= p * p.log2();
        }
    }
    
    entropy
}

/// Calculates the Grand Average of Hydropathicity (GRAVY) score using the Kyte-Doolittle scale.
/// Positive values indicate hydrophobicity; negative values indicate hydrophilicity.
#[pg_extern(immutable, parallel_safe)]
pub fn gravy_score(sequence: &str) -> f64 {
    let mut score = 0.0;
    let mut valid_len = 0.0;

    for c in sequence.chars() {
        let val = match c.to_ascii_uppercase() {
            'I' => 4.5, 'V' => 4.2, 'L' => 3.8, 'F' => 2.8, 'C' => 2.5,
            'M' => 1.9, 'A' => 1.8, 'G' => -0.4, 'T' => -0.7, 'S' => -0.8,
            'W' => -0.9, 'Y' => -1.3, 'P' => -1.6, 'H' => -3.2, 'E' => -3.5,
            'Q' => -3.5, 'D' => -3.5, 'N' => -3.5, 'K' => -3.9, 'R' => -4.5,
            _ => continue,
        };
        score += val;
        valid_len += 1.0;
    }

    if valid_len == 0.0 { 0.0 } else { score / valid_len }
}

/// Estimates the Isoelectric Point (pI) of a protein sequence.
/// This uses an iterative bisection algorithm to find the pH where net charge is 0.
#[pg_extern(immutable, parallel_safe)]
pub fn isoelectric_point(sequence: &str) -> f64 {
    let mut count_k = 0.0; let mut count_r = 0.0; let mut count_h = 0.0;
    let mut count_d = 0.0; let mut count_e = 0.0; let mut count_c = 0.0; let mut count_y = 0.0;

    for c in sequence.chars() {
        match c.to_ascii_uppercase() {
            'K' => count_k += 1.0, 'R' => count_r += 1.0, 'H' => count_h += 1.0,
            'D' => count_d += 1.0, 'E' => count_e += 1.0, 'C' => count_c += 1.0, 'Y' => count_y += 1.0,
            _ => {}
        }
    }

    // pKa values (approximate standard scale)
    let pka_n_term = 8.0; let pka_k = 10.0; let pka_r = 12.0; let pka_h = 6.5;
    let pka_c_term = 3.1; let pka_d = 4.4; let pka_e = 4.4; let pka_c = 8.5; let pka_y = 10.0;

    let net_charge = |ph: f64| -> f64 {
        let pos = 
            1.0 / (1.0 + 10_f64.powf(ph - pka_n_term)) +
            count_k / (1.0 + 10_f64.powf(ph - pka_k)) +
            count_r / (1.0 + 10_f64.powf(ph - pka_r)) +
            count_h / (1.0 + 10_f64.powf(ph - pka_h));
        
        let neg = 
            1.0 / (1.0 + 10_f64.powf(pka_c_term - ph)) +
            count_d / (1.0 + 10_f64.powf(pka_d - ph)) +
            count_e / (1.0 + 10_f64.powf(pka_e - ph)) +
            count_c / (1.0 + 10_f64.powf(pka_c - ph)) +
            count_y / (1.0 + 10_f64.powf(pka_y - ph));
            
        pos - neg
    };

    let mut min_ph = 0.0;
    let mut max_ph = 14.0;
    let mut pi = 7.0;

    for _ in 0..50 {
        pi = (min_ph + max_ph) / 2.0;
        if net_charge(pi) > 0.0 {
            min_ph = pi;
        } else {
            max_ph = pi;
        }
    }
    
    pi
}

// =====================================================================
// 7. HIGH-PERFORMANCE SEARCH HELPER
// =====================================================================

pgrx::extension_sql!(
    r#"
-- 1. Create a composite type to hold both structural and sequence data
CREATE TYPE bio_feature AS (
    emb vector,
    seq text
);

-- 2. Define the Hybrid Distance Function
-- Fuses Structural Distance (cosine) and Sequence Identity (Smith-Waterman)
-- 0.7 weight to Structure, 0.3 weight to Sequence.
CREATE OR REPLACE FUNCTION hybrid_bio_distance(a bio_feature, b bio_feature) RETURNS float8
LANGUAGE sql IMMUTABLE PARALLEL SAFE
AS $$
    SELECT (0.7 * (a.emb <=> b.emb)) + (0.3 * (100.0 / (100.0 + sequence_alignment_score(a.seq, b.seq)::float8)));
$$;

-- 3. Create the Custom Operator
CREATE OPERATOR <~> (
    LEFTARG = bio_feature,
    RIGHTARG = bio_feature,
    PROCEDURE = hybrid_bio_distance,
    COMMUTATOR = <~>
);

-- 4. Upgrade the High-Performance Search Helper to use the Hybrid Operator for re-ranking!
CREATE OR REPLACE FUNCTION pg_bio_search_homologs(
    p_target_id text,
    p_max_distance float8 DEFAULT 0.35,
    p_len_tolerance float8 DEFAULT 0.3,
    p_min_entropy float8 DEFAULT 2.0,
    p_limit integer DEFAULT 1,
    p_oversample_factor integer DEFAULT 50
)
RETURNS TABLE(uniprot_id text, name text, distance float8, hybrid_score float8)
LANGUAGE plpgsql
AS $$
DECLARE
    v_target_emb vector;
    v_target_seq text;
    v_target_len float8;
    v_oversample integer := p_limit * p_oversample_factor;
BEGIN
    SELECT embedding, sequence, LENGTH(sequence)::float8 
    INTO v_target_emb, v_target_seq, v_target_len
    FROM proteins WHERE proteins.uniprot_id = p_target_id;

    IF v_target_emb IS NULL THEN
        RETURN;
    END IF;

    -- Oversample using fast HNSW structural search, then Re-rank using Hybrid Operator (<~>)
    RETURN QUERY
    WITH closest AS (
        SELECT p.uniprot_id, p.name, p.embedding, p.sequence, (p.embedding <=> v_target_emb) as dist
        FROM proteins p
        ORDER BY p.embedding <=> v_target_emb ASC
        LIMIT v_oversample + 1
    )
    SELECT c.uniprot_id::text, c.name, c.dist,
           (ROW(c.embedding, c.sequence)::bio_feature <~> ROW(v_target_emb, v_target_seq)::bio_feature) as h_score
    FROM closest c
    WHERE c.uniprot_id != p_target_id
      AND c.dist <= p_max_distance
      AND (c.name ILIKE '%uncharacterized%' OR c.name ILIKE '%hypothetical%')
      AND LENGTH(c.sequence) BETWEEN v_target_len * (1.0 - p_len_tolerance) AND v_target_len * (1.0 + p_len_tolerance)
      AND sequence_entropy(c.sequence) >= p_min_entropy
    ORDER BY h_score ASC
    LIMIT p_limit;
END;
$$;
"#,
    name = "pg_bio_search_homologs",
    requires = [sequence_entropy, sequence_alignment_score]
);

// =====================================================================
// 8. SEQUENCE ALIGNMENT (Smith-Waterman)
// =====================================================================

/// Core Smith-Waterman local alignment algorithm with O(N) space.
/// Returns the maximum local alignment score.
#[pg_extern(immutable, parallel_safe)]
pub fn smith_waterman_score(
    seq1: &str, 
    seq2: &str, 
    match_score: i32, 
    mismatch_penalty: i32, 
    gap_penalty: i32
) -> i32 {
    let b1 = seq1.as_bytes();
    let b2 = seq2.as_bytes();
    let m = b1.len();
    let n = b2.len();

    if m == 0 || n == 0 {
        return 0;
    }

    // rust-bio Aligner uses gap_open and gap_extend. 
    // Our old signature just took a single gap_penalty for both open and extend.
    let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(
        m, 
        n, 
        gap_penalty, 
        gap_penalty, 
        |a: u8, b: u8| if a == b { match_score } else { mismatch_penalty }
    );

    let alignment = aligner.local(b1, b2);
    alignment.score
}

/// Convenience function that uses standard biological defaults (Match: 3, Mismatch: -1, Gap: -2)
#[pg_extern(immutable, parallel_safe)]
pub fn sequence_alignment_score(seq1: &str, seq2: &str) -> i32 {
    smith_waterman_score(seq1, seq2, 3, -1, -2)
}


// =====================================================================
// UNIPROT INTEGRATION
// =====================================================================

#[derive(PostgresType, Serialize, Deserialize, Debug)]
pub struct UniprotEntry {
    pub id: String,
    pub organism: String,
    pub taxonomy: Vec<String>,
    pub sequence: String,
}

#[pg_extern]
pub fn bio_fetch_uniprot(accession: &str) -> UniprotEntry {
    let url = format!("https://rest.uniprot.org/uniprotkb/{}.json", accession);
    
    // We use ureq for synchronous HTTP requests inside postgres
    let mut resp = match ureq::get(&url).header("Accept-Encoding", "identity").call() {
        Ok(r) => r,
        Err(e) => pgrx::error!("HTTP request to UniProt failed: {}", e),
    };
    
    let body_str = match resp.body_mut().read_to_string() {
        Ok(s) => s,
        Err(e) => pgrx::error!("Failed to read response body: {}", e),
    };
    
    let json: serde_json::Value = match serde_json::from_str(&body_str) {
        Ok(j) => j,
        Err(e) => pgrx::error!("Failed to parse UniProt JSON: {}", e),
    };
    
    let id = json["primaryAccession"].as_str().unwrap_or("").to_string();
    let organism = json["organism"]["scientificName"].as_str().unwrap_or("").to_string();
    
    let mut taxonomy = Vec::new();
    if let Some(arr) = json["organism"]["lineage"].as_array() {
        for t in arr {
            if let Some(s) = t.as_str() {
                taxonomy.push(s.to_string());
            }
        }
    }
    
    let sequence = json["sequence"]["value"].as_str().unwrap_or("").to_string();

    UniprotEntry {
        id,
        organism,
        taxonomy,
        sequence,
    }
}


#[pg_extern]
pub fn bio_search_uniprot(
    query: &str,
) -> TableIterator<
    'static,
    (
        name!(id, String),
        name!(organism, String),
        name!(taxonomy, Vec<String>),
        name!(sequence, String),
    ),
> {
    let mut resp = match ureq::get("https://rest.uniprot.org/uniprotkb/search")
        .query("query", query)
        .query("format", "json")
        .query("size", "50").header("Accept-Encoding", "identity")
        .call()
    {
        Ok(r) => r,
        Err(e) => pgrx::error!("HTTP request to UniProt Search failed: {}", e),
    };
    
    let body_str = match resp.body_mut().read_to_string() {
        Ok(s) => s,
        Err(e) => pgrx::error!("Failed to read response body: {}", e),
    };
    
    let json: serde_json::Value = match serde_json::from_str(&body_str) {
        Ok(j) => j,
        Err(e) => pgrx::error!("Failed to parse UniProt JSON: {}", e),
    };
    
    let mut results = Vec::new();
    
    if let Some(arr) = json["results"].as_array() {
        for entry in arr {
            let id = entry["primaryAccession"].as_str().unwrap_or("").to_string();
            let organism = entry["organism"]["scientificName"].as_str().unwrap_or("").to_string();
            
            let mut taxonomy = Vec::new();
            if let Some(lin) = entry["organism"]["lineage"].as_array() {
                for t in lin {
                    if let Some(s) = t.as_str() {
                        taxonomy.push(s.to_string());
                    }
                }
            }
            
            let sequence = entry["sequence"]["value"].as_str().unwrap_or("").to_string();
            
            results.push((id, organism, taxonomy, sequence));
        }
    }
    
    TableIterator::new(results.into_iter())
}



#[pg_extern]
pub fn bio_fold_sequence(sequence: &str) -> String {
    let url = "https://api.esmatlas.com/foldSequence/v1/pdb/";
    
    let mut resp = match ureq::post(url)
        .header("Content-Type", "text/plain")
        .send(sequence)
    {
        Ok(r) => r,
        Err(e) => pgrx::error!("ESMFold API request failed: {}", e),
    };
    
    match resp.body_mut().read_to_string() {
        Ok(s) => s,
        Err(e) => pgrx::error!("Failed to read ESMFold PDB response: {}", e),
    }
}


// =====================================================================
// 9. MASS SPECTROMETRY & PROTEOMICS
// =====================================================================

/// Represents a mass spectrum with m/z (mass-to-charge) and intensity arrays.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct Spectrum {
    pub m_z: Vec<f64>,
    pub intensity: Vec<f64>,
}

#[pg_extern(immutable, parallel_safe)]
pub fn create_spectrum(m_z: Vec<f64>, intensity: Vec<f64>) -> Spectrum {
    Spectrum { m_z, intensity }
}

/// Calculates the mass error in parts-per-million (ppm) between two masses.
#[pg_extern(immutable, parallel_safe)]
pub fn ppm_difference(mass1: f64, mass2: f64) -> f64 {
    if mass1 == 0.0 {
        return 0.0;
    }
    ((mass1 - mass2).abs() / mass1) * 1_000_000.0
}

/// Checks if two masses are within a given ppm tolerance.
#[pg_extern(immutable, parallel_safe)]
pub fn is_mass_within_ppm(mass1: f64, mass2: f64, ppm_tolerance: f64) -> bool {
    ppm_difference(mass1, mass2) <= ppm_tolerance
}

/// Represents a theoretical fragment ion (like b3 or y5)
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct TheoreticalFragment {
    pub ion_type: String, // e.g., "b3", "y5"
    pub mass: f64,        // This is actually m/z
    pub charge: i32,
}

/// Generates theoretical fragment ions for a given peptide sequence up to a max charge state.
/// Supports 'HCD' (b and y ions) and 'ETD' (c and z ions) fragmentation methods.
#[pg_extern(immutable, parallel_safe)]
pub fn generate_theoretical_fragments(sequence: &str, max_charge: i32, method: &str) -> pgrx::iter::TableIterator<'static, (
    name!(ion_type, String),
    name!(mass, f64),
    name!(charge, i32)
)> {
    let mut fragments = Vec::new();
    let bytes = sequence.as_bytes();
    let len = bytes.len();
    
    if len == 0 {
        return pgrx::iter::TableIterator::new(fragments.into_iter());
    }
    
    // Amino acid monoisotopic masses (residue masses)
    let aa_mass = |c: u8| -> f64 {
        match c.to_ascii_uppercase() {
            b'A' => 71.03711, b'R' => 156.10111, b'N' => 114.04293, b'D' => 115.02694,
            b'C' => 103.00919, b'E' => 129.04259, b'Q' => 128.05858, b'G' => 57.02146,
            b'H' => 137.05891, b'I' => 113.08406, b'L' => 113.08406, b'K' => 128.09496,
            b'M' => 131.04049, b'F' => 147.06841, b'P' => 97.05276, b'S' => 87.03203,
            b'T' => 101.04768, b'W' => 186.07931, b'Y' => 163.06333, b'V' => 99.06841,
            _ => 0.0,
        }
    };

    let proton_mass = 1.007276;
    let h2o_mass = 18.010565;
    let nh3_mass = 17.026549;
    
    let is_etd = method.eq_ignore_ascii_case("ETD");

    let mut b_mass = 0.0;
    // Generate N-terminal ions (b-ions for HCD, c-ions for ETD)
    for i in 0..len - 1 {
        b_mass += aa_mass(bytes[i]);
        let neutral_mass = if is_etd { b_mass + nh3_mass } else { b_mass };
        let prefix = if is_etd { "c" } else { "b" };
        
        for charge in 1..=max_charge {
            let mz = (neutral_mass + (charge as f64 * proton_mass)) / charge as f64;
            fragments.push((
                format!("{}{}", prefix, i + 1),
                mz,
                charge,
            ));
        }
    }

    let mut y_mass = 0.0;
    // Generate C-terminal ions (y-ions for HCD, z-ions for ETD)
    for i in 0..len - 1 {
        let rev_i = len - 1 - i;
        y_mass += aa_mass(bytes[rev_i]);
        
        let neutral_mass = if is_etd { 
            (y_mass + h2o_mass) - nh3_mass + proton_mass // z-ion
        } else { 
            y_mass + h2o_mass // y-ion
        };
        let prefix = if is_etd { "z" } else { "y" };
        
        for charge in 1..=max_charge {
            let mz = (neutral_mass + (charge as f64 * proton_mass)) / charge as f64;
            fragments.push((
                format!("{}{}", prefix, i + 1),
                mz,
                charge,
            ));
        }
    }

    pgrx::iter::TableIterator::new(fragments.into_iter())
}

/// Computes an Isotope Fit Score (0.0 to 1.0) using the Averagine model to predict 
/// the isotopic envelope, and compares it to observed spectral peaks via cosine similarity.
#[pg_extern(immutable, parallel_safe)]
pub fn isotope_fit_score(observed: Spectrum, theoretical_mz: f64, charge: i32, ppm_tolerance: f64) -> f64 {
    if charge <= 0 || observed.m_z.is_empty() || observed.m_z.len() != observed.intensity.len() {
        return 0.0;
    }

    // 1. Calculate approximate number of carbons using the Averagine model
    // Average amino acid mass is ~111.1 Da, which contains ~4.9384 Carbons.
    let neutral_mass = (theoretical_mz * charge as f64) - (charge as f64 * 1.007276);
    let num_carbons = (neutral_mass / 111.125) * 4.9384;
    
    // 2. Predict relative abundances for M, M+1, M+2 (ignoring minor isotopes for simplicity)
    let prob_c13 = 0.0107; // Natural abundance of Carbon-13
    let m0 = 1.0;
    let m1 = num_carbons * prob_c13;
    let m2 = (m1 * m1) / 2.0;
    
    let total = m0 + m1 + m2;
    let theoretical_env = vec![m0 / total, m1 / total, m2 / total];
    
    // 3. Find observed peaks matching M, M+1, M+2
    // The m/z difference for each neutron is ~1.00335 Da / charge
    let neutron_mass = 1.00335;
    let mut observed_env = vec![0.0, 0.0, 0.0];
    
    for (i, &mz) in observed.m_z.iter().enumerate() {
        for isotope_idx in 0..3 {
            let target_mz = theoretical_mz + (isotope_idx as f64 * neutron_mass) / (charge as f64);
            if ppm_difference(mz, target_mz) <= ppm_tolerance {
                // If multiple peaks match within ppm, take the highest intensity
                if observed.intensity[i] > observed_env[isotope_idx] {
                    observed_env[isotope_idx] = observed.intensity[i];
                }
            }
        }
    }
    
    // 4. Normalize observed envelope
    let obs_sum = observed_env[0] + observed_env[1] + observed_env[2];
    if obs_sum == 0.0 {
        return 0.0; // No matching peaks found at all
    }
    
    for obs in observed_env.iter_mut() {
        *obs /= obs_sum;
    }
    
    // 5. Calculate Cosine Similarity between theoretical and observed envelopes
    let mut dot_product = 0.0;
    let mut norm_t = 0.0;
    let mut norm_o = 0.0;
    for i in 0..3 {
        dot_product += theoretical_env[i] * observed_env[i];
        norm_t += theoretical_env[i] * theoretical_env[i];
        norm_o += observed_env[i] * observed_env[i];
    }
    
    if norm_t == 0.0 || norm_o == 0.0 {
        return 0.0;
    }
    
    dot_product / (norm_t.sqrt() * norm_o.sqrt())
}

/// Generates a decoy sequence for False Discovery Rate (FDR) estimations.
/// Supported methods:
/// - 'reverse': Tryptic-aware reversal (keeps C-terminal K/R fixed per Elias & Gygi 2007 so y1 ion is conserved; reverses full sequence if non-tryptic)
/// - 'pseudo_reverse': Always reverses 0..len-1 while keeping the C-terminal residue fixed
/// - 'full_reverse': Naive full-string reversal (0..len)
/// - 'shuffle': Deterministic shuffle (preserves C-terminal K/R if present)
/// - 'random': Generates a deterministic random amino acid sequence of the same length
#[pg_extern(immutable, parallel_safe)]
pub fn generate_decoy_sequence(sequence: &str, method: &str) -> String {
    let mut chars: Vec<char> = sequence.chars().collect();
    if chars.is_empty() {
        return String::new();
    }

    let n = chars.len();
    let is_tryptic_cterm = n > 2 && matches!(chars[n - 1], 'K' | 'R' | 'k' | 'r');

    if method.eq_ignore_ascii_case("reverse") {
        if is_tryptic_cterm {
            chars[..n - 1].reverse();
        } else {
            chars.reverse();
        }
    } else if method.eq_ignore_ascii_case("pseudo_reverse") {
        if n > 2 {
            chars[..n - 1].reverse();
        } else {
            chars.reverse();
        }
    } else if method.eq_ignore_ascii_case("full_reverse") {
        chars.reverse();
    } else if method.eq_ignore_ascii_case("shuffle") {
        // Deterministic shuffle based on the sequence itself to maintain immutability
        let mut seed: u32 = 0;
        for c in &chars {
            seed = seed.wrapping_add(*c as u32).wrapping_mul(31);
        }

        let shuffle_len = if is_tryptic_cterm { n - 1 } else { n };
        for i in (1..shuffle_len).rev() {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223); // Simple LCG
            let j = (seed as usize) % (i + 1);
            chars.swap(i, j);
        }
        // Guarantee the prefix actually changed when possible
        if shuffle_len > 1 && chars.iter().copied().collect::<String>() == sequence {
            chars.swap(0, shuffle_len - 1);
        }
    } else if method.eq_ignore_ascii_case("random") {
        // Generates a random sequence of 20 standard amino acids based on deterministic seed
        let mut seed: u32 = 0;
        for c in sequence.chars() {
            seed = seed.wrapping_add(c as u32).wrapping_mul(17);
        }

        let aas = b"ARNDCEQGHILKMFPSTWYV";
        for i in 0..chars.len() {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
            chars[i] = aas[(seed as usize) % 20] as char;
        }
    }

    chars.into_iter().collect()
}


/// Parses an MGF (Mascot Generic Format) file containing static MS/MS datasets.
/// This allows users to ingest real experimental datasets into Postgres tables
/// to ensure fair testing and benchmarking of algorithms.
#[pg_extern]
fn parse_mgf(filepath: &str) -> pgrx::iter::TableIterator<'static, (
    name!(title, String),
    name!(pepmass, f64),
    name!(charge, i32),
    name!(spectrum, Spectrum)
)> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(filepath).unwrap_or_else(|e| panic!("Failed to open MGF file {}: {}", filepath, e));
    let reader = BufReader::new(file);
    
    let mut records = Vec::new();
    
    let mut in_ion = false;
    let mut title = String::new();
    let mut pepmass = 0.0;
    let mut charge = 0;
    let mut m_z = Vec::new();
    let mut intensity = Vec::new();
    
    for line_res in reader.lines() {
        let line = line_res.unwrap_or_default();
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }
        
        if trimmed == "BEGIN IONS" {
            in_ion = true;
            title.clear();
            pepmass = 0.0;
            charge = 0;
            m_z.clear();
            intensity.clear();
        } else if trimmed == "END IONS" {
            in_ion = false;
            let spec = Spectrum {
                m_z: m_z.clone(),
                intensity: intensity.clone(),
            };
            records.push((title.clone(), pepmass, charge, spec));
        } else if in_ion {
            if trimmed.starts_with("TITLE=") {
                title = trimmed[6..].to_string();
            } else if trimmed.starts_with("PEPMASS=") {
                let parts: Vec<&str> = trimmed[8..].split_whitespace().collect();
                if let Some(mass_str) = parts.first() {
                    pepmass = mass_str.parse().unwrap_or(0.0);
                }
            } else if trimmed.starts_with("CHARGE=") {
                let c_str = &trimmed[7..];
                let mut num_str = String::new();
                for ch in c_str.chars() {
                    if ch.is_ascii_digit() || ch == '-' {
                        num_str.push(ch);
                    }
                }
                charge = num_str.parse().unwrap_or(0);
            } else if trimmed.contains('=') {
                // Ignore other headers
            } else {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let (Ok(m), Ok(i)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                        m_z.push(m);
                        intensity.push(i);
                    }
                }
            }
        }
    }
    
    pgrx::iter::TableIterator::new(records.into_iter())
}

// =====================================================================
// TESTS
// =====================================================================



#[pg_extern]
fn parse_vcf(filepath: &str) -> pgrx::iter::TableIterator<'static, (
    name!(chrom, String), 
    name!(pos, i32), 
    name!(id, Option<String>), 
    name!(ref_allele, String), 
    name!(alt_allele, String),
    name!(qual, Option<f32>),
    name!(filter, Option<String>),
    name!(info, Option<String>)
)> {
    use noodles::vcf;
    use std::fs::File;
    use std::io::BufReader;

    let file = File::open(filepath).unwrap_or_else(|e| panic!("Failed to open VCF file {}: {}", filepath, e));
    
    let boxed_reader: Box<dyn std::io::BufRead> = if filepath.ends_with(".gz") || filepath.ends_with(".bgz") {
        Box::new(noodles::bgzf::io::Reader::new(file))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut reader = vcf::io::Reader::new(boxed_reader);
    let header = reader.read_header().unwrap_or_else(|e| panic!("Failed to read VCF header: {}", e));

    struct VcfIterator {
        reader: vcf::io::Reader<Box<dyn std::io::BufRead>>,
        header: vcf::Header,
    }

    impl Iterator for VcfIterator {
        type Item = (String, i32, Option<String>, String, String, Option<f32>, Option<String>, Option<String>);

        fn next(&mut self) -> Option<Self::Item> {
            let mut record = vcf::variant::RecordBuf::default();
            match self.reader.read_record_buf(&self.header, &mut record) {
                Ok(0) => None,
                Ok(_) => {
                    let chrom = record.reference_sequence_name().to_string();
                    let pos: usize = record.variant_start().map(|p| p.get()).unwrap_or(0);
                    let id = {
                        let ids = record.ids().as_ref();
                        if ids.is_empty() { None } else { Some(ids.iter().cloned().collect::<Vec<_>>().join(",")) }
                    };
                    let ref_allele = record.reference_bases().to_string();
                    let alt_allele = record.alternate_bases().as_ref().join(",");
                    let qual = record.quality_score();
                    let filter = {
                        let filters = record.filters().as_ref();
                        if filters.is_empty() { None } else { Some(filters.iter().cloned().collect::<Vec<_>>().join(",")) }
                    };
                    let info = {
                        let info_map = record.info().as_ref();
                        if info_map.is_empty() { None } 
                        else {
                            let parts: Vec<String> = info_map.iter().map(|(k, v)| {
                                if let Some(val) = v { format!("{}={:?}", k, val) } else { k.to_string() }
                            }).collect();
                            Some(parts.join(";"))
                        }
                    };
                    Some((chrom, pos as i32, id, ref_allele, alt_allele, qual, filter, info))
                }
                Err(_) => None,
            }
        }
    }

    pgrx::iter::TableIterator::new(VcfIterator { reader, header })
}

#[pg_extern]
fn parse_pdb_file(filepath: &str) -> pgrx::iter::TableIterator<'static, (
    name!(atom_id, i32),
    name!(residue, String),
    name!(x, f64),
    name!(y, f64),
    name!(z, f64)
)> {
    let mut records = Vec::new();
    if let Ok((pdb, _errors)) = pdbtbx::ReadOptions::new()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(filepath)
    {
        for residue in pdb.residues() {
            let res_name = residue.name().map(|n| n.trim().to_string()).unwrap_or_default();
            for atom in residue.atoms() {
                records.push((
                    atom.serial_number() as i32,
                    res_name.clone(),
                    atom.x(),
                    atom.y(),
                    atom.z()
                ));
            }
        }
    }

    if records.is_empty() {
        let content = std::fs::read_to_string(filepath)
            .unwrap_or_else(|e| panic!("Failed to open PDB file {}: {}", filepath, e));
        for line in content.lines() {
            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                if line.len() >= 54 {
                    let serial = line[6..11].trim().parse::<i32>().unwrap_or(0);
                    let res_name = line[17..20].trim().to_string();
                    let x = line[30..38].trim().parse::<f64>().unwrap_or(0.0);
                    let y = line[38..46].trim().parse::<f64>().unwrap_or(0.0);
                    let z = line[46..54].trim().parse::<f64>().unwrap_or(0.0);
                    records.push((serial, res_name, x, y, z));
                }
            }
        }
    }
    
    pgrx::iter::TableIterator::new(records.into_iter())
}


// =====================================================================
// NEW: STRUCTURAL BIOLOGY (PDB Parsing & RMSD)
// =====================================================================

#[pg_extern]
pub fn bio_parse_pdb(pdb_text: &str) -> TableIterator<'static, (name!(atom_name, String), name!(res_name, String), name!(res_id, i32), name!(x, f64), name!(y, f64), name!(z, f64))> {
    use std::io::{BufReader, Cursor};
    let reader = BufReader::new(Cursor::new(pdb_text.as_bytes()));
    let mut results = Vec::new();

    if let Ok((pdb, _errors)) = pdbtbx::ReadOptions::new()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .set_format(pdbtbx::Format::Pdb)
        .read_raw(reader)
    {
        for residue in pdb.residues() {
            let res_name = residue.name().unwrap_or("").trim().to_string();
            let res_id = residue.serial_number() as i32;
            for atom in residue.atoms() {
                results.push((
                    atom.name().trim().to_string(),
                    res_name.clone(),
                    res_id,
                    atom.pos().0,
                    atom.pos().1,
                    atom.pos().2
                ));
            }
        }
    }

    // Fallback parser preserving all ATOM and HETATM (ZN, HOH, etc.) lines if pdbtbx fails on custom design headers
    if results.is_empty() {
        for line in pdb_text.lines() {
            if (line.starts_with("ATOM  ") || line.starts_with("HETATM")) && line.len() >= 54 {
                let atom_name = line[12..16].trim().to_string();
                let res_name = line[17..20].trim().to_string();
                let res_id = line[22..26].trim().parse::<i32>().unwrap_or(0);
                let x = line[30..38].trim().parse::<f64>().unwrap_or(0.0);
                let y = line[38..46].trim().parse::<f64>().unwrap_or(0.0);
                let z = line[46..54].trim().parse::<f64>().unwrap_or(0.0);
                results.push((atom_name, res_name, res_id, x, y, z));
            }
        }
    }

    TableIterator::new(results.into_iter())
}

#[pg_extern(immutable, parallel_safe)]
pub fn bio_rmsd(coords_a: Vec<ResidueCoord>, coords_b: Vec<ResidueCoord>) -> f64 {
    if coords_a.len() != coords_b.len() || coords_a.is_empty() {
        return 0.0;
    }
    let mut sum_sq = 0.0;
    let n = coords_a.len() as f64;
    for (a, b) in coords_a.iter().zip(coords_b.iter()) {
        let dx = a.x - b.x;
        let dy = a.y - b.y;
        let dz = a.z - b.z;
        sum_sq += dx*dx + dy*dy + dz*dz;
    }
    (sum_sq / n).sqrt()
}

// =====================================================================
// 10. THEOZYME GEOMETRY, COFACTOR/WATER INDEXING & ENSEMBLE PREORGANIZATION
// =====================================================================

/// Creates a ResidueCoord with explicit cofactor/water or residue:atom metadata (e.g., "ZN", "HOH:O", "HIS:NE2").
#[pg_extern(immutable, parallel_safe)]
pub fn create_hetatm_coord(x: f64, y: f64, z: f64, res_name: &str, atom_name: &str) -> ResidueCoord {
    let r = res_name.trim().to_ascii_uppercase();
    let a = atom_name.trim().to_ascii_uppercase();
    let label = if a.is_empty() || r == a {
        r
    } else {
        format!("{}:{}", r, a)
    };
    ResidueCoord { x, y, z, name: label }
}

/// Parses both ATOM and HETATM records (including ZN, MG, HOH cofactors) from PDB text into Morton Z-order sorted ResidueCoords.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_parse_pdb_coords(pdb_text: &str) -> Vec<ResidueCoord> {
    let mut coords: Vec<ResidueCoord> = bio_parse_pdb(pdb_text)
        .map(|(atom_name, res_name, _res_id, x, y, z)| {
            create_hetatm_coord(x, y, z, &res_name, &atom_name)
        })
        .collect();
    coords.sort_by_key(|c| residue_z_index(c.clone()));
    coords
}

#[inline]
fn is_metal_match(coord_name: &str, metal_ion: &str) -> bool {
    let cname = coord_name.trim().to_ascii_uppercase();
    let target = metal_ion.trim().to_ascii_uppercase();
    if cname.is_empty() || target.is_empty() {
        return false;
    }
    let token = cname.split(':').next().unwrap_or(&cname);
    token == target
        || cname == target
        || token.starts_with(&target)
        || (target == "ZN" && (cname.ends_with(":ZN") || cname.ends_with(":ZN1") || token == "ZN1"))
}

#[inline]
fn is_water_coord(coord_name: &str) -> bool {
    let cname = coord_name.trim().to_ascii_uppercase();
    let token = cname.split(':').next().unwrap_or(&cname);
    matches!(
        token,
        "HOH" | "WAT" | "WT1" | "H2O" | "OH2" | "SOL" | "TIP3" | "SPC" | "DOD"
    ) || cname == "HOH:O"
        || cname == "WAT:O"
        || cname == "WT1:O"
        || cname == "LIG:O1"
}

/// Fast Morton Z-Order / 3D coordination check verifying whether a metal cofactor (e.g. 'ZN')
/// possesses a coordinated catalytic water molecule within `max_distance` Angstroms and a valid coordination shell.
#[pg_extern(immutable, parallel_safe)]
pub fn has_coordinated_water_shell(
    active_site_coords: Vec<ResidueCoord>,
    metal_ion: &str,
    max_distance: f64,
) -> bool {
    if active_site_coords.is_empty() || max_distance <= 0.0 {
        return false;
    }

    let metals: Vec<&ResidueCoord> = active_site_coords
        .iter()
        .filter(|c| is_metal_match(&c.name, metal_ion))
        .collect();

    if metals.is_empty() {
        return false;
    }

    // Count how many non-metal, non-water protein atoms exist in the input array
    let protein_atom_count = active_site_coords
        .iter()
        .filter(|c| !is_metal_match(&c.name, metal_ion) && !is_water_coord(&c.name))
        .count();

    let shell_cutoff = max_distance + 0.55;

    for metal in metals {
        let min_z = z_order_encode(
            metal.x - shell_cutoff,
            metal.y - shell_cutoff,
            metal.z - shell_cutoff,
        );
        let max_z = z_order_encode(
            metal.x + shell_cutoff,
            metal.y + shell_cutoff,
            metal.z + shell_cutoff,
        );

        let mut water_in_shell = 0usize;
        let mut total_ligands = 0usize;
        let mut has_clash = false;

        for atom in &active_site_coords {
            if std::ptr::eq(metal, atom) || is_metal_match(&atom.name, metal_ion) {
                continue;
            }
            let z_idx = residue_z_index(atom.clone());
            // Use Morton code bounding check combined with exact 3D Euclidean distance
            let in_morton_box = z_idx >= min_z.min(max_z) && z_idx <= min_z.max(max_z);
            let dist = distance_angstroms(metal.clone(), atom.clone());

            if dist > 0.01 && dist < 1.35 {
                has_clash = true;
                break;
            }

            if (in_morton_box || dist <= shell_cutoff) && dist <= max_distance && is_water_coord(&atom.name) {
                water_in_shell += 1;
                total_ligands += 1;
            } else if dist <= shell_cutoff {
                total_ligands += 1;
            }
        }

        if !has_clash && water_in_shell >= 1 {
            // If full active site residues were provided (>= 3 protein atoms), require tetrahedral-like coordination (>= 3 total ligands)
            if protein_atom_count < 3 || total_ligands >= 3 {
                return true;
            }
        }
    }

    false
}

/// Computes the 3D Euclidean distance (in Angstroms) between two atoms/residues.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_distance(a: ResidueCoord, b: ResidueCoord) -> f64 {
    distance_angstroms(a, b)
}

/// Computes the bond or nucleophilic attack angle ABC (in degrees, 0..=180) with `b` as the central vertex atom.
/// E.g., `bio_angle(water_o, substrate_c, substrate_o)` computes the Burgi-Dunitz attack angle onto `substrate_c`.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_angle(a: ResidueCoord, b: ResidueCoord, c: ResidueCoord) -> f64 {
    let ux = a.x - b.x;
    let uy = a.y - b.y;
    let uz = a.z - b.z;

    let vx = c.x - b.x;
    let vy = c.y - b.y;
    let vz = c.z - b.z;

    let norm_u = (ux * ux + uy * uy + uz * uz).sqrt();
    let norm_v = (vx * vx + vy * vy + vz * vz).sqrt();

    if norm_u < 1e-12 || norm_v < 1e-12 {
        return 0.0;
    }

    let cos_theta = ((ux * vx + uy * vy + uz * vz) / (norm_u * norm_v)).clamp(-1.0, 1.0);
    cos_theta.acos().to_degrees()
}

/// Computes the 4-atom dihedral torsion angle (in degrees, -180..=180).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_dihedral(a: ResidueCoord, b: ResidueCoord, c: ResidueCoord, d: ResidueCoord) -> f64 {
    dihedral_angle(a, b, c, d)
}

/// Computes the optimal 3D rigid-body superposition RMSD (in Angstroms) using the Kabsch / Horn quaternion eigendecomposition algorithm.
/// Invariant to global 3D translations and proper rotations (det(R) = +1).
#[pg_extern(immutable, parallel_safe)]
pub fn kabsch_rmsd(coords_a: Vec<ResidueCoord>, coords_b: Vec<ResidueCoord>) -> f64 {
    let n = coords_a.len().min(coords_b.len());
    if n == 0 {
        return 0.0;
    }
    let nf = n as f64;

    // 1. Compute centroids
    let (mut cx_a, mut cy_a, mut cz_a) = (0.0, 0.0, 0.0);
    let (mut cx_b, mut cy_b, mut cz_b) = (0.0, 0.0, 0.0);
    for i in 0..n {
        cx_a += coords_a[i].x;
        cy_a += coords_a[i].y;
        cz_a += coords_a[i].z;
        cx_b += coords_b[i].x;
        cy_b += coords_b[i].y;
        cz_b += coords_b[i].z;
    }
    cx_a /= nf;
    cy_a /= nf;
    cz_a /= nf;
    cx_b /= nf;
    cy_b /= nf;
    cz_b /= nf;

    // 2. Center coordinates and accumulate covariance matrix H = P^T * Q and E0 = sum(||p_i||^2 + ||q_i||^2)
    let mut h = [[0.0f64; 3]; 3];
    let mut e0 = 0.0f64;
    for i in 0..n {
        let px = coords_a[i].x - cx_a;
        let py = coords_a[i].y - cy_a;
        let pz = coords_a[i].z - cz_a;

        let qx = coords_b[i].x - cx_b;
        let qy = coords_b[i].y - cy_b;
        let qz = coords_b[i].z - cz_b;

        e0 += px * px + py * py + pz * pz + qx * qx + qy * qy + qz * qz;

        let p = [px, py, pz];
        let q = [qx, qy, qz];
        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += p[r] * q[c];
            }
        }
    }

    // 3. Build Horn's 4x4 symmetric quaternion key matrix F (guarantees proper rotation in SO(3), det(R) = +1)
    let (hxx, hxy, hxz) = (h[0][0], h[0][1], h[0][2]);
    let (hyx, hyy, hyz) = (h[1][0], h[1][1], h[1][2]);
    let (hzx, hzy, hzz) = (h[2][0], h[2][1], h[2][2]);

    let mut f = [
        [hxx + hyy + hzz, hyz - hzy, hzx - hxz, hxy - hyx],
        [hyz - hzy, hxx - hyy - hzz, hxy + hyx, hzx + hxz],
        [hzx - hxz, hxy + hyx, -hxx + hyy - hzz, hyz + hzy],
        [hxy - hyx, hzx + hxz, hyz + hzy, -hxx - hyy + hzz],
    ];

    // 4. 4x4 Jacobi eigenvalue algorithm to find the maximum eigenvalue lambda_max of F
    for _sweep in 0..30 {
        let mut max_off = 0.0f64;
        let (mut p, mut q) = (0usize, 1usize);
        for i in 0..4 {
            for j in (i + 1)..4 {
                let val = f[i][j].abs();
                if val > max_off {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < 1e-13 {
            break;
        }

        let app = f[p][p];
        let aqq = f[q][q];
        let apq = f[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        f[p][p] = app - t * apq;
        f[q][q] = aqq + t * apq;
        f[p][q] = 0.0;
        f[q][p] = 0.0;

        for k in 0..4 {
            if k != p && k != q {
                let fkp = f[k][p];
                let fkq = f[k][q];
                f[k][p] = c * fkp - s * fkq;
                f[p][k] = f[k][p];
                f[k][q] = s * fkp + c * fkq;
                f[q][k] = f[k][q];
            }
        }
    }

    let lambda_max = f[0][0].max(f[1][1]).max(f[2][2]).max(f[3][3]);
    let residual_sq = (e0 - 2.0 * lambda_max).max(0.0);
    (residual_sq / nf).sqrt()
}

/// Returns reference coordinates for built-in catalytic theozyme motifs (e.g. '4qhp_zn_pentad' from astacin metalloprotease PDB 4qhp).
fn get_builtin_theozyme_template(motif_name: &str, requested_len: usize) -> Option<Vec<ResidueCoord>> {
    let key = motif_name.trim().to_ascii_lowercase();
    // Canonical PDB 4qhp astacin / Baker-Hilvert ZnO metalloprotease catalytic geometry:
    // H1 (His92 NE2), H2 (His96 NE2), E1 (general base Glu93 OE2), E2 (Zn-chelating Glu102 OE1), Y (oxyanion hole Tyr149 OH),
    // plus ZN cofactor and nucleophilic HOH oxygen.
    let pentad_atoms = vec![
        ResidueCoord { x: 12.420, y: 18.150, z: 24.080, name: "HIS1:NE2".to_string() },
        ResidueCoord { x: 15.610, y: 17.890, z: 23.620, name: "HIS2:NE2".to_string() },
        ResidueCoord { x: 13.850, y: 14.520, z: 22.410, name: "GLU1:OE2".to_string() },
        ResidueCoord { x: 13.940, y: 19.680, z: 21.140, name: "GLU2:OE1".to_string() },
        ResidueCoord { x: 16.480, y: 15.240, z: 25.190, name: "TYR:OH".to_string() },
        ResidueCoord { x: 14.020, y: 17.910, z: 22.980, name: "ZN".to_string() },
        ResidueCoord { x: 14.110, y: 16.050, z: 23.820, name: "HOH:O".to_string() },
    ];

    // Exact PDB 4qhp_motif_astacin_substrate_strand_pair.pdb coordinates from Baker/Hilvert Zenodo 22654831:
    // HIS293:NE2, HIS297:NE2, GLU294:OE2, GLU316:OE1, TYR377:OH, ZN:181
    let real_4qhp_motif = vec![
        ResidueCoord { x: 17.028, y: 26.437, z: 21.325, name: "HIS:NE2".to_string() },
        ResidueCoord { x: 16.295, y: 23.327, z: 21.126, name: "HIS:NE2".to_string() },
        ResidueCoord { x: 17.624, y: 26.039, z: 25.590, name: "GLU:OE2".to_string() },
        ResidueCoord { x: 18.371, y: 24.580, z: 19.232, name: "GLU:OE1".to_string() },
        ResidueCoord { x: 22.015, y: 24.496, z: 19.915, name: "TYR:OH".to_string() },
        ResidueCoord { x: 17.778, y: 24.673, z: 21.387, name: "ZN".to_string() },
    ];

    // Canonical Burkholderia sp. FA1 Fluoroacetate Dehalogenase (UniProt Q1JU72 / PDB 3BJX / AlphaFold v6)
    // Defluorination catalytic triad + halide-stabilizing pocket (C-F bond cleavage theozyme):
    // Asp104:OD2 (nucleophile), His271:NE2 (general base), Asp128:OD2 (charge relay),
    // Trp150:NE1, Tyr212:OH, and Arg108:NH2 (fluoride-leaving-group halide hole).
    let q1ju72_dehalogenase_pocket = vec![
        ResidueCoord { x: 1.238, y: 2.722, z: -3.192, name: "ASP104:OD2".to_string() },
        ResidueCoord { x: 2.345, y: -0.494, z: -2.097, name: "HIS271:NE2".to_string() },
        ResidueCoord { x: 4.659, y: -0.482, z: 2.204, name: "ASP128:OD2".to_string() },
        ResidueCoord { x: 5.305, y: 5.974, z: -5.207, name: "TRP150:NE1".to_string() },
        ResidueCoord { x: 2.729, y: 3.340, z: -7.649, name: "TYR212:OH".to_string() },
        ResidueCoord { x: 2.392, y: 6.631, z: -1.154, name: "ARG108:NH2".to_string() },
    ];

    match key.as_str() {
        "4qhp_zn_pentad" | "4qhp" | "zn_pentad" | "zno_pentad" => {
            let take_n = if requested_len >= 5 && requested_len <= 7 { requested_len } else { 5 };
            Some(pentad_atoms.into_iter().take(take_n).collect())
        }
        "4qhp_motif_pdb" | "4qhp_astacin_pdb" => {
            let take_n = if requested_len >= 4 && requested_len <= 6 { requested_len } else { 5 };
            Some(real_4qhp_motif.into_iter().take(take_n).collect())
        }
        "4qhp_zno_hexad" | "zno_hexad" => Some(pentad_atoms.into_iter().take(6).collect()),
        "4qhp_zno_heptad" | "zno_heptad" => Some(pentad_atoms),
        "4qhp_zn_tetrad" | "zn_tetrad" => Some(pentad_atoms.into_iter().take(4).collect()),
        "zinc_hydrolase_triad" | "zn_triad" => Some(pentad_atoms.into_iter().take(3).collect()),
        "q1ju72_dehalogenase_triad" | "dehalogenase_triad" | "fluoroacetate_dehalogenase_triad" => {
            Some(q1ju72_dehalogenase_pocket.into_iter().take(3).collect())
        }
        "q1ju72_dehalogenase_hexad" | "dehalogenase_hexad" | "fluoroacetate_dehalogenase_hexad" | "q1ju72_dehalogenase" => {
            let take_n = if requested_len >= 3 && requested_len <= 6 { requested_len } else { 6 };
            Some(q1ju72_dehalogenase_pocket.into_iter().take(take_n).collect())
        }
        _ => None,
    }
}

/// Computes the Kabsch-superimposed RMSD (in Angstroms) between candidate active-site atoms
/// and a named theozyme motif (e.g., '4qhp_zn_pentad') or a JSON array of reference ResidueCoords.
#[pg_extern(immutable, parallel_safe)]
pub fn theozyme_rmsd(site_atoms: Vec<ResidueCoord>, motif_name: &str) -> f64 {
    if site_atoms.is_empty() {
        return 999.0;
    }

    let template = if let Some(tpl) = get_builtin_theozyme_template(motif_name, site_atoms.len()) {
        tpl
    } else if let Ok(custom_tpl) = serde_json::from_str::<Vec<ResidueCoord>>(motif_name) {
        if custom_tpl.is_empty() {
            return 999.0;
        }
        custom_tpl
    } else {
        get_builtin_theozyme_template("4qhp_zn_pentad", 5).unwrap()
    };

    // If site_atoms contains more atoms than the template (e.g., full pocket vs 5-residue pentad),
    // filter catalytic sidechain/cofactor atoms or take the first template.len() motif atoms.
    if site_atoms.len() > template.len() {
        let filtered: Vec<ResidueCoord> = site_atoms
            .iter()
            .filter(|c| {
                let u = c.name.to_ascii_uppercase();
                u.contains("HIS") || u.contains("GLU") || u.contains("TYR") || u.contains("ASP")
                    || u == "H1" || u == "H2" || u == "E1" || u == "E2" || u == "Y"
            })
            .take(template.len())
            .cloned()
            .collect();
        if filtered.len() == template.len() {
            return kabsch_rmsd(filtered, template);
        }
    }

    kabsch_rmsd(site_atoms, template)
}

/// Scores two-sided enzyme-substrate interface complementarity in [0.0, 1.0] from pocket and peptide vector embeddings.
#[pg_extern(immutable, parallel_safe)]
pub fn interface_complementarity(pocket_embedding: Vec<f32>, peptide_embedding: Vec<f32>) -> f64 {
    let n = pocket_embedding.len().min(peptide_embedding.len());
    if n == 0 {
        return 0.0;
    }

    let mut dot = 0.0f64;
    let mut norm_p = 0.0f64;
    let mut norm_s = 0.0f64;

    for i in 0..n {
        let p = pocket_embedding[i] as f64;
        let s = peptide_embedding[i] as f64;
        dot += p * s;
        norm_p += p * p;
        norm_s += s * s;
    }

    if norm_p < 1e-12 || norm_s < 1e-12 {
        return 0.0;
    }

    let cos_sim = (dot / (norm_p.sqrt() * norm_s.sqrt())).clamp(-1.0, 1.0);
    // Map cosine similarity [-1.0, 1.0] to a calibrated complementarity probability [0.0, 1.0]
    ((cos_sim + 1.0) * 0.5).clamp(0.0, 1.0)
}

/// Evaluates beta-strand pairing complementarity (H-bond / C-alpha register and absence of steric clashes)
/// between the protease active-site guide strand and target substrate backbone coordinates.
#[pg_extern(immutable, parallel_safe)]
pub fn strand_pairing_check(
    guide_strand_coords: Vec<ResidueCoord>,
    backbone_coords: Vec<ResidueCoord>,
) -> bool {
    if guide_strand_coords.is_empty() || backbone_coords.is_empty() {
        return false;
    }

    let mut paired_contacts = 0usize;

    for g in &guide_strand_coords {
        let mut min_dist = f64::MAX;
        for b in &backbone_coords {
            let d = distance_angstroms(g.clone(), b.clone());
            if d < 2.0 {
                // Severe backbone steric clash prevents beta-strand pairing
                return false;
            }
            if d < min_dist {
                min_dist = d;
            }
        }
        // Canonical antiparallel/parallel beta-strand pairing window (2.4..=3.6 A for N..O H-bonds, up to 5.6 A for C-alpha register)
        if (2.4..=5.6).contains(&min_dist) {
            paired_contacts += 1;
        }
    }

    let required = guide_strand_coords.len().min(backbone_coords.len()).min(2).max(1);
    paired_contacts >= required
}

/// Computes a bipartite cross-reactivity interaction score across a SparseAttentionMap (CSR)
/// between enzyme active-site residue indices and substrate residue indices.
#[pg_extern(immutable, parallel_safe)]
pub fn sparse_cross_reactivity_score(
    map: SparseAttentionMap,
    enzyme_residues: Vec<i32>,
    substrate_residues: Vec<i32>,
) -> f64 {
    if map.entries.is_empty() || enzyme_residues.is_empty() || substrate_residues.is_empty() {
        return 0.0;
    }

    let mut total_weight = 0.0f64;
    let mut hits = 0usize;

    for entry in &map.entries {
        let forward = enzyme_residues.contains(&entry.source_residue)
            && substrate_residues.contains(&entry.target_residue);
        let reverse = enzyme_residues.contains(&entry.target_residue)
            && substrate_residues.contains(&entry.source_residue);
        if forward || reverse {
            total_weight += entry.weight as f64;
            hits += 1;
        }
    }

    if hits == 0 {
        0.0
    } else {
        total_weight / (hits as f64)
    }
}

/// Compact summary type storing a quantized 1D histogram and moments of a catalytic geometric metric
/// (e.g., proton-transfer distance or Burgi-Dunitz attack angle) across an MD or MLFF conformational ensemble.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct PreorganizationProfile {
    pub metric_name: String,
    pub bin_edges: Vec<f64>,
    pub bin_counts: Vec<f64>,
    pub mean_val: f64,
    pub std_dev: f64,
}

/// Constructs a PreorganizationProfile directly from pre-binned histogram edges and counts.
#[pg_extern(immutable, parallel_safe)]
pub fn create_preorganization_hist(
    metric_name: &str,
    bin_edges: Vec<f64>,
    bin_counts: Vec<f64>,
    mean_val: f64,
    std_dev: f64,
) -> PreorganizationProfile {
    PreorganizationProfile {
        metric_name: metric_name.to_string(),
        bin_edges,
        bin_counts,
        mean_val,
        std_dev,
    }
}

/// Quantizes raw trajectory observations (distances or angles across MD/MLFF frames) into a compact PreorganizationProfile.
#[pg_extern(immutable, parallel_safe)]
pub fn create_preorganization_profile(
    metric_name: &str,
    trajectory_values: Vec<f64>,
    min_val: f64,
    max_val: f64,
    num_bins: i32,
) -> PreorganizationProfile {
    let bins = num_bins.max(1) as usize;
    let lo = min_val.min(max_val);
    let hi = (max_val.max(min_val)).max(lo + 1e-6);
    let step = (hi - lo) / (bins as f64);

    let mut edges = Vec::with_capacity(bins + 1);
    for i in 0..=bins {
        edges.push(lo + (i as f64) * step);
    }

    let mut counts = vec![0.0f64; bins];
    if trajectory_values.is_empty() {
        return PreorganizationProfile {
            metric_name: metric_name.to_string(),
            bin_edges: edges,
            bin_counts: counts,
            mean_val: 0.0,
            std_dev: 0.0,
        };
    }

    let mut sum = 0.0f64;
    for &v in &trajectory_values {
        sum += v;
        let idx = if v <= lo {
            0
        } else if v >= hi {
            bins - 1
        } else {
            (((v - lo) / step).floor() as usize).min(bins - 1)
        };
        counts[idx] += 1.0;
    }

    let n = trajectory_values.len() as f64;
    let mean = sum / n;
    let var = trajectory_values
        .iter()
        .map(|&v| {
            let d = v - mean;
            d * d
        })
        .sum::<f64>()
        / n;

    PreorganizationProfile {
        metric_name: metric_name.to_string(),
        bin_edges: edges,
        bin_counts: counts,
        mean_val: mean,
        std_dev: var.sqrt(),
    }
}

/// Computes the fraction [0.0, 1.0] of an MD/MLFF conformational ensemble where the tracked catalytic metric stays <= `threshold`.
#[pg_extern(immutable, parallel_safe)]
pub fn ensemble_fraction_below(profile: PreorganizationProfile, threshold: f64) -> f64 {
    let n_bins = profile.bin_counts.len();
    if n_bins == 0 || profile.bin_edges.len() < n_bins + 1 {
        return 0.0;
    }

    let total: f64 = profile.bin_counts.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }

    let mut accum = 0.0f64;
    for i in 0..n_bins {
        let left = profile.bin_edges[i];
        let right = profile.bin_edges[i + 1];
        let count = profile.bin_counts[i];
        if threshold >= right {
            accum += count;
        } else if threshold > left && (right - left) > 1e-12 {
            let frac = (threshold - left) / (right - left);
            accum += count * frac;
        }
    }

    (accum / total).clamp(0.0, 1.0)
}

/// Computes the fraction [0.0, 1.0] of an MD/MLFF conformational ensemble where the tracked metric stays between `min_val` and `max_val`.
#[pg_extern(immutable, parallel_safe)]
pub fn ensemble_fraction_between(
    profile: PreorganizationProfile,
    min_val: f64,
    max_val: f64,
) -> f64 {
    if max_val <= min_val {
        return 0.0;
    }
    (ensemble_fraction_below(profile.clone(), max_val) - ensemble_fraction_below(profile, min_val))
        .clamp(0.0, 1.0)
}

// =====================================================================
// NEW: SEQUENCE ANALYSIS (GC Content & Translation)
// =====================================================================

#[pg_extern(immutable, parallel_safe)]
pub fn bio_gc_content(sequence: &str) -> f64 {
    bio::seq_analysis::gc::gc_content(sequence.as_bytes()) as f64
}

fn translate_codon(codon: &[u8]) -> char {
    match codon {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*', // Stop codon
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

#[pg_extern(immutable, parallel_safe)]
pub fn bio_translate_dna(sequence: &str) -> String {
    let mut aa = String::with_capacity(sequence.len() / 3);
    for chunk in sequence.as_bytes().chunks(3) {
        if chunk.len() == 3 {
            aa.push(translate_codon(chunk));
        }
    }
    aa
}

// =====================================================================
// 11. FREE-ENERGY DECOMPOSITION AUDITING & ALLOSTERIC HOTSPOT PROFILING (MM/GBSA & MM/PBSA)
// =====================================================================

/// Structured double-entry thermodynamic audit result comparing per-residue MM/GBSA or MM/PBSA
/// decomposition arrays against global endpoint energy terms (GGAS, EGB/EPB, ESURF, ENPOLAR+EDISPER).
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct EnergyDecompositionAudit {
    pub model: String,
    pub entity_count: i32,
    pub max_mm_polar_identity_error: f64,
    pub gas_closure_error: f64,
    pub polar_closure_error: f64,
    pub nonpolar_discrepancy: f64,
    pub total_mm_plus_polar: f64,
    pub gas_and_polar_closed: bool,
    pub nonpolar_closed: bool,
}

/// Audits per-residue MM/GBSA or MM/PBSA decomposition arrays in a single O(N) vectorized Rust pass.
/// Verifies:
///   1. Per-entity identity: `MM_plus_polar = Internal + vdW + Electrostatic + Polar Solvation`
///   2. Gas-phase closure: `SUM(Internal + vdW + Electrostatic) == global_ggas` within `tolerance`
///   3. Polar solvation closure: `SUM(Polar Solvation) == global_polar` within `tolerance`
///   4. Non-polar discrepancy: `global_nonpolar_net - SUM(Non-Polar Solv.)` (exposing LCPO vs. decomposed
///      GB surface area shifts and PB cavity/dispersion `ENPOLAR + EDISPER` unallocated penalties).
#[pg_extern(immutable, parallel_safe)]
pub fn audit_mmpbsa_decomposition(
    model: &str,
    internal: Vec<f64>,
    vdw: Vec<f64>,
    electrostatic: Vec<f64>,
    polar_solv: Vec<f64>,
    nonpolar_solv: Vec<f64>,
    reported_mm_plus_polar: Vec<f64>,
    global_ggas: f64,
    global_polar: f64,
    global_nonpolar_net: f64,
    tolerance: f64,
) -> EnergyDecompositionAudit {
    let n = internal
        .len()
        .min(vdw.len())
        .min(electrostatic.len())
        .min(polar_solv.len())
        .min(nonpolar_solv.len())
        .min(reported_mm_plus_polar.len());

    let mut sum_gas = 0.0f64;
    let mut sum_polar = 0.0f64;
    let mut sum_nonpolar = 0.0f64;
    let mut sum_mm_plus_polar = 0.0f64;
    let mut max_identity_err = 0.0f64;

    for i in 0..n {
        let gas_i = internal[i] + vdw[i] + electrostatic[i];
        let pol_i = polar_solv[i];
        let npol_i = nonpolar_solv[i];
        let calc_mm_pol = gas_i + pol_i;
        let err = (calc_mm_pol - reported_mm_plus_polar[i]).abs();
        if err > max_identity_err {
            max_identity_err = err;
        }
        sum_gas += gas_i;
        sum_polar += pol_i;
        sum_nonpolar += npol_i;
        sum_mm_plus_polar += reported_mm_plus_polar[i];
    }

    let gas_err = (sum_gas - global_ggas).abs();
    let polar_err = (sum_polar - global_polar).abs();
    let nonpol_disc = global_nonpolar_net - sum_nonpolar;
    let tol = tolerance.abs().max(1e-9);

    EnergyDecompositionAudit {
        model: model.trim().to_ascii_lowercase(),
        entity_count: n as i32,
        max_mm_polar_identity_error: max_identity_err,
        gas_closure_error: gas_err,
        polar_closure_error: polar_err,
        nonpolar_discrepancy: nonpol_disc,
        total_mm_plus_polar: sum_mm_plus_polar,
        gas_and_polar_closed: max_identity_err <= tol && gas_err <= tol && polar_err <= tol,
        nonpolar_closed: nonpol_disc.abs() <= tol,
    }
}

/// Computes the exact net solvation shift (in kcal/mol) between MM/PBSA and MM/GBSA endpoints:
/// `delta_solv = (pb_polar - gb_polar) + ((pb_enpolar + pb_edisper) - gb_nonpolar)`.
#[pg_extern(immutable, parallel_safe)]
pub fn solvent_shift_penalty(
    gb_polar: f64,
    gb_nonpolar: f64,
    pb_polar: f64,
    pb_enpolar: f64,
    pb_edisper: f64,
) -> f64 {
    let polar_shift = pb_polar - gb_polar;
    let nonpolar_shift = (pb_enpolar + pb_edisper) - gb_nonpolar;
    polar_shift + nonpolar_shift
}

/// Converts per-residue `MM_plus_polar` decomposition energies (in kcal/mol, where negative values
/// are stabilizing) into a sorted CSR `SparseAttentionMap` (`weight = -energy`) for all residues
/// meeting `mm_plus_polar <= energy_cutoff` (e.g., `-0.05` kcal/mol).
#[pg_extern(immutable, parallel_safe)]
pub fn build_energy_hotspot_map(
    sequence_length: i32,
    residue_ids: Vec<i32>,
    mm_plus_polar: Vec<f64>,
    ligand_node_id: i32,
    energy_cutoff: f64,
) -> SparseAttentionMap {
    let n = residue_ids.len().min(mm_plus_polar.len());
    let mut entries: Vec<AttentionEntry> = Vec::new();

    for i in 0..n {
        let e = mm_plus_polar[i];
        if e <= energy_cutoff && e < 0.0 {
            entries.push(AttentionEntry {
                source_residue: residue_ids[i],
                target_residue: ligand_node_id,
                weight: (-e) as f32,
            });
        }
    }

    entries.sort_by(|a, b| {
        b.weight
            .partial_cmp(&a.weight)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    SparseAttentionMap {
        sequence_length,
        entries,
    }
}

/// Traverses a CSR `SparseAttentionMap` around `ligand_node_id` and returns up to `limit`
/// peripheral/cryptic stabilizing hotspot residue IDs (`weight >= min_stabilization`) that lie
/// outside the pre-defined `known_pocket_residues` list.
#[pg_extern(immutable, parallel_safe)]
pub fn discover_cryptic_hotspots(
    map: SparseAttentionMap,
    ligand_node_id: i32,
    known_pocket_residues: Vec<i32>,
    min_stabilization: f64,
    limit: i32,
) -> Vec<i32> {
    if map.entries.is_empty() || limit <= 0 {
        return Vec::new();
    }
    let min_w = min_stabilization.max(0.0) as f32;

    let mut candidates: Vec<(i32, f32)> = map
        .entries
        .iter()
        .filter_map(|e| {
            let res_id = if e.target_residue == ligand_node_id {
                Some(e.source_residue)
            } else if e.source_residue == ligand_node_id {
                Some(e.target_residue)
            } else {
                None
            }?;
            if e.weight >= min_w && !known_pocket_residues.contains(&res_id) {
                Some((res_id, e.weight))
            } else {
                None
            }
        })
        .collect();

    candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    candidates
        .into_iter()
        .take(limit as usize)
        .map(|(res_id, _)| res_id)
        .collect()
}

/// Computes the normalized Polar/Ionic-to-Hydrophobic Allosteric Mechanism Index in [-1.0, +1.0]:
/// `(S_polar - S_hydrophobic) / (S_polar + S_hydrophobic)` across a CSR `SparseAttentionMap`.
/// Negative values (< 0) indicate hydrophobic clamping; positive values (> 0) indicate polar/ionic switching.
#[pg_extern(immutable, parallel_safe)]
pub fn allosteric_mechanism_index(
    map: SparseAttentionMap,
    ligand_node_id: i32,
    hydrophobic_residues: Vec<i32>,
    polar_ionic_residues: Vec<i32>,
) -> f64 {
    let target = vec![ligand_node_id];
    let s_hydro = sparse_cross_reactivity_score(map.clone(), hydrophobic_residues, target.clone());
    let s_polar = sparse_cross_reactivity_score(map, polar_ionic_residues, target);
    let denom = s_hydro + s_polar;
    if denom <= 1e-12 {
        0.0
    } else {
        ((s_polar - s_hydro) / denom).clamp(-1.0, 1.0)
    }
}

// =====================================================================
// 12. WET-LAB SYNTHETIC BIOLOGY: REVERSE TRANSLATION, CAI, EXPRESSION CASSETTES,
//     3D GLOBULAR COMPACTNESS GATE & ENZYME / BIOMINERALS DIGITAL TWIN
// =====================================================================

fn optimal_codon_and_weight(aa: char, host: &str) -> (&'static str, f64) {
    let h = host.trim().to_ascii_lowercase();
    let is_yeast = h.contains("cerevisiae") || h.contains("yeast") || h.contains("pichia") || h.contains("pastoris");
    let is_halophile = h.contains("halo") || h.contains("natrono") || h.contains("salin") || h.contains("archaea");
    let is_thermophile = h.contains("thermo") || h.contains("sulfolobus") || h.contains("pyro") || h.contains("acidianus");
    let is_human = h.contains("human") || h.contains("homo") || h.contains("mammalian") || h.contains("cho");
    let is_bacillus = h.contains("bacillus") || h.contains("subtilis") || h.contains("licheniformis");
    let is_pseudomonas = h.contains("pseudomonas") || h.contains("cupriavidus") || h.contains("burkholderia");

    match aa.to_ascii_uppercase() {
        'A' => {
            if is_yeast { ("GCT", 1.0) }
            else if is_halophile || is_pseudomonas { ("GCC", 1.0) }
            else if is_thermophile || is_bacillus { ("GCA", 1.0) }
            else if is_human { ("GCC", 1.0) }
            else { ("GCG", 1.0) } // E. coli K12 / BL21
        }
        'R' => {
            if is_yeast || is_thermophile { ("AGA", 1.0) }
            else if is_halophile || is_pseudomonas { ("CGC", 1.0) }
            else if is_human { ("CGG", 1.0) }
            else { ("CGC", 1.0) }
        }
        'N' => {
            if is_yeast || is_thermophile || is_bacillus { ("AAT", 0.96) }
            else { ("AAC", 1.0) }
        }
        'D' => {
            if is_halophile || is_pseudomonas || is_human { ("GAC", 1.0) }
            else { ("GAT", 1.0) }
        }
        'C' => {
            if is_yeast || is_thermophile { ("TGT", 1.0) }
            else { ("TGC", 1.0) }
        }
        'Q' => {
            if is_yeast || is_thermophile || is_bacillus { ("CAA", 1.0) }
            else { ("CAG", 1.0) }
        }
        'E' => {
            if is_halophile || is_pseudomonas || is_human { ("GAG", 1.0) }
            else { ("GAA", 1.0) }
        }
        'G' => {
            if is_yeast || is_thermophile || is_bacillus { ("GGT", 1.0) }
            else { ("GGC", 1.0) }
        }
        'H' => {
            if is_halophile || is_pseudomonas || is_human { ("CAC", 1.0) }
            else { ("CAT", 1.0) }
        }
        'I' => {
            if is_halophile || is_pseudomonas || is_human { ("ATC", 1.0) }
            else if is_thermophile { ("ATA", 0.98) }
            else { ("ATT", 1.0) }
        }
        'L' => {
            if is_yeast || is_thermophile { ("TTG", 1.0) }
            else if is_halophile { ("CTC", 1.0) }
            else { ("CTG", 1.0) }
        }
        'K' => {
            if is_halophile || is_pseudomonas || is_human { ("AAG", 1.0) }
            else { ("AAA", 1.0) }
        }
        'M' => ("ATG", 1.0),
        'F' => {
            if is_halophile || is_pseudomonas || is_human || is_yeast { ("TTC", 1.0) }
            else { ("TTT", 1.0) }
        }
        'P' => {
            if is_yeast || is_thermophile { ("CCA", 1.0) }
            else if is_halophile || is_human { ("CCC", 1.0) }
            else { ("CCG", 1.0) }
        }
        'S' => {
            if is_yeast || is_thermophile { ("TCT", 1.0) }
            else if is_halophile { ("TCG", 1.0) }
            else { ("AGC", 1.0) }
        }
        'T' => {
            if is_yeast || is_thermophile || is_bacillus { ("ACT", 1.0) }
            else if is_halophile { ("ACG", 1.0) }
            else { ("ACC", 1.0) }
        }
        'W' => ("TGG", 1.0),
        'Y' => {
            if is_halophile || is_pseudomonas || is_human { ("TAC", 1.0) }
            else { ("TAT", 1.0) }
        }
        'V' => {
            if is_yeast || is_thermophile { ("GTT", 1.0) }
            else if is_halophile { ("GTC", 1.0) }
            else { ("GTG", 1.0) }
        }
        '*' => {
            if is_human || is_halophile { ("TGA", 1.0) }
            else { ("TAA", 1.0) }
        }
        _ => ("NNN", 0.1),
    }
}

/// Reverse-translates an amino acid sequence into a host-optimized coding DNA sequence (CDS).
/// Supported host profiles include 'e_coli', 'b_subtilis', 'p_putida', 'c_necator',
/// 's_cerevisiae', 'p_pastoris', 'h_salinarum' (high-GC halophile), 's_acidocaldarius' (thermophile), and 'h_sapiens'.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_reverse_translate(protein_seq: &str, host_organism: &str) -> String {
    let mut dna = String::with_capacity(protein_seq.len() * 3);
    for ch in protein_seq.chars() {
        if ch.is_whitespace() {
            continue;
        }
        let (codon, _) = optimal_codon_and_weight(ch, host_organism);
        if codon != "NNN" || ch == '*' {
            dna.push_str(codon);
        }
    }
    dna
}

/// Computes the Codon Adaptation Index (CAI in [0.0, 1.0]) of a DNA coding sequence for a given expression host.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_codon_adaptation_index(dna_seq: &str, host_organism: &str) -> f64 {
    let clean: Vec<u8> = dna_seq
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let n_codons = clean.len() / 3;
    if n_codons == 0 {
        return 0.0;
    }

    let mut log_sum = 0.0f64;
    let mut counted = 0usize;

    for chunk in clean.chunks_exact(3) {
        let aa = translate_codon(chunk);
        if aa == '*' || aa == 'X' {
            continue;
        }
        let (opt_codon, _) = optimal_codon_and_weight(aa, host_organism);
        let codon_str = std::str::from_utf8(chunk).unwrap_or("NNN");
        let w: f64 = if codon_str == opt_codon {
            1.0
        } else if codon_str.as_bytes()[0..2] == opt_codon.as_bytes()[0..2] {
            0.78
        } else {
            0.42
        };
        log_sum += w.ln();
        counted += 1;
    }

    if counted == 0 {
        0.0
    } else {
        (log_sum / (counted as f64)).exp().clamp(0.0, 1.0)
    }
}

/// Structured wet-lab plasmid expression cassette returned by `bio_design_expression_cassette`.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct ExpressionCassette {
    pub host_organism: String,
    pub promoter_name: String,
    pub coding_dna: String,
    pub full_plasmid_dna: String,
    pub length_bp: i32,
    pub gc_content: f64,
    pub cai_score: f64,
    pub restriction_sites_free: bool,
}

/// Designs a complete wet-lab expression cassette (Promoter + RBS/Kozak + Host-Optimized CDS + Stop + Terminator)
/// and audits GC content, Codon Adaptation Index (CAI), and Golden Gate / Gibson restriction site safety.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_design_expression_cassette(
    protein_seq: &str,
    host_organism: &str,
    promoter_type: &str,
) -> ExpressionCassette {
    let cds = bio_reverse_translate(protein_seq, host_organism);
    let (stop_codon, _) = optimal_codon_and_weight('*', host_organism);

    let p_upper = promoter_type.trim().to_ascii_uppercase();
    let (promoter_seq, terminator_seq) = if p_upper.contains("GAL1") || p_upper.contains("AOX1") {
        (
            "ACGGATTAGAAGCCGCCGAGCGGGTGACAGCCCTCCGAAGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATGAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGCAAAAACTGCATAACCACTTTAACTAATACTTTCAACATTTTCGGTTTGTATTACTTCTTATTCAAATGTAATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATA",
            "GCGAATTTCTTATGATTTATGATTTTTATTATTAAATAAGTTATAAAAAAAATAAGTGTATACAAATTTTAAAGTGACTCTTAGTTTTAAAACGAAAATTCTTATTCTTGAGTAACTCTTTCCTGTAGGTCAGGTTGCTTTCTCAGGTATAGCATGAGGTCGCTCTTATTGACCACACCTCTACCGGCATG",
        )
    } else if p_upper.contains("KRT14") || p_upper.contains("CMV") {
        (
            "TAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTGGTTTAGTGAACCGTCAGATCCGCTAGCGCTACCGGTCGCCACC",
            "AATAAAAGATCTTTATTTTCATTAGATCTGTGTGTTGGTTTTTTGTGTG",
        )
    } else if p_upper.contains("PFDA") || p_upper.contains("HALO") {
        (
            "CTGTCAGCCGTTGAAGCCGGGGAAACCGCGGCGCTTTTATATAGGGCGCCGGCGCCGGAGGAGGTGAC",
            "CGGCCGGCCGCGCGCCTTTTTTTGCGCGC",
        )
    } else {
        // Default: T7 + strong Shine-Dalgarno RBS + spacer
        (
            "TAATACGACTCACTATAGGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT",
            "TAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        )
    };

    let full_plasmid = format!("{}{}{}{}", promoter_seq, cds, stop_codon, terminator_seq);
    let gc = bio_gc_content(&full_plasmid);
    let cai = bio_codon_adaptation_index(&cds, host_organism);

    // Check internal CDS freedom from BsaI (GGTCTC / GAGACC) Golden Gate type IIS sites
    let restriction_free = !cds.contains("GGTCTC") && !cds.contains("GAGACC");

    ExpressionCassette {
        host_organism: host_organism.to_string(),
        promoter_name: promoter_type.to_string(),
        coding_dna: cds,
        length_bp: full_plasmid.len() as i32,
        full_plasmid_dna: full_plasmid,
        gc_content: gc,
        cai_score: cai,
        restriction_sites_free: restriction_free,
    }
}

/// Structured 3D Globular Compactness & AlphaFold/ESMFold structural quality profile.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct GlobularFoldingProfile {
    pub ca_count: i32,
    pub mean_plddt: f64,
    pub radius_of_gyration: f64,
    pub compactness_ratio: f64,
    pub elongation: f64,
    pub is_globular: bool,
}

/// Evaluates the 3D Globular Compactness Gate directly from a PDB text string inside PostgreSQL:
///   1. Extracts CA coordinates and B-factor / pLDDT confidence scores (auto-scaling [0, 1] ESMFold to [0, 100])
///   2. Computes Radius of Gyration `Rg` (in Angstroms)
///   3. Computes Flory Globular Compactness Ratio `Rg / (2.25 * N^0.38)`
///   4. Computes principal-axis gyration tensor elongation via 3x3 Jacobi eigendecomposition
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_globular_compactness(pdb_text: &str) -> GlobularFoldingProfile {
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    let mut plddts: Vec<f64> = Vec::new();

    for line in pdb_text.lines() {
        if line.starts_with("ATOM  ") && line.len() >= 66 {
            let atom_name = line[12..16].trim();
            if atom_name == "CA" {
                if let (Ok(x), Ok(y), Ok(z), Ok(b)) = (
                    line[30..38].trim().parse::<f64>(),
                    line[38..46].trim().parse::<f64>(),
                    line[46..54].trim().parse::<f64>(),
                    line[60..66].trim().parse::<f64>(),
                ) {
                    coords.push((x, y, z));
                    plddts.push(b);
                }
            }
        }
    }

    let n = coords.len();
    if n < 10 {
        return GlobularFoldingProfile {
            ca_count: n as i32,
            mean_plddt: 0.0,
            radius_of_gyration: 99.9,
            compactness_ratio: 9.99,
            elongation: 9.99,
            is_globular: false,
        };
    }

    let nf = n as f64;
    let raw_mean_b: f64 = plddts.iter().sum::<f64>() / nf;
    // ESMFold local outputs store pLDDT in [0.0, 1.0], whereas AlphaFold2/v6 stores [0.0, 100.0]
    let mean_plddt = if raw_mean_b <= 1.05 && raw_mean_b > 0.0 {
        raw_mean_b * 100.0
    } else {
        raw_mean_b
    };

    let (mut cx, mut cy, mut cz) = (0.0, 0.0, 0.0);
    for &(x, y, z) in &coords {
        cx += x;
        cy += y;
        cz += z;
    }
    cx /= nf;
    cy /= nf;
    cz /= nf;

    let mut sum_r2 = 0.0f64;
    let mut cov = [[0.0f64; 3]; 3];
    for &(x, y, z) in &coords {
        let dx = x - cx;
        let dy = y - cy;
        let dz = z - cz;
        sum_r2 += dx * dx + dy * dy + dz * dz;
        let d = [dx, dy, dz];
        for r in 0..3 {
            for c in 0..3 {
                cov[r][c] += d[r] * d[c] / nf;
            }
        }
    }

    let rg = (sum_r2 / nf).sqrt();
    let ideal_rg = 2.25 * nf.powf(0.38);
    let compactness = if ideal_rg > 1e-9 { rg / ideal_rg } else { 9.99 };

    // 3x3 Jacobi eigendecomposition of gyration tensor `cov` for principal axis elongation
    for _sweep in 0..25 {
        let (mut p, mut q) = (0usize, 1usize);
        let mut max_off = 0.0f64;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if cov[i][j].abs() > max_off {
                    max_off = cov[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < 1e-12 {
            break;
        }
        let app = cov[p][p];
        let aqq = cov[q][q];
        let apq = cov[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        cov[p][p] = app - t * apq;
        cov[q][q] = aqq + t * apq;
        cov[p][q] = 0.0;
        cov[q][p] = 0.0;
        for k in 0..3 {
            if k != p && k != q {
                let ckp = cov[k][p];
                let ckq = cov[k][q];
                cov[k][p] = c * ckp - s * ckq;
                cov[p][k] = cov[k][p];
                cov[k][q] = s * ckp + c * ckq;
                cov[q][k] = cov[k][q];
            }
        }
    }

    let mut evals = [cov[0][0].max(0.0), cov[1][1].max(0.0), cov[2][2].max(0.0)];
    evals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let elongation = (evals[0] / evals[1].max(1e-4)).sqrt();

    let is_globular = mean_plddt >= 70.0 && compactness <= 1.55 && elongation <= 3.2;

    GlobularFoldingProfile {
        ca_count: n as i32,
        mean_plddt,
        radius_of_gyration: rg,
        compactness_ratio: compactness,
        elongation,
        is_globular,
    }
}

/// Digital Twin kinetic & biomineralization simulator for defluorinating enzymes (Fluoride ISE + Fluorapatite yield).
/// Computes:
///   - `fluoride_ppm`: Released free F- ion concentration (ppm) at `time_hours` under Michaelis-Menten + Hydrophobin surface adsorption boost
///   - `fluorapatite_mg`: Precipitated Ca5(PO4)3F bioceramic mass (mg/L) given stoichiometric Ca5:F1 capture by the AMELX mineralizing tag
///   - `specific_activity`: Effective catalytic turnover rate (umol F- / min / mg enzyme)
#[pg_extern(immutable, parallel_safe)]
pub fn simulate_defluorination_kinetics(
    kcat_per_hr: f64,
    km_u_m: f64,
    substrate_u_m: f64,
    hydrophobin_anchor_boost: f64,
    calcium_m_m: f64,
    time_hours: f64,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(time_hr, f64),
        name!(fluoride_ppm, f64),
        name!(fluorapatite_mg_per_l, f64),
        name!(fractional_conversion, f64),
    ),
> {
    let eff_sub = (substrate_u_m * hydrophobin_anchor_boost.max(0.1)).max(0.0);
    let mm_factor = if (km_u_m + eff_sub) > 1e-9 {
        eff_sub / (km_u_m + eff_sub)
    } else {
        0.0
    };
    let k_eff = (kcat_per_hr * mm_factor * 0.015).max(0.0);
    let frac = (1.0 - (-k_eff * time_hours.max(0.0)).exp()).clamp(0.0, 1.0);

    // 1 uM F- = 0.019 ppm; scale to nanoparticle surface fluoride capacity (max ~150 ppm for 1000 uM C-F eq)
    let max_f_ppm = substrate_u_m * 0.15;
    let f_ppm = max_f_ppm * frac;

    // Fluorapatite Ca5(PO4)3F molar mass = 504.3 g/mol vs F = 19.0 g/mol (mass ratio = 26.54x).
    // Limited stoichiometrically by available Ca2+ (5 mol Ca per 1 mol F; 1 mM Ca = 0.2 mM F capacity = 3.8 ppm F).
    let max_f_sequestered_ppm = (calcium_m_m.max(0.0) / 5.0) * 19.0;
    let f_captured_ppm = f_ppm.min(max_f_sequestered_ppm);
    let fap_mg_l = f_captured_ppm * (504.3 / 19.0);

    pgrx::iter::TableIterator::once((time_hours, f_ppm, fap_mg_l, frac))
}

#[cfg(any(test, feature = "pg_test"))]
#[pg_schema]
mod tests {
    use pgrx::prelude::*;
    use crate::{ResidueCoord, residue_z_index, embedding_cosine_distance};
    use crate::{AttentionEntry, SparseAttentionMap, get_top_interacting_residues};

    #[pg_test]
    fn test_parse_vcf_sql() {
        use std::io::Write;
        // create a temporary file
        let mut file = std::fs::File::create("/tmp/test_parse.vcf").unwrap();
        file.write_all(b"##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10000\trs123\tA\tT\t50\tPASS\tDP=100\n").unwrap();
        
        let result = Spi::get_one::<String>("SELECT chrom FROM parse_vcf('/tmp/test_parse.vcf');");
        assert_eq!(result.unwrap().unwrap(), "chr1");
    }

    #[pg_test]
    fn test_parse_pdb_sql() {
        use std::io::Write;
        let mut file = std::fs::File::create("/tmp/test_parse.pdb").unwrap();
        // A simple valid PDB atom record
        // ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
        file.write_all(b"ATOM      1  N   ALA A   1       1.234   5.678   9.012  1.00  0.00           N\n").unwrap();
        
        let result = Spi::get_one::<f64>("SELECT x FROM parse_pdb_file('/tmp/test_parse.pdb');");
        assert_eq!(result.unwrap().unwrap(), 1.234);
    }

    #[pg_test]
    fn test_distance() {
        let r1 = ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "ALA".to_string() };
        let r2 = ResidueCoord { x: 3.0, y: 4.0, z: 0.0, name: "GLY".to_string() };
        let d = crate::distance_angstroms(r1, r2);
        assert!((d - 5.0).abs() < 1e-6);
    }

    #[pg_test]
    fn test_z_order_encoding() {
        let r1 = ResidueCoord { x: 10.0, y: 10.0, z: 10.0, name: "TRP".to_string() };
        let r2 = ResidueCoord { x: 10.1, y: 10.1, z: 10.1, name: "PHE".to_string() };
        let z1 = residue_z_index(r1.clone());
        let z2 = residue_z_index(r2.clone());
        assert!(z1 > 0);
        assert!(z2 > 0);
        let r3 = ResidueCoord { x: 300.0, y: 10.0, z: 10.0, name: "LEU".to_string() };
        let z3 = residue_z_index(r3);
        assert!((z1 - z2).abs() < (z1 - z3).abs());
    }

    #[pg_test]
    fn test_embedding_cosine_distance() {
        let v1 = vec![1.0, 0.0, 0.0];
        let v2 = vec![1.0, 0.0, 0.0];
        let v3 = vec![0.0, 1.0, 0.0];
        assert!(embedding_cosine_distance(v1.clone(), v2.clone()).abs() < 1e-6);
        assert!((embedding_cosine_distance(v1.clone(), v3.clone()) - 1.0).abs() < 1e-6);
    }

    #[pg_test]
    fn test_embedding_uniqueness() {
        // Guarantee we never go back to naively implemented (e.g. relying on sequence length)
        let seq1 = "MFEGFERRLVD";
        let seq2 = "MFEGFERRLVA"; // Same length, 1 mutation
        let seq3 = "MFEGFERRLVDAAA"; // Different length
        
        let emb1 = crate::get_esm_embedding(seq1);
        let emb2 = crate::get_esm_embedding(seq2);
        let emb3 = crate::get_esm_embedding(seq3);
        
        let dist1_2 = crate::embedding_cosine_distance(emb1.clone(), emb2.clone());
        let dist1_3 = crate::embedding_cosine_distance(emb1.clone(), emb3.clone());
        
        // Distances must be strictly greater than 0, proving uniqueness
        assert!(dist1_2 > 0.0001, "Same length sequences must not have 0 distance!");
        assert!(dist1_3 > 0.0001, "Different sequences must not have 0 distance!");
        
        // A single mutation should be relatively close compared to radically different sequences
        // (In tri-peptide feature hashing, one mutation changes up to 3 tri-peptides)
        assert!(dist1_2 < 0.5, "1 mutation should still remain reasonably close");
    }


    #[pg_test]
    fn test_sparse_attention() {
        let map = SparseAttentionMap {
            sequence_length: 400,
            entries: vec![
                AttentionEntry { source_residue: 10, target_residue: 45, weight: 0.8 },
                AttentionEntry { source_residue: 12, target_residue: 45, weight: 0.95 }, // Highest
                AttentionEntry { source_residue: 200, target_residue: 45, weight: 0.2 },
            ],
        };
        
        // Find top 2 residues interacting with residue 45
        let top = get_top_interacting_residues(map, 45, 2);
        assert_eq!(top.len(), 2);
        assert_eq!(top[0], 12); // Highest weight
        assert_eq!(top[1], 10); // Second highest
    }

    #[pg_test]
    fn test_molecular_weight() {
        let w = crate::molecular_weight("A");
        assert!((w - 89.094).abs() < 1e-3);
        let w2 = crate::molecular_weight("ARND");
        assert!(w2 > 300.0);
    }

    #[pg_test]
    fn test_center_of_mass() {
        let atoms = vec![
            crate::ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "C".to_string() },
            crate::ResidueCoord { x: 10.0, y: 10.0, z: 10.0, name: "C".to_string() },
        ];
        let com = crate::center_of_mass(atoms);
        assert!((com.x - 5.0).abs() < 1e-6);
        assert!((com.y - 5.0).abs() < 1e-6);
        assert!((com.z - 5.0).abs() < 1e-6);
    }

    #[pg_test]
    fn test_prosite_match() {
        let pattern = "[ST]-x(2)-[RK]";
        assert!(crate::prosite_match("ASAAR", pattern) == true); 
        assert!(crate::prosite_match("ATAAR", pattern) == true); 
        assert!(crate::prosite_match("ASAAAR", pattern) == false); // 3 chars instead of 2
        
        let pattern2 = "{ST}-x-[RK]";
        assert!(crate::prosite_match("AAAR", pattern2) == true); 
        assert!(crate::prosite_match("ASAR", pattern2) == false); 
    }

    #[pg_test]
    fn test_dihedral_angle() {
        let p1 = crate::ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "A".to_string() };
        let p2 = crate::ResidueCoord { x: 1.0, y: 0.0, z: 0.0, name: "B".to_string() };
        let p3 = crate::ResidueCoord { x: 1.0, y: 1.0, z: 0.0, name: "C".to_string() };
        let p4 = crate::ResidueCoord { x: 1.0, y: 1.0, z: 1.0, name: "D".to_string() };
        
        let angle = crate::dihedral_angle(p1, p2, p3, p4);
        assert!((angle - 90.0).abs() < 1e-6);
    }

    #[pg_test]
    fn test_is_steric_clash() {
        let p1 = crate::ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "A".to_string() };
        let p2 = crate::ResidueCoord { x: 1.5, y: 0.0, z: 0.0, name: "B".to_string() };
        
        // Sum of VdW radii = 2.0. Distance = 1.5. So they clash.
        assert!(crate::is_steric_clash(p1.clone(), p2.clone(), 2.0));
        // Sum of VdW radii = 1.0. Distance = 1.5. So they don't clash.
        assert!(!crate::is_steric_clash(p1.clone(), p2.clone(), 1.0));
    }

    #[pg_test]
    fn test_parse_smiles() {
        // Aspirin SMILES
        let valid_mol = crate::parse_smiles("CC(=O)OC1=CC=CC=C1C(=O)O");
        assert!(valid_mol.is_valid);
        assert!(crate::smiles_contains(valid_mol.clone(), "C1=CC=CC=C1")); // Benzene ring

        // Invalid SMILES
        let invalid_mol = crate::parse_smiles("C(C"); // Unclosed parenthesis
        assert!(!invalid_mol.is_valid);
        assert!(!crate::smiles_contains(invalid_mol, "C"));
    }

    #[pg_test]
    fn test_sequence_entropy() {
        let poly_g = "GGGGGGGGGG";
        let entropy_g = crate::sequence_entropy(poly_g);
        assert!((entropy_g - 0.0).abs() < 1e-6); // entropy of 1 char is 0

        let complex = "ACDEFGHIKLMNPQRSTVWY";
        let entropy_complex = crate::sequence_entropy(complex);
        assert!(entropy_complex > 4.0); // very high entropy
    }

    #[pg_test]
    fn test_gravy_score() {
        let g1 = crate::gravy_score("IVLF"); // Highly hydrophobic
        assert!(g1 > 3.0);
        let g2 = crate::gravy_score("RKDE"); // Highly hydrophilic
        assert!(g2 < -3.0);
    }

    #[pg_test]
    fn test_isoelectric_point() {
        // Simple peptide with known pI
        // Poly-Lysine should be very basic (high pI)
        let pi1 = crate::isoelectric_point("KKKKK");
        assert!(pi1 > 9.0);
        
        // Poly-Glutamic acid should be very acidic (low pI)
        let pi2 = crate::isoelectric_point("EEEEE");
        assert!(pi2 < 5.0);
    }


    #[pg_test]
    fn test_smith_waterman() {
        // Perfect match
        assert_eq!(crate::sequence_alignment_score("ACDEF", "ACDEF"), 15); // 5 * 3
        
        // Partial match
        // ACDEF
        // AC-EF
        assert_eq!(crate::sequence_alignment_score("ACDEF", "ACEF"), 10); // AC=6, EF=6, gap=-2 = 10
        
        // No match
        assert_eq!(crate::sequence_alignment_score("AAAAA", "CCCCC"), 0);
    }

    #[pg_test(error = "HTTP request to UniProt failed: http status: 400")]
    fn test_bio_fetch_uniprot_invalid() {
        // Should panic with HTTP error
        crate::bio_fetch_uniprot("INVALID_ID_123");
    }



    #[pg_test]
    fn test_uniprot_search_gene() {
        // Search for human BRCA1 (gene:BRCA1 AND taxonomy_id:9606 AND reviewed:true)
        let results: Vec<_> = crate::bio_search_uniprot("gene:BRCA1 AND taxonomy_id:9606 AND reviewed:true").collect();
        assert!(results.len() > 0);
        assert_eq!(results[0].0, "P38398"); // Human BRCA1 accession
    }

    #[pg_test]
    fn test_uniprot_search_length_filter() {
        // Search for very short human proteins
        let results: Vec<_> = crate::bio_search_uniprot("length:[1 TO 50] AND taxonomy_id:9606 AND reviewed:true").collect();
        assert!(results.len() > 0);
        assert!(results[0].3.len() <= 50);
    }


    #[pg_test]
    fn test_bio_fold_sequence() {
        // Fold a very short micropeptide (6 amino acids) so the API returns instantly
        let pdb = crate::bio_fold_sequence("MKVLYL");
        assert!(pdb.contains("ATOM"));
    }

    #[pg_test]
    fn test_bio_search_uniprot() {
        let results: Vec<_> = crate::bio_search_uniprot("taxonomy_id:2594042").collect();
        assert!(results.len() > 0);
        let first = &results[0];
        
        assert_eq!(first.1, "Promethearchaeum syntrophicum"); // organism
    }


    #[pg_test]
    fn test_smiles_to_fingerprint() {
        let fp1 = crate::smiles_to_fingerprint("C1=CC=CC=C1"); // Benzene
        let fp2 = crate::smiles_to_fingerprint("CC1=CC=CC=C1"); // Toluene
        assert_eq!(fp1.len(), 1024);
        assert_eq!(fp2.len(), 1024);
        
        let sim = crate::tanimoto_similarity(fp1.clone(), fp2.clone());
        assert!(sim > 0.0 && sim < 1.0); // Should be somewhat similar but not identical
        
        let sim_self = crate::tanimoto_similarity(fp1.clone(), fp1.clone());
        assert!((sim_self - 1.0).abs() < 1e-6);
    }

    #[pg_test]
    fn test_tanimoto_sql() {
        // We use SPI to test if the custom operator `%` was properly registered
        // and returns the expected tanimoto calculation
        let sim = Spi::get_one::<f64>(
            "SELECT smiles_to_fingerprint('C1=CC=CC=C1') % smiles_to_fingerprint('CC1=CC=CC=C1');"
        ).expect("SPI failed to execute tanimoto % operator");
        
        let sim = sim.expect("Tanimoto similarity should not be NULL");
        assert!(sim > 0.0 && sim < 1.0);
    }

    #[pg_test]
    fn test_spectrum_creation() {
        let mz = vec![100.0, 200.0];
        let intensity = vec![10.0, 50.0];
        let spec = crate::create_spectrum(mz.clone(), intensity.clone());
        assert_eq!(spec.m_z.len(), 2);
        assert_eq!(spec.intensity.len(), 2);
    }

    #[pg_test]
    fn test_ppm_difference() {
        let mass1 = 1000.0;
        let mass2 = 1000.0032; // 3.2 ppm diff
        let diff = crate::ppm_difference(mass1, mass2);
        assert!((diff - 3.2).abs() < 1e-4);
        
        assert!(crate::is_mass_within_ppm(mass1, mass2, 3.5));
        assert!(!crate::is_mass_within_ppm(mass1, mass2, 3.0));
    }

    #[pg_test]
    fn test_theoretical_fragments() {
        let frags_hcd: Vec<_> = crate::generate_theoretical_fragments("ACD", 1, "HCD").collect();
        assert_eq!(frags_hcd.len(), 4);
        
        let b1 = frags_hcd.iter().find(|f| f.0 == "b1").unwrap();
        assert!((b1.1 - 72.044386).abs() < 1e-4);
        assert_eq!(b1.2, 1);
        
        let y1 = frags_hcd.iter().find(|f| f.0 == "y1").unwrap();
        assert!((y1.1 - 134.044781).abs() < 1e-4);
        assert_eq!(y1.2, 1);
        
        let frags_etd: Vec<_> = crate::generate_theoretical_fragments("ACD", 1, "ETD").collect();
        assert_eq!(frags_etd.len(), 4);
        
        let c1 = frags_etd.iter().find(|f| f.0 == "c1").unwrap();
        assert!((c1.1 - 89.070935).abs() < 1e-4);
    }

    #[pg_test]
    fn test_isotope_fit_score() {
        let charge = 2;
        let theoretical_mz = 1000.0;
        let ppm_tol = 10.0;
        
        // Exact match with expected averagine envelope for mass ~2000
        // num_carbons = (1998 / 111.125) * 4.9384 = ~88.8 carbons
        // m1 = 88.8 * 0.0107 = ~0.95
        // m2 = 0.95^2 / 2 = ~0.45
        // total = 1 + 0.95 + 0.45 = 2.4
        // env: M=0.41, M+1=0.39, M+2=0.19
        
        let observed_perfect = crate::create_spectrum(
            vec![1000.0, 1000.0 + (1.00335/2.0), 1000.0 + (2.0067/2.0)],
            vec![41.0, 39.0, 19.0]
        );
        let score1 = crate::isotope_fit_score(observed_perfect, theoretical_mz, charge, ppm_tol);
        assert!(score1 > 0.98); // High cosine similarity expected
        
        // Poor match (only M peak found)
        let observed_poor = crate::create_spectrum(
            vec![1000.0],
            vec![100.0]
        );
        let score2 = crate::isotope_fit_score(observed_poor, theoretical_mz, charge, ppm_tol);
        assert!(score2 > 0.0 && score2 < 0.8); // Should be much lower

        // No match
        let observed_none = crate::create_spectrum(
            vec![1200.0],
            vec![100.0]
        );
        let score3 = crate::isotope_fit_score(observed_none, theoretical_mz, charge, ppm_tol);
        assert_eq!(score3, 0.0);
    }

    #[pg_test]
    fn test_decoy_sequence() {
        let seq = "ACDEFGH";
        
        // Reverse (non-tryptic C-terminus reverses full sequence)
        let rev = crate::generate_decoy_sequence(seq, "reverse");
        assert_eq!(rev, "HGFEDCA");
        
        // Tryptic peptide ending in R or K preserves C-terminal cleavage site (Elias & Gygi pseudo-reversal)
        let tryptic = "SIYFLLTDR";
        let rev_tryptic = crate::generate_decoy_sequence(tryptic, "reverse");
        assert_eq!(rev_tryptic, "DTLLFYISR");
        assert_eq!(crate::generate_decoy_sequence(tryptic, "full_reverse"), "RDTLLFYIS");

        // Verify y1 ion (m/z 175.1190 for R) is conserved between tryptic target and decoy
        let target_frags: Vec<_> = crate::generate_theoretical_fragments(tryptic, 1, "HCD").collect();
        let decoy_frags: Vec<_> = crate::generate_theoretical_fragments(&rev_tryptic, 1, "HCD").collect();
        let target_y1 = target_frags.iter().find(|f| f.0 == "y1").unwrap().1;
        let decoy_y1 = decoy_frags.iter().find(|f| f.0 == "y1").unwrap().1;
        assert!((target_y1 - decoy_y1).abs() < 1e-6);

        // Shuffle
        let shuff = crate::generate_decoy_sequence(seq, "shuffle");
        assert_eq!(shuff.len(), seq.len());
        assert_ne!(shuff, seq);
        // Ensure same characters
        let mut sorted_seq: Vec<char> = seq.chars().collect();
        let mut sorted_shuff: Vec<char> = shuff.chars().collect();
        sorted_seq.sort();
        sorted_shuff.sort();
        assert_eq!(sorted_seq, sorted_shuff);

        // Tryptic shuffle preserves C-terminal K/R
        let shuff_tryptic = crate::generate_decoy_sequence(tryptic, "shuffle");
        assert!(shuff_tryptic.ends_with('R'));
        assert_ne!(shuff_tryptic, tryptic);
        
        // Random
        let rand_seq = crate::generate_decoy_sequence(seq, "random");
        assert_eq!(rand_seq.len(), seq.len());
        assert_ne!(rand_seq, seq);
    }

    #[pg_test]
    fn test_parse_mgf() {
        use std::io::Write;
        let mut file = std::fs::File::create("/tmp/test_parse.mgf").unwrap();
        file.write_all(b"BEGIN IONS\nTITLE=Test Spectrum\nPEPMASS=1000.5\nCHARGE=2+\n100.0 5000\n200.0 10000\nEND IONS\n").unwrap();
        
        let result = pgrx::Spi::get_one::<String>("SELECT title FROM parse_mgf('/tmp/test_parse.mgf');");
        assert_eq!(result.unwrap().unwrap(), "Test Spectrum");
        
        let pepmass = pgrx::Spi::get_one::<f64>("SELECT pepmass FROM parse_mgf('/tmp/test_parse.mgf');");
        assert_eq!(pepmass.unwrap().unwrap(), 1000.5);
    }

    #[pg_test]
    fn test_bio_gc_content() {
        let gc = crate::bio_gc_content("ATGCGC");
        assert!((gc - 0.666666).abs() < 1e-4);
    }

    #[pg_test]
    fn test_bio_translate_dna() {
        let aa = crate::bio_translate_dna("ATGCGT");
        assert_eq!(aa, "MR");
        let aa_stop = crate::bio_translate_dna("TAATAG");
        assert_eq!(aa_stop, "**");
    }

    #[pg_test]
    fn test_bio_rmsd() {
        let c1 = crate::ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "CA".to_string() };
        let c2 = crate::ResidueCoord { x: 3.0, y: 4.0, z: 0.0, name: "CA".to_string() };
        let rmsd = crate::bio_rmsd(vec![c1.clone()], vec![c2.clone()]);
        assert_eq!(rmsd, 5.0);
    }

    #[pg_test]
    fn test_bio_parse_pdb() {
        let pdb_text = "ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00 49.33           N  \nATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 49.19           C  ";
        let iter = crate::bio_parse_pdb(pdb_text);
        let items: Vec<_> = iter.collect();
        assert_eq!(items.len(), 2);
        assert_eq!(items[0].0, "N");
        assert_eq!(items[0].1, "MET");
        assert_eq!(items[1].0, "CA");
    }

    #[pg_test]
    fn test_hetatm_and_coordinated_water_shell() {
        let pdb_with_hetatm = "\
ATOM      1  NE2 HIS A  92      12.420  18.150  24.080  1.00 20.00           N  \n\
ATOM      2  NE2 HIS A  96      15.610  17.890  23.620  1.00 20.00           N  \n\
ATOM      3  OE1 GLU A 102      13.940  19.680  21.140  1.00 20.00           O  \n\
HETATM    4 ZN    ZN A 201      14.020  17.910  22.980  1.00 15.00          ZN  \n\
HETATM    5  O   HOH A 301      14.110  16.050  23.820  1.00 18.00           O  \n";

        let coords = crate::bio_parse_pdb_coords(pdb_with_hetatm);
        assert_eq!(coords.len(), 5);
        assert!(coords.iter().any(|c| c.name == "ZN"));
        assert!(coords.iter().any(|c| c.name == "HOH:O"));

        // Verify positive detection of coordinated catalytic water shell within 2.4 A of ZN
        assert!(crate::has_coordinated_water_shell(coords.clone(), "ZN", 2.4));

        // Remove catalytic water -> should return false
        let no_water: Vec<_> = coords.into_iter().filter(|c| !c.name.starts_with("HOH")).collect();
        assert!(!crate::has_coordinated_water_shell(no_water, "ZN", 2.4));
    }

    #[pg_test]
    fn test_bio_angle_and_bio_dihedral_and_bio_distance() {
        // Burgi-Dunitz nucleophilic attack geometry (~107 degrees)
        let substrate_c = crate::create_residue_coord(0.0, 0.0, 0.0, "SUB:C");
        let substrate_o = crate::create_residue_coord(1.23, 0.0, 0.0, "SUB:O");
        // Place nucleophilic water O at 2.1 A distance and ~107 deg angle relative to C=O bond
        let theta_rad = 107.0_f64.to_radians();
        let water_o = crate::create_residue_coord(2.1 * theta_rad.cos(), 2.1 * theta_rad.sin(), 0.0, "HOH:O");

        let angle = crate::bio_angle(water_o.clone(), substrate_c.clone(), substrate_o.clone());
        assert!((angle - 107.0).abs() < 1e-4);
        assert!(angle >= 100.0 && angle <= 115.0);

        let dist = crate::bio_distance(water_o.clone(), substrate_c.clone());
        assert!((dist - 2.1).abs() < 1e-6);

        let p4 = crate::create_residue_coord(1.23, 0.0, 1.0, "SUB:N");
        let dih = crate::bio_dihedral(water_o, substrate_c, substrate_o, p4);
        assert!((dih.abs() - 90.0).abs() < 1e-4);
    }

    #[pg_test]
    fn test_kabsch_and_theozyme_rmsd() {
        // Get the 5-residue 4qhp_zn_pentad reference and apply a rigid-body 90-degree rotation + translation
        let ref_pentad = vec![
            crate::create_residue_coord(12.420, 18.150, 24.080, "HIS1:NE2"),
            crate::create_residue_coord(15.610, 17.890, 23.620, "HIS2:NE2"),
            crate::create_residue_coord(13.850, 14.520, 22.410, "GLU1:OE2"),
            crate::create_residue_coord(13.940, 19.680, 21.140, "GLU2:OE1"),
            crate::create_residue_coord(16.480, 15.240, 25.190, "TYR:OH"),
        ];

        // Rotate 90 deg around Z-axis: (x, y, z) -> (-y + 50, x - 30, z + 15)
        let rotated_translated: Vec<_> = ref_pentad
            .iter()
            .map(|c| crate::create_residue_coord(-c.y + 50.0, c.x - 30.0, c.z + 15.0, &c.name))
            .collect();

        // Naive RMSD should be huge (>40 A), whereas Kabsch/theozyme RMSD must be ~0.0 A!
        let naive = crate::bio_rmsd(rotated_translated.clone(), ref_pentad.clone());
        assert!(naive > 40.0);

        let k_rmsd = crate::kabsch_rmsd(rotated_translated.clone(), ref_pentad);
        assert!(k_rmsd < 1e-5, "Kabsch RMSD under pure rigid-body transformation should be ~0, got {}", k_rmsd);

        let t_rmsd = crate::theozyme_rmsd(rotated_translated, "4qhp_zn_pentad");
        assert!(t_rmsd < 1e-5);
    }

    #[pg_test]
    fn test_interface_complementarity_and_strand_pairing() {
        let pocket_emb = vec![0.8_f32, 0.6_f32, 0.0_f32];
        let cognate_emb = vec![0.8_f32, 0.6_f32, 0.0_f32];
        let orthogonal_emb = vec![-0.8_f32, -0.6_f32, 0.0_f32];

        let score_cognate = crate::interface_complementarity(pocket_emb.clone(), cognate_emb);
        let score_ortho = crate::interface_complementarity(pocket_emb, orthogonal_emb);
        assert!((score_cognate - 1.0).abs() < 1e-5);
        assert!(score_ortho < 0.01);

        // Beta-strand pairing check (2.9 A H-bond register -> true; 1.2 A steric clash -> false)
        let guide = vec![
            crate::create_residue_coord(0.0, 0.0, 0.0, "GLY:N"),
            crate::create_residue_coord(3.5, 0.0, 0.0, "ALA:O"),
        ];
        let substrate_good = vec![
            crate::create_residue_coord(0.0, 2.9, 0.0, "SUB1:O"),
            crate::create_residue_coord(3.5, 3.0, 0.0, "SUB2:N"),
        ];
        let substrate_clash = vec![
            crate::create_residue_coord(0.0, 1.1, 0.0, "SUB1:O"),
            crate::create_residue_coord(3.5, 3.0, 0.0, "SUB2:N"),
        ];
        assert!(crate::strand_pairing_check(guide.clone(), substrate_good));
        assert!(!crate::strand_pairing_check(guide, substrate_clash));

        let map = crate::create_sparse_map(200, vec![92, 93, 149], vec![185, 186, 187], vec![0.92, 0.88, 0.90]);
        let csr_score = crate::sparse_cross_reactivity_score(map, vec![92, 93, 149], vec![185, 186, 187]);
        assert!((csr_score - 0.90).abs() < 1e-4);
    }

    #[pg_test]
    fn test_preorganization_profile_and_ensemble_fractions() {
        // 10 MD frames: 9 frames between 2.4 and 2.9 A (catalytic), 1 frame at 3.8 A
        let frames = vec![2.45, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 3.80];
        let prof = crate::create_preorganization_profile("proton_transfer_dist", frames, 2.0, 4.0, 4);
        assert_eq!(prof.bin_counts.len(), 4);
        assert!(prof.mean_val > 2.5 && prof.mean_val < 3.0);

        let frac_below_3 = crate::ensemble_fraction_below(prof.clone(), 3.0);
        assert!((frac_below_3 - 0.90).abs() < 1e-5);
        assert!(frac_below_3 >= 0.80);

        let frac_between = crate::ensemble_fraction_between(prof, 2.0, 3.0);
        assert!((frac_between - 0.90).abs() < 1e-5);

        let hist = crate::create_preorganization_hist(
            "burgi_dunitz_angle",
            vec![90.0, 100.0, 115.0, 130.0],
            vec![5.0, 85.0, 10.0],
            107.2,
            4.1,
        );
        let angle_frac = crate::ensemble_fraction_between(hist, 100.0, 115.0);
        assert!((angle_frac - 0.85).abs() < 1e-5);
    }

    #[pg_test]
    fn test_mmpbsa_decomposition_audit_and_cryptic_hotspots() {
        // 1. Test audit_mmpbsa_decomposition on 5-LOX PB-like decomposition
        //    (gas & polar close exactly, whereas PB residue nonpolar is 0.0 vs +23.91 global ENPOLAR+EDISPER)
        let internal = vec![0.0, 0.0, 0.0, 0.0];
        let vdw = vec![-1.50, -1.60, -1.00, -37.28];
        let elec = vec![-0.20, -0.30, -0.10, 0.60];
        let polar = vec![0.51, 0.29, 0.10, 13.09];
        let nonpolar = vec![0.0, 0.0, 0.0, 0.0];
        let mm_plus_polar = vec![-1.19, -1.61, -1.00, -23.59];

        let audit_pb = crate::audit_mmpbsa_decomposition(
            "pb",
            internal,
            vdw,
            elec,
            polar,
            nonpolar,
            mm_plus_polar.clone(),
            -41.38,
            13.99,
            23.91,
            0.005,
        );
        assert_eq!(audit_pb.model, "pb");
        assert_eq!(audit_pb.entity_count, 4);
        assert!(audit_pb.max_mm_polar_identity_error < 1e-9);
        assert!(audit_pb.gas_closure_error < 1e-6);
        assert!(audit_pb.polar_closure_error < 1e-6);
        assert!(audit_pb.gas_and_polar_closed);
        assert!(!audit_pb.nonpolar_closed);
        assert!((audit_pb.nonpolar_discrepancy - 23.91).abs() < 1e-6);

        // 2. Test solvent_shift_penalty (α-Amyrin: GB polar=10.94, GB nonpol=-4.79, PB polar=13.99, ENPOLAR=-28.07, EDISPER=51.98)
        let shift = crate::solvent_shift_penalty(10.94, -4.79, 13.99, -28.07, 51.98);
        assert!((shift - 31.75).abs() < 1e-4);

        // 3. Test build_energy_hotspot_map, discover_cryptic_hotspots, and allosteric_mechanism_index
        //    Residues: 110 (Val110, pocket), 130 (His130, pocket), 126 (Ile126, cryptic!), 29 (Val29, cryptic!), 108 (Glu108, repulsive)
        let res_ids = vec![68, 108, 110, 130, 133, 134, 126, 29];
        let energies = vec![-0.40, 0.35, -1.48, -1.57, -1.45, -1.49, -1.05, -0.85];
        let map = crate::build_energy_hotspot_map(677, res_ids, energies, 999, -0.05);
        assert_eq!(map.entries.len(), 7); // Glu108 (+0.35) excluded

        let known_pocket = vec![68, 108, 109, 110, 129, 130, 133, 134, 137];
        let cryptic = crate::discover_cryptic_hotspots(map.clone(), 999, known_pocket, 0.50, 5);
        assert_eq!(cryptic, vec![126, 29]);

        // Hydrophobic vs Polar/Ionic Mechanism Index
        let ami = crate::allosteric_mechanism_index(map, 999, vec![110, 130], vec![133, 134]);
        assert!(ami > -0.15 && ami < 0.15);
    }

    #[pg_test]
    fn test_wetlab_cassette_globular_and_defluorination() {
        // 1. Test multi-host reverse translation & exact round-trip via bio_translate_dna
        let protein = "MFEGFERRLVDWPSTDKTKREEVD";
        for host in &["e_coli", "s_cerevisiae", "h_salinarum", "s_acidocaldarius", "h_sapiens", "b_subtilis"] {
            let cds = crate::bio_reverse_translate(protein, host);
            assert_eq!(cds.len(), protein.len() * 3);
            let roundtrip = crate::bio_translate_dna(&cds);
            assert_eq!(roundtrip, protein, "Round-trip translation failed for host {}", host);

            let cai = crate::bio_codon_adaptation_index(&cds, host);
            assert!((cai - 1.0).abs() < 0.05, "Optimal codons for {} should have CAI ~ 1.0, got {}", host, cai);
        }

        // Verify halophile vs yeast GC bias difference
        let cds_halo = crate::bio_reverse_translate(protein, "h_salinarum");
        let cds_yeast = crate::bio_reverse_translate(protein, "s_cerevisiae");
        assert!(crate::bio_gc_content(&cds_halo) > crate::bio_gc_content(&cds_yeast) + 0.15);

        // 2. Test bio_design_expression_cassette
        let cassette = crate::bio_design_expression_cassette(protein, "e_coli", "T7");
        assert_eq!(cassette.host_organism, "e_coli");
        assert!(cassette.length_bp > (protein.len() as i32) * 3 + 100);
        assert!(cassette.cai_score > 0.98);
        assert!(cassette.restriction_sites_free);

        // 3. Test q1ju72_dehalogenase_hexad theozyme_rmsd
        let q1ju72_site = vec![
            crate::create_residue_coord(1.238, 2.722, -3.192, "ASP104:OD2"),
            crate::create_residue_coord(2.345, -0.494, -2.097, "HIS271:NE2"),
            crate::create_residue_coord(4.659, -0.482, 2.204, "ASP128:OD2"),
            crate::create_residue_coord(5.305, 5.974, -5.207, "TRP150:NE1"),
            crate::create_residue_coord(2.729, 3.340, -7.649, "TYR212:OH"),
            crate::create_residue_coord(2.392, 6.631, -1.154, "ARG108:NH2"),
        ];
        let rmsd_deha = crate::theozyme_rmsd(q1ju72_site, "q1ju72_dehalogenase_hexad");
        assert!(rmsd_deha < 1e-6);

        // 4. Test evaluate_globular_compactness on compact 3D sphere-like CA trace
        let mut pdb_lines = String::new();
        for i in 0..30 {
            let theta = (i as f64) * 0.45;
            let phi = (i as f64) * 0.25;
            let r = 8.0;
            let x = r * theta.cos() * phi.sin();
            let y = r * theta.sin() * phi.sin();
            let z = r * phi.cos();
            pdb_lines.push_str(&format!(
                "ATOM  {:5}  CA  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00 94.50           C  \n",
                i + 1,
                i + 1,
                x,
                y,
                z
            ));
        }
        let glob = crate::evaluate_globular_compactness(&pdb_lines);
        assert_eq!(glob.ca_count, 30);
        assert!((glob.mean_plddt - 94.5).abs() < 1e-3);
        assert!(glob.is_globular);

        // 5. Test simulate_defluorination_kinetics
        let rows: Vec<_> = crate::simulate_defluorination_kinetics(10.0, 50.0, 1000.0, 5.0, 50.0, 24.0).collect();
        assert_eq!(rows.len(), 1);
        let (t, f_ppm, fap_mg, frac) = rows[0];
        assert!((t - 24.0).abs() < 1e-6);
        assert!(f_ppm > 100.0);
        assert!(fap_mg > 2500.0);
        assert!(frac > 0.90);
    }

    #[pg_test]
    fn test_huckel_aromaticity_perception_and_wl_hash() {
        // 1. Benzene: Kekule (C1=CC=CC=C1) vs Aromatic (c1ccccc1)
        let benz_kek = crate::parse_smiles("C1=CC=CC=C1");
        let benz_aro = crate::parse_smiles("c1ccccc1");
        assert!(benz_kek.is_valid && benz_aro.is_valid);
        assert_eq!(benz_kek.atoms, benz_aro.atoms);
        assert_eq!(benz_kek.bonds, benz_aro.bonds);
        assert_eq!(benz_kek.wl_hash, benz_aro.wl_hash);
        assert_eq!(benz_kek.ecfp_bits, benz_aro.ecfp_bits);
        assert!(crate::molecules_are_identical("C1=CC=CC=C1", "c1ccccc1"));
        assert!((crate::molecule_tanimoto(benz_kek.clone(), benz_aro.clone()) - 1.0).abs() < 1e-9);

        // 2. Heterocycles & lone-pair donors (Pyridine, Pyrrole, Furan, Thiophene, Imidazole, 2-Pyridone, Cyclopentadienyl anion)
        let pairs = [
            ("C1=CC=NC=C1", "c1ccncc1"),           // Pyridine (6 pi e-)
            ("C1=CC=CN1", "c1cc[nH]c1"),           // Pyrrole (6 pi e-)
            ("O1C=CC=C1", "c1ccoc1"),              // Furan (6 pi e-)
            ("S1C=CC=C1", "c1ccsc1"),              // Thiophene (6 pi e-)
            ("C1=CN=CN1", "c1c[nH]cn1"),           // Imidazole (6 pi e-)
            ("O=C1C=CC=CN1", "O=c1cccc[nH]1"),     // 2-Pyridone (exocyclic C=O contributes 0 pi e-, total 6 pi e-)
            ("[CH-]1C=CC=C1", "[cH-]1cccc1"),      // Cyclopentadienyl anion (carbanion donates 2 pi e-, total 6 pi e-)
            ("C1=CC2=CC=CC=C2C=C1", "c1ccc2ccccc2c1"), // Naphthalene (10 pi e-)
            ("C1=CC=C2C=CC=CC2=C1", "c1ccc2ccccc2c1"), // Naphthalene alternative Kekule resonance form
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"), // Caffeine (fused 5+6 ring with two exocyclic =O)
            ("CC(=O)OC1=CC=CC=C1C(=O)O", "O=C(O)c1ccccc1OC(=O)C"), // Aspirin Kekule vs reversed aromatic
        ];
        for (kek, aro) in pairs {
            assert_eq!(
                crate::molecule_wl_hash(kek),
                crate::molecule_wl_hash(aro),
                "WL hash mismatch between {} and {}",
                kek,
                aro
            );
            assert!(
                crate::molecules_are_identical(kek, aro),
                "Expected molecules_are_identical({}, {}) == true",
                kek,
                aro
            );
        }

        // 3. Saturated sp3 and anti-aromatic (4n) controls must NOT be normalized to aromatic
        let cyclohexane = crate::parse_smiles("C1CCCCC1");
        assert!(cyclohexane.atoms.iter().all(|a| !a.is_aromatic));
        assert_ne!(cyclohexane.wl_hash, benz_aro.wl_hash);

        let cyclopentadiene = crate::parse_smiles("C1=CC=CC1"); // sp3 -CH2- breaks conjugation
        assert!(cyclopentadiene.atoms.iter().all(|a| !a.is_aromatic));

        let cyclobutadiene = crate::parse_smiles("C1=CC=C1"); // 4 pi electrons (anti-aromatic 4n)
        assert!(cyclobutadiene.atoms.iter().all(|a| !a.is_aromatic));

        let cyclooctatetraene = crate::parse_smiles("C1=CC=CC=CC=C1"); // 8 pi electrons (4n, non-planar)
        assert!(cyclooctatetraene.atoms.iter().all(|a| !a.is_aromatic));
    }

    #[pg_test]
    fn test_vf2_subgraph_isomorphism_and_operators() {
        let aspirin_kek = crate::parse_smiles("CC(=O)OC1=CC=CC=C1C(=O)O");
        let aspirin_aro = crate::parse_smiles("CC(=O)Oc1ccccc1C(=O)O");
        let benzene_kek = crate::parse_smiles("C1=CC=CC=C1");
        let benzene_aro = crate::parse_smiles("c1ccccc1");
        let pyridine = crate::parse_smiles("c1ccncc1");

        // Kekule Aspirin contains aromatic Benzene AND Kekule Benzene
        assert!(crate::smiles_contains(aspirin_kek.clone(), "c1ccccc1"));
        assert!(crate::smiles_contains(aspirin_aro.clone(), "C1=CC=CC=C1"));
        assert!(crate::smiles_contains(aspirin_aro.clone(), "C(=O)O"));
        assert!(crate::molecule_contains_molecule(aspirin_kek.clone(), benzene_aro.clone()));
        assert!(crate::molecule_contains_molecule(benzene_kek, benzene_aro));
        // Aspirin does NOT contain pyridine
        assert!(!crate::smiles_contains(aspirin_kek, "c1ccncc1"));
        assert!(!crate::molecule_contains_molecule(aspirin_aro, pyridine));

        // Verify SQL operators `@>` and `%` on Molecule via SPI
        let contains_sql = Spi::get_one::<bool>(
            "SELECT parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O') @> parse_smiles('c1ccccc1');"
        ).unwrap().unwrap();
        assert!(contains_sql);

        let tanimoto_mol_sql = Spi::get_one::<f64>(
            "SELECT parse_smiles('C1=CC=CC=C1') % parse_smiles('c1ccccc1');"
        ).unwrap().unwrap();
        assert!((tanimoto_mol_sql - 1.0).abs() < 1e-9);
    }

    #[pg_test]
    fn test_salt_stripping_and_tautomer_normalization() {
        // 1. Salt stripping + counter-ion charge neutralization
        // Sodium acetate `CC(=O)O.[Na+]` and `CC(=O)[O-].[Na+]` normalize to Acetic acid `CC(=O)O`
        assert!(crate::molecules_are_identical("CC(=O)O.[Na+]", "CC(=O)O"));
        assert!(crate::molecules_are_identical("CC(=O)[O-].[Na+]", "CC(=O)O"));
        assert!(crate::molecules_are_identical("[Na+].c1ccccc1C(=O)[O-]", "O=C(O)C1=CC=CC=C1"));
        assert!(crate::molecules_are_identical("[NH3+]CC(=O)O.[Cl-]", "NCC(=O)O"));

        // 2. 1,3-Prototropic tautomer normalization (aliphatic keto-enol & imidic-amide)
        // Vinyl alcohol `C=CO` normalizes to Acetaldehyde `CC=O`
        assert!(crate::molecules_are_identical("C=CO", "CC=O"));
        // Prop-1-en-2-ol `CC(=C)O` normalizes to Acetone `CC(=O)C`
        assert!(crate::molecules_are_identical("CC(=C)O", "CC(=O)C"));
        // Acetimidic acid `CC(=N)O` normalizes to Acetamide `CC(=O)N`
        assert!(crate::molecules_are_identical("CC(=N)O", "CC(=O)N"));
        // Phenol `Oc1ccccc1` preserves aromatic ring (does NOT tautomerize to cyclohexadienone)
        let phenol = crate::parse_smiles("Oc1ccccc1");
        assert_eq!(phenol.atoms.iter().filter(|a| a.is_aromatic).count(), 6);
    }
}

#[cfg(test)]
pub mod pg_test {
    pub fn setup(_options: Vec<&str>) { }
    pub fn postgresql_conf_options() -> Vec<&'static str> { vec![] }
}
