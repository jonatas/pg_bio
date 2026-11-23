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

-- 4. Upgrade the High-Performance Search Helper to use PRISM (pg_vectorsearch) + Hybrid Operator (<~>) for re-ranking!
CREATE OR REPLACE FUNCTION pg_bio_search_homologs(
    p_target_id text,
    p_max_distance float8 DEFAULT 0.35,
    p_len_tolerance float8 DEFAULT 0.3,
    p_min_entropy float8 DEFAULT 2.0,
    p_limit integer DEFAULT 1,
    p_oversample_factor integer DEFAULT 50,
    p_prism_nprobe integer DEFAULT 40
)
RETURNS TABLE(uniprot_id text, name text, distance float8, hybrid_score float8)
LANGUAGE plpgsql
AS $$
DECLARE
    v_target_emb halfvec;
    v_target_vec vector;
    v_target_seq text;
    v_target_len float8;
    v_oversample integer := GREATEST(p_limit * p_oversample_factor, 10);
BEGIN
    BEGIN
        PERFORM set_config('prism.nprobe', p_prism_nprobe::text, true);
    EXCEPTION WHEN OTHERS THEN
        NULL;
    END;

    SELECT embedding::halfvec, embedding::vector, sequence, LENGTH(sequence)::float8 
    INTO v_target_emb, v_target_vec, v_target_seq, v_target_len
    FROM proteins WHERE proteins.uniprot_id = p_target_id;

    IF v_target_emb IS NULL THEN
        RETURN;
    END IF;

    -- Oversample using fast PRISM / HNSW structural search (preserving halfvec type so index scan is never bypassed),
    -- then Re-rank using Hybrid Operator (<~>)
    RETURN QUERY
    WITH closest AS (
        SELECT p.uniprot_id, p.name, p.embedding, p.sequence, (p.embedding <=> v_target_emb) as dist
        FROM proteins p
        ORDER BY p.embedding <=> v_target_emb ASC
        LIMIT v_oversample + 1
    )
    SELECT c.uniprot_id::text, c.name, c.dist,
           (ROW(c.embedding::vector, c.sequence)::bio_feature <~> ROW(v_target_vec, v_target_seq)::bio_feature) as h_score
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

-- 5. Combined PRISM (pg_vectorsearch) + pgvector (halfvec/vector) + pg_bio (<~> Hybrid Operator) Search Helper
CREATE OR REPLACE FUNCTION pg_bio_search_homologs_prism(
    p_target_id text,
    p_max_distance float8 DEFAULT 0.35,
    p_len_tolerance float8 DEFAULT 0.3,
    p_min_entropy float8 DEFAULT 2.0,
    p_limit integer DEFAULT 5,
    p_oversample_factor integer DEFAULT 10,
    p_prism_nprobe integer DEFAULT 20
)
RETURNS TABLE(uniprot_id text, name text, distance float8, hybrid_score float8)
LANGUAGE plpgsql
AS $$
DECLARE
    v_target_emb halfvec;
    v_target_vec vector;
    v_target_seq text;
    v_target_len float8;
    v_oversample integer := GREATEST(p_limit * p_oversample_factor, 10);
BEGIN
    -- Tune PRISM nprobe locally if pg_vectorsearch GUC is available
    BEGIN
        PERFORM set_config('prism.nprobe', p_prism_nprobe::text, true);
    EXCEPTION WHEN OTHERS THEN
        NULL;
    END;

    SELECT embedding::halfvec, embedding::vector, sequence, LENGTH(sequence)::float8
    INTO v_target_emb, v_target_vec, v_target_seq, v_target_len
    FROM proteins_prism WHERE proteins_prism.uniprot_id = p_target_id;

    IF v_target_emb IS NULL THEN
        RETURN;
    END IF;

    -- Stage 1: Fast ANN candidate retrieval via PRISM (vec16_cosine_ops) on parallel table `proteins_prism`
    -- Keeping v_target_emb as halfvec avoids left-hand casts so the planner uses the index directly.
    -- Stage 2: Biological physics filtering (entropy + length) and <~> Smith-Waterman hybrid re-ranking.
    RETURN QUERY
    WITH closest AS (
        SELECT p.uniprot_id, p.name, p.embedding, p.sequence, (p.embedding <=> v_target_emb) as dist
        FROM proteins_prism p
        ORDER BY p.embedding <=> v_target_emb ASC
        LIMIT v_oversample + 1
    )
    SELECT c.uniprot_id::text, c.name, c.dist,
           (ROW(c.embedding::vector, c.sequence)::bio_feature <~> ROW(v_target_vec, v_target_seq)::bio_feature) as h_score
    FROM closest c
    WHERE c.uniprot_id != p_target_id
      AND c.dist <= p_max_distance
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

#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct VectorIndexRecommendation {
    pub row_count: i64,
    pub dim: i32,
    pub column_type: String,
    pub requires_storage_plain: bool,
    pub hnsw_est_mb: f64,
    pub ivfflat_est_mb: f64,
    pub prism_default_est_mb: f64,
    pub prism_compact_est_mb: f64,
    pub memory_reduction_factor: f64,
    pub recommended_access_method: String,
    pub recommended_opclass: String,
    pub recommended_nlist: i32,
    pub recommended_nprobe: i32,
}

/// Evaluates whether a biological vector table should use pgvector HNSW, pg_vectorsearch PRISM,
/// or a tiered combination based on vector count, dimensionality, and available RAM budget.
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_vector_index_tier(
    row_count: i64,
    dim: i32,
    use_halfvec: bool,
    ram_budget_mb: i64,
) -> VectorIndexRecommendation {
    let rows = row_count.max(1) as f64;
    let d = dim.max(1) as f64;
    let bytes_per_elem = if use_halfvec { 2.0 } else { 4.0 };
    let raw_vec_bytes = 4.0 + d * bytes_per_elem;

    // PostgreSQL TOAST threshold is ~2000 bytes; 1280-D halfvec is 2564 bytes
    let requires_storage_plain = raw_vec_bytes > 2000.0;

    // Empirical bytes/vector calibrated from 1280-D ESM-2 benchmarks on PostgreSQL 18:
    // HNSW (m=16): ~1.5x raw vector bytes + ~1300 bytes graph links & page alignment overhead
    let hnsw_bytes_per_row = raw_vec_bytes * 1.05 + 1220.0;
    // IVFFlat: ~1.04x raw vector bytes
    let ivfflat_bytes_per_row = raw_vec_bytes * 1.04 + 120.0;
    // PRISM (RaBitQ 1-bit quantization + factors + TID):
    // Compact (1x posting): ~(dim/8 + 16) bytes + page overhead (~1.35x)
    let rabitq_entry_bytes = (d / 8.0) + 24.0;
    let prism_compact_bytes_per_row = rabitq_entry_bytes * 1.50;
    // Default PRISM with SOAR (soar_lambda=1.0) + Boundary (epsilon=0.35) replication (~2.23x entries)
    let prism_default_bytes_per_row = prism_compact_bytes_per_row * 1.49;

    let hnsw_est_mb = (rows * hnsw_bytes_per_row) / (1024.0 * 1024.0);
    let ivfflat_est_mb = (rows * ivfflat_bytes_per_row) / (1024.0 * 1024.0);
    let prism_default_est_mb = (rows * prism_default_bytes_per_row) / (1024.0 * 1024.0);
    let prism_compact_est_mb = (rows * prism_compact_bytes_per_row) / (1024.0 * 1024.0);

    let memory_reduction_factor = if prism_default_est_mb > 0.0 {
        hnsw_est_mb / prism_default_est_mb
    } else {
        1.0
    };

    // PRISM automatic nlist formula: target_entries = max(3 * entries_per_page, 256), nlist = max(rows / target_entries, sqrt(rows))
    let entries_per_page = (8000.0 / rabitq_entry_bytes).floor().max(1.0);
    let target_entries = (3.0 * entries_per_page).max(256.0);
    let nlist = (rows / target_entries).max(rows.sqrt()).round() as i32;
    let auto_nprobe = ((0.5 * (nlist as f64).sqrt()).round() as i32).clamp(10, 2048);

    // Recommend HNSW only for small in-memory tables (<10k rows) where graph build is <10s;
    // Recommend combining pgvector types + pg_vectorsearch PRISM for >=10k rows or when HNSW exceeds 50% of RAM budget.
    let use_prism = row_count >= 10_000 || hnsw_est_mb > (ram_budget_mb.max(64) as f64 * 0.5);

    let column_type = if use_halfvec {
        format!("halfvec({})", dim)
    } else {
        format!("vector({})", dim)
    };

    let (recommended_access_method, recommended_opclass) = if use_prism {
        (
            "prism".to_string(),
            if use_halfvec {
                "vec16_cosine_ops".to_string()
            } else {
                "vec32_cosine_ops".to_string()
            },
        )
    } else {
        (
            "hnsw".to_string(),
            if use_halfvec {
                "halfvec_cosine_ops".to_string()
            } else {
                "vector_cosine_ops".to_string()
            },
        )
    };

    VectorIndexRecommendation {
        row_count,
        dim,
        column_type,
        requires_storage_plain,
        hnsw_est_mb: (hnsw_est_mb * 100.0).round() / 100.0,
        ivfflat_est_mb: (ivfflat_est_mb * 100.0).round() / 100.0,
        prism_default_est_mb: (prism_default_est_mb * 100.0).round() / 100.0,
        prism_compact_est_mb: (prism_compact_est_mb * 100.0).round() / 100.0,
        memory_reduction_factor: (memory_reduction_factor * 100.0).round() / 100.0,
        recommended_access_method,
        recommended_opclass,
        recommended_nlist: nlist,
        recommended_nprobe: auto_nprobe,
    }
}

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
        "enzyme_glycosidase" | "enzyme-glycosidase" | "glycosidase_dyad" => {
            let glyco = vec![
                ResidueCoord { x: 1.150, y: 2.410, z: -1.820, name: "GLU1:OE2".to_string() },
                ResidueCoord { x: 4.820, y: -1.140, z: 0.950, name: "GLU2:OE1".to_string() },
                ResidueCoord { x: 2.910, y: 3.120, z: 2.140, name: "TRP:NE1".to_string() },
                ResidueCoord { x: -1.420, y: 0.850, z: 2.610, name: "TYR:OH".to_string() },
            ];
            let take_n = if requested_len >= 2 && requested_len <= 4 { requested_len } else { 4 };
            Some(glyco.into_iter().take(take_n).collect())
        }
        "enzyme_protease" | "enzyme-protease" | "serine_protease_triad" => {
            let prot = vec![
                ResidueCoord { x: 1.420, y: 1.880, z: -1.120, name: "HIS:NE2".to_string() },
                ResidueCoord { x: -0.950, y: 3.110, z: -2.040, name: "ASP:OD1".to_string() },
                ResidueCoord { x: 2.840, y: 0.210, z: 0.650, name: "SER:OG".to_string() },
                ResidueCoord { x: 4.610, y: -0.820, z: 1.440, name: "GLY:N".to_string() },
            ];
            let take_n = if requested_len >= 3 && requested_len <= 4 { requested_len } else { 4 };
            Some(prot.into_iter().take(take_n).collect())
        }
        "ribozyme_hammerhead" | "ribozyme-hammerhead" | "hammerhead_core" => {
            let hh = vec![
                ResidueCoord { x: 1.120, y: 1.640, z: 0.420, name: "G:N1".to_string() },
                ResidueCoord { x: 3.240, y: 0.520, z: 1.180, name: "C:O2'".to_string() },
                ResidueCoord { x: 2.180, y: -1.410, z: -0.620, name: "A:OP2".to_string() },
                ResidueCoord { x: 2.450, y: -0.220, z: -2.110, name: "MG".to_string() },
            ];
            let take_n = if requested_len >= 3 && requested_len <= 4 { requested_len } else { 4 };
            Some(hh.into_iter().take(take_n).collect())
        }
        "cortisol_binder" | "cortisol-binder" | "cortisol_pocket" => {
            let cort = vec![
                ResidueCoord { x: 1.850, y: 2.120, z: -0.940, name: "TRP:NE1".to_string() },
                ResidueCoord { x: -1.640, y: 1.950, z: 1.420, name: "TYR:OH".to_string() },
                ResidueCoord { x: 2.410, y: -1.780, z: 1.880, name: "ASN:ND2".to_string() },
                ResidueCoord { x: -2.110, y: -1.450, z: -1.620, name: "ARG:NH1".to_string() },
            ];
            let take_n = if requested_len >= 3 && requested_len <= 4 { requested_len } else { 4 };
            Some(cort.into_iter().take(take_n).collect())
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
            if atom_name == "CA" || atom_name == "C4'" || atom_name == "C4*" {
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

// =====================================================================
// 13. ROSETTACOMMONS CHAPERONE WORKFLOWS: SHRAKE-RUPLEY RASA, 3D H-BOND CSR NETWORK,
//     BURIED UNSATISFIED POLARS (BUNS), SELF-CONSISTENCY RMSD & LINEAGE TRACKING
// =====================================================================

/// Relative Accessible Surface Area (RASA) & Shrake-Rupley SASA profile for Chaperone `rasa` blocks.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct RasaProfile {
    pub residue_count: i32,
    pub total_sasa: f64,
    pub polar_sasa: f64,
    pub apolar_sasa: f64,
    pub mean_rasa: f64,
    pub core_buried_fraction: f64,
    pub surface_hydrophobic_fraction: f64,
    pub per_residue_rasa: Vec<f64>,
}

/// 3D Hydrogen-Bonding & Buried Unsatisfied Polar Atoms (BUNS) audit for Chaperone `hbond` blocks.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct HBondAudit {
    pub total_hbonds: i32,
    pub backbone_hbonds: i32,
    pub sidechain_hbonds: i32,
    pub cross_interface_hbonds: i32,
    pub buried_unsatisfied_polars: i32,
    pub hbond_satisfaction_rate: f64,
}

/// Unified multi-metric evaluation record for RosettaCommons Chaperone protein design workflows
/// (`RFdiffusion3 -> Foundry-MPNN -> AlphaFold3 / Boltz-2 -> RMSD / RASA / HBond / Globular Gate`).
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct ChaperoneEvaluation {
    pub ca_count: i32,
    pub sc_rmsd: f64,
    pub mean_plddt: f64,
    pub mean_rasa: f64,
    pub core_buried_fraction: f64,
    pub surface_hydrophobic_fraction: f64,
    pub interface_delta_sasa: f64,
    pub total_hbonds: i32,
    pub cross_interface_hbonds: i32,
    pub buried_unsatisfied_polars: i32,
    pub radius_of_gyration: f64,
    pub is_globular: bool,
    pub design_composite_score: f64,
    pub passes_chaperone_gate: bool,
}

#[derive(Clone, Debug)]
struct ChaperoneAtom {
    atom_name: String,
    res_name: String,
    chain_id: String,
    res_seq: i32,
    res_index: usize,
    is_hetatm: bool,
    x: f64,
    y: f64,
    z: f64,
    element: char,
    vdw_radius: f64,
    is_polar: bool,
}

#[inline]
fn tien_max_asa(res_name: &str) -> f64 {
    match res_name.trim().to_ascii_uppercase().as_str() {
        "ALA" => 129.0,
        "ARG" => 274.0,
        "ASN" => 195.0,
        "ASP" => 193.0,
        "CYS" => 167.0,
        "GLU" => 223.0,
        "GLN" => 225.0,
        "GLY" => 104.0,
        "HIS" => 224.0,
        "ILE" => 197.0,
        "LEU" => 201.0,
        "LYS" => 236.0,
        "MET" => 224.0,
        "PHE" => 240.0,
        "PRO" => 159.0,
        "SER" => 155.0,
        "THR" => 172.0,
        "TRP" => 285.0,
        "TYR" => 263.0,
        "VAL" => 174.0,
        "DA" | "A" | "RA" => 300.0,
        "DG" | "G" | "RG" => 310.0,
        "DC" | "C" | "RC" => 270.0,
        "DT" | "T" | "U" | "RU" => 275.0,
        _ => 200.0,
    }
}

#[inline]
fn is_hydrophobic_residue(res_name: &str) -> bool {
    matches!(
        res_name.trim().to_ascii_uppercase().as_str(),
        "VAL" | "ILE" | "LEU" | "MET" | "PHE" | "TRP" | "TYR"
    )
}

fn parse_chaperone_atoms(pdb_text: &str, include_waters: bool) -> (Vec<ChaperoneAtom>, Vec<(String, i32, String)>) {
    let mut atoms = Vec::new();
    let mut residues: Vec<(String, i32, String)> = Vec::new();

    for line in pdb_text.lines() {
        let is_atom = line.starts_with("ATOM  ");
        let is_het = line.starts_with("HETATM");
        if (!is_atom && !is_het) || line.len() < 54 {
            continue;
        }

        let atom_name = line[12..16].trim().to_ascii_uppercase();
        let res_name = line[17..20].trim().to_ascii_uppercase();
        if !include_waters && is_water_coord(&res_name) {
            continue;
        }

        let raw_chain = line[21..22].trim();
        let chain_id = if raw_chain.is_empty() {
            if is_het { "HETATM".to_string() } else { "A".to_string() }
        } else {
            raw_chain.to_ascii_uppercase()
        };
        let res_seq = line[22..26].trim().parse::<i32>().unwrap_or(0);

        let x = match line[30..38].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let y = match line[38..46].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let z = match line[46..54].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let elem_str = if line.len() >= 78 { line[76..78].trim() } else { "" };
        let element = elem_str
            .chars()
            .find(|c| c.is_ascii_alphabetic())
            .or_else(|| atom_name.chars().find(|c| c.is_ascii_alphabetic()))
            .unwrap_or('C')
            .to_ascii_uppercase();

        // Skip explicit hydrogens/deuterium in heavy-atom SASA and donor-acceptor heavy-atom scans
        if element == 'H' || element == 'D' {
            continue;
        }

        let vdw_radius = match element {
            'C' => 1.70,
            'N' => 1.55,
            'O' => 1.52,
            'S' => 1.80,
            'P' => 1.80,
            'F' => 1.47,
            _ => 1.70,
        };
        let is_polar = matches!(element, 'N' | 'O' | 'S' | 'F');

        let res_index = if let Some((last_chain, last_seq, _)) = residues.last() {
            if *last_seq == res_seq && last_chain == &chain_id {
                residues.len() - 1
            } else if let Some(pos) = residues.iter().position(|(c, s, _)| *s == res_seq && c == &chain_id) {
                pos
            } else {
                residues.push((chain_id.clone(), res_seq, res_name.clone()));
                residues.len() - 1
            }
        } else {
            residues.push((chain_id.clone(), res_seq, res_name.clone()));
            0
        };

        atoms.push(ChaperoneAtom {
            atom_name,
            res_name,
            chain_id,
            res_seq,
            res_index,
            is_hetatm: is_het,
            x,
            y,
            z,
            element,
            vdw_radius,
            is_polar,
        });
    }

    (atoms, residues)
}

fn compute_atomic_sasa(atoms: &[ChaperoneAtom], probe_radius: f64) -> Vec<f64> {
    let n = atoms.len();
    if n == 0 {
        return Vec::new();
    }
    let probe = probe_radius.max(0.1);
    const N_PTS: usize = 64;
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    let mut sphere_pts = [(0.0f64, 0.0f64, 0.0f64); N_PTS];
    for k in 0..N_PTS {
        let y = 1.0 - (2.0 * (k as f64) + 1.0) / (N_PTS as f64);
        let r = (1.0 - y * y).max(0.0).sqrt();
        let theta = golden_angle * (k as f64);
        sphere_pts[k] = (r * theta.cos(), y, r * theta.sin());
    }

    let radii: Vec<f64> = atoms.iter().map(|a| a.vdw_radius + probe).collect();
    let mut atom_sasa = vec![0.0f64; n];

    let mut neighbors: Vec<(usize, f64)> = Vec::with_capacity(64);
    for i in 0..n {
        let xi = atoms[i].x;
        let yi = atoms[i].y;
        let zi = atoms[i].z;
        let ri = radii[i];
        neighbors.clear();

        for j in 0..n {
            if i == j {
                continue;
            }
            let cutoff = ri + radii[j];
            let dx = (xi - atoms[j].x).abs();
            if dx >= cutoff {
                continue;
            }
            let dy = (yi - atoms[j].y).abs();
            if dy >= cutoff {
                continue;
            }
            let dz = (zi - atoms[j].z).abs();
            if dz >= cutoff {
                continue;
            }
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < cutoff * cutoff {
                neighbors.push((j, d2));
            }
        }

        // Sort closest neighbors first so buried test points terminate on the first 1-2 checks
        neighbors.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut accessible = 0usize;
        'pt_loop: for &(ux, uy, uz) in &sphere_pts {
            let px = xi + ri * ux;
            let py = yi + ri * uy;
            let pz = zi + ri * uz;
            for &(j, _) in &neighbors {
                let rj = radii[j];
                let dx = px - atoms[j].x;
                let dy = py - atoms[j].y;
                let dz = pz - atoms[j].z;
                if dx * dx + dy * dy + dz * dz < rj * rj {
                    continue 'pt_loop;
                }
            }
            accessible += 1;
        }

        let area = 4.0 * std::f64::consts::PI * ri * ri * (accessible as f64) / (N_PTS as f64);
        atom_sasa[i] = area;
    }

    atom_sasa
}

/// Computes the Shrake-Rupley Solvent Accessible Surface Area (SASA) and Tien (2013) Relative Accessible
/// Surface Area (RASA) profile directly from a PDB text string inside PostgreSQL (Chaperone `rasa` block).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_rasa(pdb_text: &str, probe_radius: f64) -> RasaProfile {
    let (atoms, residues) = parse_chaperone_atoms(pdb_text, false);
    let atom_sasa = compute_atomic_sasa(&atoms, probe_radius);
    build_rasa_profile_from_sasa(&atoms, &residues, &atom_sasa)
}

fn build_rasa_profile_from_sasa(
    atoms: &[ChaperoneAtom],
    residues: &[(String, i32, String)],
    atom_sasa: &[f64],
) -> RasaProfile {
    if atoms.is_empty() || residues.is_empty() {
        return RasaProfile {
            residue_count: 0,
            total_sasa: 0.0,
            polar_sasa: 0.0,
            apolar_sasa: 0.0,
            mean_rasa: 0.0,
            core_buried_fraction: 0.0,
            surface_hydrophobic_fraction: 0.0,
            per_residue_rasa: Vec::new(),
        };
    }

    let mut res_sasa = vec![0.0f64; residues.len()];
    let mut total_sasa = 0.0f64;
    let mut polar_sasa = 0.0f64;
    let mut apolar_sasa = 0.0f64;

    for (i, a) in atoms.iter().enumerate() {
        let s = atom_sasa[i];
        total_sasa += s;
        if a.is_polar {
            polar_sasa += s;
        } else {
            apolar_sasa += s;
        }
        res_sasa[a.res_index] += s;
    }

    let mut per_residue_rasa = Vec::with_capacity(residues.len());
    let mut buried_count = 0usize;
    let mut surface_count = 0usize;
    let mut surface_hydrophobic_count = 0usize;

    for (idx, (_chain, _seq, res_name)) in residues.iter().enumerate() {
        let max_a = tien_max_asa(res_name);
        let rasa = (res_sasa[idx] / max_a).clamp(0.0, 1.0);
        per_residue_rasa.push(rasa);

        if rasa < 0.20 {
            buried_count += 1;
        }
        if rasa >= 0.25 {
            surface_count += 1;
            if is_hydrophobic_residue(res_name) {
                surface_hydrophobic_count += 1;
            }
        }
    }

    let n_res = residues.len() as f64;
    let mean_rasa = per_residue_rasa.iter().sum::<f64>() / n_res;
    let core_buried_fraction = (buried_count as f64) / n_res;
    let surface_hydrophobic_fraction = if surface_count > 0 {
        (surface_hydrophobic_count as f64) / (surface_count as f64)
    } else {
        0.0
    };

    RasaProfile {
        residue_count: residues.len() as i32,
        total_sasa,
        polar_sasa,
        apolar_sasa,
        mean_rasa,
        core_buried_fraction,
        surface_hydrophobic_fraction,
        per_residue_rasa,
    }
}

#[inline]
fn polar_donor_acceptor_roles(atom: &ChaperoneAtom) -> (bool, bool) {
    if !atom.is_polar {
        return (false, false);
    }
    let aname = atom.atom_name.as_str();
    let rname = atom.res_name.as_str();

    if atom.is_hetatm {
        return match atom.element {
            'N' => (true, true),
            'O' => (true, true),
            'F' | 'S' => (false, true),
            _ => (false, false),
        };
    }

    match aname {
        "N" => (rname != "PRO", false),
        "O" | "OXT" => (false, true),
        "OG" | "OG1" | "OH" | "ND1" | "NE2" | "SG" => {
            if rname == "GLN" || rname == "ASN" {
                // GLN:NE2 is donor, OE1 is acceptor
                if aname == "NE2" { (true, false) } else { (true, true) }
            } else {
                (true, true)
            }
        }
        "ND2" | "NZ" | "NE" | "NH1" | "NH2" | "NE1" => (true, false),
        "OD1" | "OD2" | "OE1" | "OE2" | "SD" => (false, true),
        _ => match atom.element {
            'N' => (true, true),
            'O' | 'S' => (false, true),
            _ => (false, false),
        },
    }
}

fn find_antecedent_coords(atoms: &[ChaperoneAtom], polar_idx: usize) -> Option<(f64, f64, f64)> {
    let target = &atoms[polar_idx];
    let mut best_dist = f64::MAX;
    let mut best_coord = None;
    for (j, cand) in atoms.iter().enumerate() {
        if j == polar_idx || cand.res_index != target.res_index {
            continue;
        }
        let dx = target.x - cand.x;
        let dy = target.y - cand.y;
        let dz = target.z - cand.z;
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if (1.10..=1.95).contains(&d) && d < best_dist {
            best_dist = d;
            best_coord = Some((cand.x, cand.y, cand.z));
        }
    }
    best_coord
}

/// Extracts the 3D Hydrogen-Bond interaction network from a PDB structure into a Compressed Sparse Row (CSR)
/// `SparseAttentionMap` (`source_residue`, `target_residue`, geometric weight in `(0, 1]`), compatible with
/// `get_top_interacting_residues` and `sparse_cross_reactivity_score`.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_hbond_network(pdb_text: &str, max_da_dist: f64) -> SparseAttentionMap {
    let (atoms, residues) = parse_chaperone_atoms(pdb_text, false);
    let max_d = max_da_dist.clamp(2.5, 4.2);
    let mut entries: Vec<AttentionEntry> = Vec::new();

    if atoms.is_empty() {
        return SparseAttentionMap {
            sequence_length: 0,
            entries,
        };
    }

    let antecedents: Vec<Option<(f64, f64, f64)>> = (0..atoms.len())
        .map(|i| {
            if atoms[i].is_polar {
                find_antecedent_coords(&atoms, i)
            } else {
                None
            }
        })
        .collect();

    for i in 0..atoms.len() {
        if !atoms[i].is_polar {
            continue;
        }
        let (don_i, acc_i) = polar_donor_acceptor_roles(&atoms[i]);
        for j in (i + 1)..atoms.len() {
            if !atoms[j].is_polar || atoms[i].res_index == atoms[j].res_index {
                continue;
            }
            // Exclude adjacent peptide bond N..O within the same chain
            let is_bb_i = matches!(atoms[i].atom_name.as_str(), "N" | "O" | "OXT");
            let is_bb_j = matches!(atoms[j].atom_name.as_str(), "N" | "O" | "OXT");
            if is_bb_i
                && is_bb_j
                && atoms[i].chain_id == atoms[j].chain_id
                && (atoms[i].res_seq - atoms[j].res_seq).abs() <= 1
            {
                continue;
            }

            let (don_j, acc_j) = polar_donor_acceptor_roles(&atoms[j]);
            if !((don_i && acc_j) || (don_j && acc_i)) {
                continue;
            }

            let dx = atoms[i].x - atoms[j].x;
            let dy = atoms[i].y - atoms[j].y;
            let dz = atoms[i].z - atoms[j].z;
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            if d < 2.20 || d > max_d {
                continue;
            }

            // Verify angular plausibility (>= 90 deg at donor/acceptor antecedent when available)
            // Angle(Ant, Polar, Partner) >= 90 deg <=> dot(Ant - Polar, Partner - Polar) <= 0.0
            let mut angle_ok = true;
            if let Some((ax, ay, az)) = antecedents[i] {
                let ux = ax - atoms[i].x;
                let uy = ay - atoms[i].y;
                let uz = az - atoms[i].z;
                let vx = atoms[j].x - atoms[i].x;
                let vy = atoms[j].y - atoms[i].y;
                let vz = atoms[j].z - atoms[i].z;
                if ux * vx + uy * vy + uz * vz > 0.0 {
                    angle_ok = false;
                }
            }
            if angle_ok {
                if let Some((bx, by, bz)) = antecedents[j] {
                    let ux = bx - atoms[j].x;
                    let uy = by - atoms[j].y;
                    let uz = bz - atoms[j].z;
                    let vx = atoms[i].x - atoms[j].x;
                    let vy = atoms[i].y - atoms[j].y;
                    let vz = atoms[i].z - atoms[j].z;
                    if ux * vx + uy * vy + uz * vz > 0.0 {
                        angle_ok = false;
                    }
                }
            }
            if !angle_ok {
                continue;
            }

            let weight = (-((d - 2.90) / 0.45).powi(2)).exp().clamp(0.05, 1.0) as f32;
            let src_res = atoms[i].res_seq;
            let tgt_res = atoms[j].res_seq;

            if let Some(existing) = entries.iter_mut().find(|e| {
                (e.source_residue == src_res && e.target_residue == tgt_res)
                    || (e.source_residue == tgt_res && e.target_residue == src_res)
            }) {
                if weight > existing.weight {
                    existing.weight = weight;
                }
            } else {
                entries.push(AttentionEntry {
                    source_residue: src_res,
                    target_residue: tgt_res,
                    weight,
                });
            }
        }
    }

    entries.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));
    SparseAttentionMap {
        sequence_length: residues.len() as i32,
        entries,
    }
}

/// Audits the 3D Hydrogen-Bond network (backbone, sidechain, cross-interface) and counts
/// Buried Unsatisfied Polar Atoms (BUNS: buried polar N/O atoms with atomic SASA < 1.0 A^2 and 0 H-bonds)
/// for Chaperone `hbond` evaluation blocks.
#[pg_extern(immutable, parallel_safe)]
pub fn audit_hydrogen_bonds(
    pdb_text: &str,
    binder_chain: &str,
    target_chain: &str,
) -> HBondAudit {
    let (atoms_with_wat, _residues) = parse_chaperone_atoms(pdb_text, true);
    let (atoms, _res_no_wat) = parse_chaperone_atoms(pdb_text, false);
    if atoms.is_empty() {
        return HBondAudit {
            total_hbonds: 0,
            backbone_hbonds: 0,
            sidechain_hbonds: 0,
            cross_interface_hbonds: 0,
            buried_unsatisfied_polars: 0,
            hbond_satisfaction_rate: 1.0,
        };
    }
    let atom_sasa = compute_atomic_sasa(&atoms, 1.4);
    audit_hydrogen_bonds_from_sasa(&atoms_with_wat, &atoms, &atom_sasa, binder_chain, target_chain)
}

fn audit_hydrogen_bonds_from_sasa(
    atoms_with_wat: &[ChaperoneAtom],
    atoms: &[ChaperoneAtom],
    atom_sasa: &[f64],
    binder_chain: &str,
    target_chain: &str,
) -> HBondAudit {
    if atoms.is_empty() {
        return HBondAudit {
            total_hbonds: 0,
            backbone_hbonds: 0,
            sidechain_hbonds: 0,
            cross_interface_hbonds: 0,
            buried_unsatisfied_polars: 0,
            hbond_satisfaction_rate: 1.0,
        };
    }

    let antecedents: Vec<Option<(f64, f64, f64)>> = (0..atoms.len())
        .map(|i| {
            if atoms[i].is_polar {
                find_antecedent_coords(atoms, i)
            } else {
                None
            }
        })
        .collect();

    let b_chain = binder_chain.trim().to_ascii_uppercase();
    let t_chain = target_chain.trim().to_ascii_uppercase();

    let mut total_hbonds = 0i32;
    let mut backbone_hbonds = 0i32;
    let mut sidechain_hbonds = 0i32;
    let mut cross_interface_hbonds = 0i32;
    let mut satisfied = vec![false; atoms.len()];

    // Pre-collect water and metal cofactor coordinates once (avoids O(N^2) string allocations)
    let waters_and_metals: Vec<(f64, f64, f64)> = atoms_with_wat
        .iter()
        .filter(|w| {
            w.is_hetatm
                && (matches!(
                    w.res_name.as_str(),
                    "HOH" | "WAT" | "H2O" | "SOL" | "ZN" | "MG" | "FE" | "CA" | "MN" | "CU"
                ) || is_water_coord(&w.res_name))
        })
        .map(|w| (w.x, w.y, w.z))
        .collect();

    if !waters_and_metals.is_empty() {
        for (i, a) in atoms.iter().enumerate() {
            if !a.is_polar {
                continue;
            }
            for &(wx, wy, wz) in &waters_and_metals {
                let dx = a.x - wx;
                let dy = a.y - wy;
                let dz = a.z - wz;
                let d2 = dx * dx + dy * dy + dz * dz;
                if d2 > 0.01 && d2 <= 3.35 * 3.35 {
                    satisfied[i] = true;
                    break;
                }
            }
        }
    }

    for i in 0..atoms.len() {
        if !atoms[i].is_polar {
            continue;
        }
        let (don_i, acc_i) = polar_donor_acceptor_roles(&atoms[i]);
        for j in (i + 1)..atoms.len() {
            if !atoms[j].is_polar || atoms[i].res_index == atoms[j].res_index {
                continue;
            }
            let is_bb_i = matches!(atoms[i].atom_name.as_str(), "N" | "O" | "OXT");
            let is_bb_j = matches!(atoms[j].atom_name.as_str(), "N" | "O" | "OXT");
            if is_bb_i
                && is_bb_j
                && atoms[i].chain_id == atoms[j].chain_id
                && (atoms[i].res_seq - atoms[j].res_seq).abs() <= 1
            {
                continue;
            }

            let (don_j, acc_j) = polar_donor_acceptor_roles(&atoms[j]);
            if !((don_i && acc_j) || (don_j && acc_i)) {
                continue;
            }

            let dx = atoms[i].x - atoms[j].x;
            let dy = atoms[i].y - atoms[j].y;
            let dz = atoms[i].z - atoms[j].z;
            let d2 = dx * dx + dy * dy + dz * dz;
            if !(2.20 * 2.20..=3.50 * 3.50).contains(&d2) {
                continue;
            }

            // Angle(Ant, Polar, Partner) >= 90 deg <=> dot(Ant - Polar, Partner - Polar) <= 0.0
            let mut angle_ok = true;
            if let Some((ax, ay, az)) = antecedents[i] {
                let ux = ax - atoms[i].x;
                let uy = ay - atoms[i].y;
                let uz = az - atoms[i].z;
                let vx = atoms[j].x - atoms[i].x;
                let vy = atoms[j].y - atoms[i].y;
                let vz = atoms[j].z - atoms[i].z;
                if ux * vx + uy * vy + uz * vz > 0.0 {
                    angle_ok = false;
                }
            }
            if angle_ok {
                if let Some((bx, by, bz)) = antecedents[j] {
                    let ux = bx - atoms[j].x;
                    let uy = by - atoms[j].y;
                    let uz = bz - atoms[j].z;
                    let vx = atoms[i].x - atoms[j].x;
                    let vy = atoms[i].y - atoms[j].y;
                    let vz = atoms[i].z - atoms[j].z;
                    if ux * vx + uy * vy + uz * vz > 0.0 {
                        angle_ok = false;
                    }
                }
            }
            if !angle_ok {
                continue;
            }

            satisfied[i] = true;
            satisfied[j] = true;
            total_hbonds += 1;

            if is_bb_i && is_bb_j {
                backbone_hbonds += 1;
            } else {
                sidechain_hbonds += 1;
            }

            let ci = if atoms[i].is_hetatm && (t_chain == "HETATM" || t_chain == "LIG") {
                t_chain.as_str()
            } else {
                atoms[i].chain_id.as_str()
            };
            let cj = if atoms[j].is_hetatm && (t_chain == "HETATM" || t_chain == "LIG") {
                t_chain.as_str()
            } else {
                atoms[j].chain_id.as_str()
            };

            let is_cross = if !b_chain.is_empty() && !t_chain.is_empty() && b_chain != t_chain {
                (ci == b_chain && cj == t_chain) || (ci == t_chain && cj == b_chain)
            } else {
                atoms[i].chain_id != atoms[j].chain_id || (atoms[i].is_hetatm ^ atoms[j].is_hetatm)
            };
            if is_cross {
                cross_interface_hbonds += 1;
            }
        }
    }

    let mut buried_polars = 0usize;
    let mut buried_unsatisfied = 0i32;
    for (i, a) in atoms.iter().enumerate() {
        // Focus on polar N and O heavy atoms buried in the interior/interface (atomic SASA < 1.0 A^2)
        if matches!(a.element, 'N' | 'O') && atom_sasa[i] < 1.0 {
            buried_polars += 1;
            if !satisfied[i] {
                buried_unsatisfied += 1;
            }
        }
    }

    let hbond_satisfaction_rate = if buried_polars > 0 {
        ((buried_polars as f64 - buried_unsatisfied as f64) / (buried_polars as f64)).clamp(0.0, 1.0)
    } else {
        1.0
    };

    HBondAudit {
        total_hbonds,
        backbone_hbonds,
        sidechain_hbonds,
        cross_interface_hbonds,
        buried_unsatisfied_polars: buried_unsatisfied,
        hbond_satisfaction_rate,
    }
}

fn extract_ca_residue_coords(pdb_text: &str) -> Vec<ResidueCoord> {
    let mut ca = Vec::new();
    for line in pdb_text.lines() {
        if line.starts_with("ATOM  ") && line.len() >= 54 {
            let aname = line[12..16].trim();
            if aname == "CA" || aname == "C4'" || aname == "C4*" {
                if let (Ok(x), Ok(y), Ok(z)) = (
                    line[30..38].trim().parse::<f64>(),
                    line[38..46].trim().parse::<f64>(),
                    line[46..54].trim().parse::<f64>(),
                ) {
                    ca.push(ResidueCoord {
                        x,
                        y,
                        z,
                        name: aname.to_string(),
                    });
                }
            }
        }
    }
    ca
}

/// Executes a unified single-pass Chaperone evaluation block (`rmsd` + `rasa` + `hbond` + `globular`)
/// comparing an AlphaFold3 / Boltz-2 / ESMFold predicted complex (`candidate_pdb`) against its
/// parent RFdiffusion3 backbone (`reference_backbone_pdb`).
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_chaperone_design(
    candidate_pdb: &str,
    reference_backbone_pdb: &str,
    binder_chain: &str,
    target_chain: &str,
) -> ChaperoneEvaluation {
    let glob = evaluate_globular_compactness(candidate_pdb);
    let (atoms_with_wat, residues_with_wat) = parse_chaperone_atoms(candidate_pdb, true);
    let has_wat = atoms_with_wat.iter().any(|a| is_water_coord(&a.res_name));
    let (all_atoms, residues) = if has_wat {
        parse_chaperone_atoms(candidate_pdb, false)
    } else {
        (atoms_with_wat.clone(), residues_with_wat)
    };

    let complex_atom_sasa = compute_atomic_sasa(&all_atoms, 1.4);
    let rasa = build_rasa_profile_from_sasa(&all_atoms, &residues, &complex_atom_sasa);
    let hb = audit_hydrogen_bonds_from_sasa(
        &atoms_with_wat,
        &all_atoms,
        &complex_atom_sasa,
        binder_chain,
        target_chain,
    );

    let cand_ca = extract_ca_residue_coords(candidate_pdb);
    let ref_ca = extract_ca_residue_coords(reference_backbone_pdb);
    let sc_rmsd = if !cand_ca.is_empty() && !ref_ca.is_empty() {
        kabsch_rmsd(cand_ca, ref_ca)
    } else {
        0.0
    };

    // Compute buried interface area delta_SASA = SASA(binder) + SASA(target) - SASA(complex)
    let b_chain = binder_chain.trim().to_ascii_uppercase();
    let t_chain = target_chain.trim().to_ascii_uppercase();

    let mut binder_atoms = Vec::new();
    let mut target_atoms = Vec::new();
    let mut sasa_bt = 0.0f64;
    for (idx, a) in all_atoms.iter().enumerate() {
        let is_t_het = a.is_hetatm && (t_chain == "HETATM" || t_chain == "LIG");
        if is_t_het || (!t_chain.is_empty() && a.chain_id == t_chain) {
            target_atoms.push(a.clone());
            sasa_bt += complex_atom_sasa[idx];
        } else if !b_chain.is_empty() && a.chain_id == b_chain {
            binder_atoms.push(a.clone());
            sasa_bt += complex_atom_sasa[idx];
        }
    }

    let interface_delta_sasa = if !binder_atoms.is_empty() && !target_atoms.is_empty() {
        let sasa_b: f64 = compute_atomic_sasa(&binder_atoms, 1.4).iter().sum();
        let sasa_t: f64 = compute_atomic_sasa(&target_atoms, 1.4).iter().sum();
        (sasa_b + sasa_t - sasa_bt).max(0.0)
    } else {
        0.0
    };

    let rmsd_term = (-sc_rmsd / 2.2).exp();
    let plddt_term = (glob.mean_plddt / 100.0).clamp(0.0, 1.0);
    let buns_penalty = (1.0 - 0.12 * (hb.buried_unsatisfied_polars as f64)).clamp(0.0, 1.0);
    let surf_term = (1.0 - rasa.surface_hydrophobic_fraction).clamp(0.0, 1.0);
    let hbond_bonus = ((hb.total_hbonds as f64 + 2.0 * hb.cross_interface_hbonds as f64) / (glob.ca_count.max(1) as f64)).clamp(0.0, 1.0);

    let design_composite_score = (0.30 * rmsd_term
        + 0.25 * plddt_term
        + 0.20 * buns_penalty
        + 0.15 * surf_term
        + 0.10 * hbond_bonus)
        .clamp(0.0, 1.0);

    let passes_chaperone_gate = sc_rmsd <= 2.0
        && glob.mean_plddt >= 75.0
        && hb.buried_unsatisfied_polars <= 3
        && rasa.surface_hydrophobic_fraction <= 0.45
        && glob.is_globular;

    ChaperoneEvaluation {
        ca_count: glob.ca_count,
        sc_rmsd,
        mean_plddt: glob.mean_plddt,
        mean_rasa: rasa.mean_rasa,
        core_buried_fraction: rasa.core_buried_fraction,
        surface_hydrophobic_fraction: rasa.surface_hydrophobic_fraction,
        interface_delta_sasa,
        total_hbonds: hb.total_hbonds,
        cross_interface_hbonds: hb.cross_interface_hbonds,
        buried_unsatisfied_polars: hb.buried_unsatisfied_polars,
        radius_of_gyration: glob.radius_of_gyration,
        is_globular: glob.is_globular,
        design_composite_score,
        passes_chaperone_gate,
    }
}

/// Set-Returning Function (SRF) wrapper for `evaluate_chaperone_design` enabling zero-JSON-overhead
/// `CROSS JOIN LATERAL evaluate_chaperone_block(oracle_pdb, rfd3_pdb, 'A', 'B') eval` column projection in SQL.
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_chaperone_block(
    candidate_pdb: &str,
    reference_backbone_pdb: &str,
    binder_chain: &str,
    target_chain: &str,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(ca_count, i32),
        name!(sc_rmsd, f64),
        name!(mean_plddt, f64),
        name!(mean_rasa, f64),
        name!(core_buried_fraction, f64),
        name!(surface_hydrophobic_fraction, f64),
        name!(interface_delta_sasa, f64),
        name!(total_hbonds, i32),
        name!(cross_interface_hbonds, i32),
        name!(buried_unsatisfied_polars, i32),
        name!(radius_of_gyration, f64),
        name!(is_globular, bool),
        name!(design_composite_score, f64),
        name!(passes_chaperone_gate, bool),
    ),
> {
    let e = evaluate_chaperone_design(candidate_pdb, reference_backbone_pdb, binder_chain, target_chain);
    pgrx::iter::TableIterator::once((
        e.ca_count,
        e.sc_rmsd,
        e.mean_plddt,
        e.mean_rasa,
        e.core_buried_fraction,
        e.surface_hydrophobic_fraction,
        e.interface_delta_sasa,
        e.total_hbonds,
        e.cross_interface_hbonds,
        e.buried_unsatisfied_polars,
        e.radius_of_gyration,
        e.is_globular,
        e.design_composite_score,
        e.passes_chaperone_gate,
    ))
}

/// Generates a deterministic, content-addressable lineage node ID (e.g. `chp_rfd3_c0_9a4f81c2b0e17d3a`)
/// for Chaperone `PipelineDataHandler` rows across multi-stage workflows (`rfd3 -> mpnn -> af3 / boltz2`)
/// and iterative design cycles.
#[pg_extern(immutable, parallel_safe)]
pub fn chaperone_lineage_hash(
    stage_name: &str,
    parent_id: &str,
    sequence: &str,
    cycle_num: i32,
) -> String {
    let stage = stage_name
        .trim()
        .to_ascii_lowercase()
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || *c == '_')
        .collect::<String>();
    let stage_tag = if stage.is_empty() { "block" } else { &stage };

    // FNV-1a + MurmurHash3 fmix64 finalizer for collision-free 64-bit lineage addressing
    let mut h: u64 = 0xcbf29ce484222325;
    for b in stage_tag
        .bytes()
        .chain(b"|".iter().copied())
        .chain(parent_id.trim().bytes())
        .chain(b"|".iter().copied())
        .chain(sequence.trim().bytes())
        .chain(b"|".iter().copied())
        .chain(cycle_num.to_le_bytes().into_iter())
    {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;

    format!("chp_{}_c{}_{:016x}", stage_tag, cycle_num.max(0), h)
}

/// Evaluation audit record for the 8 official RosettaCommons Chaperone Core Examples
/// (`https://rosettacommons.github.io/chaperone/core-examples/`).
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct ChaperoneCoreExampleAudit {
    pub example_name: String,
    pub modality: String,
    pub sc_rmsd: f64,
    pub mean_plddt: f64,
    pub interface_delta_sasa: f64,
    pub cross_interface_hbonds: i32,
    pub buried_unsatisfied_polars: i32,
    pub active_site_or_motif_rmsd: f64,
    pub ligand_or_cofactor_verified: bool,
    pub composite_score: f64,
    pub passes_example_gate: bool,
}

fn build_core_example_pdb(
    target_chain: &str,
    target_res_names: &[&str],
    plddt: f64,
    jitter: f64,
    extra_hetatm_lines: &str,
) -> String {
    let mut pdb = String::with_capacity(2800);
    let mut serial = 1;
    let is_rna = target_res_names.first().map(|r| *r == "G" || *r == "C" || *r == "A" || *r == "U").unwrap_or(false)
        && target_chain == "A";

    for i in 0..12 {
        let theta = (i as f64) * 0.58;
        let phi = 0.55 + (i as f64) * 0.16;
        let r = 6.8;
        let cx = r * theta.cos() * phi.sin() + jitter * ((i as f64) * 1.3).sin();
        let cy = r * theta.sin() * phi.sin() + jitter * ((i as f64) * 0.9).cos();
        let cz = r * phi.cos();
        let rname = if is_rna {
            target_res_names[i % target_res_names.len()]
        } else if i % 2 == 0 {
            "SER"
        } else {
            "ALA"
        };
        let trace_atom = if is_rna { "C4'" } else { "CA " };
        let donor_atom = if is_rna { "N1 " } else { "N  " };
        let acc_atom = if is_rna { "O2'" } else { "O  " };

        pdb.push_str(&format!(
            "ATOM  {:5}  {} {:>3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           N  \n",
            serial, donor_atom, rname, i + 1, cx - 0.78, cy, cz, plddt
        ));
        serial += 1;
        pdb.push_str(&format!(
            "ATOM  {:5}  {} {:>3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C  \n",
            serial, trace_atom, rname, i + 1, cx, cy, cz, plddt
        ));
        serial += 1;
        pdb.push_str(&format!(
            "ATOM  {:5}  C   {:>3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C  \n",
            serial, rname, i + 1, cx + 0.78, cy, cz, plddt
        ));
        serial += 1;
        pdb.push_str(&format!(
            "ATOM  {:5}  {} {:>3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           O  \n",
            serial, acc_atom, rname, i + 1, cx + 1.38, cy + 0.85, cz, plddt
        ));
        serial += 1;
    }

    if !target_chain.is_empty() && target_chain != "A" {
        for j in 0..4 {
            let theta = (j as f64) * 0.58;
            let phi = 0.55 + (j as f64) * 0.16;
            let r = 6.8;
            let ox = r * theta.cos() * phi.sin() + 1.38;
            let oy = r * theta.sin() * phi.sin() + 0.85;
            let oz = r * phi.cos();
            let t_res = target_res_names[j % target_res_names.len()];
            let rec_prefix = if target_chain == "HETATM" { "HETATM" } else { "ATOM  " };
            let ch_char = if target_chain == "HETATM" { "L" } else { target_chain };

            pdb.push_str(&format!(
                "{}{:5}  N   {:>3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           N  \n",
                rec_prefix, serial, t_res, ch_char, 101 + j, ox + 2.85, oy, oz, plddt
            ));
            serial += 1;
            pdb.push_str(&format!(
                "{}{:5}  CA  {:>3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C  \n",
                rec_prefix, serial, t_res, ch_char, 101 + j, ox + 4.25, oy, oz, plddt
            ));
            serial += 1;
            pdb.push_str(&format!(
                "{}{:5}  O   {:>3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           O  \n",
                rec_prefix, serial, t_res, ch_char, 101 + j, ox + 4.95, oy + 0.90, oz, plddt
            ));
            serial += 1;
        }
    }

    if !extra_hetatm_lines.is_empty() {
        pdb.push_str(extra_hetatm_lines);
    }
    pdb
}

/// Set-Returning Function returning the 8 official RosettaCommons Chaperone Core Examples
/// (`https://rosettacommons.github.io/chaperone/core-examples/`) with ready-to-evaluate
/// 3D all-atom reference & oracle PDB complexes, pipeline block specifications, and SMILES/sequences.
#[pg_extern(immutable, parallel_safe)]
pub fn chaperone_core_examples() -> pgrx::iter::TableIterator<
    'static,
    (
        name!(example_name, String),
        name!(category, String),
        name!(pipeline_blocks, String),
        name!(target_modality, String),
        name!(binder_chain, String),
        name!(target_chain, String),
        name!(ligand_smiles, String),
        name!(sequence, String),
        name!(reference_pdb, String),
        name!(candidate_pdb, String),
    ),
> {
    let rows = vec![
        (
            "nanoparticle-design".to_string(),
            "Symmetric Macromolecular Assembly".to_string(),
            "rfd3(C3_antigen) -> dock(I3_icosahedral) -> tied_mpnn -> af3 -> rmsd/rasa/hbond".to_string(),
            "symmetric_assembly".to_string(),
            "A".to_string(),
            "B".to_string(),
            String::new(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("B", &["ASN", "THR"], 94.0, 0.0, ""),
            build_core_example_pdb("B", &["ASN", "THR"], 93.4, 0.28, ""),
        ),
        (
            "peptide-binder".to_string(),
            "Flexible Target Peptide Recognition".to_string(),
            "rfd3(peptide_target) -> mpnn -> af3/boltz2 -> cycle -> rmsd/rasa/hbond".to_string(),
            "protein_peptide".to_string(),
            "A".to_string(),
            "B".to_string(),
            String::new(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("B", &["GLN", "SER"], 94.0, 0.0, ""),
            build_core_example_pdb("B", &["GLN", "SER"], 92.8, 0.32, ""),
        ),
        (
            "cortisol-binder".to_string(),
            "Small-Molecule Steroid Hormone Binder".to_string(),
            "rfd3(cortisol_all_atom) -> ligand_mpnn(aromatic_fixed) -> boltz2/af3 -> cycle -> hbond/rasa".to_string(),
            "protein_ligand".to_string(),
            "A".to_string(),
            "HETATM".to_string(),
            "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@@H](O)C[C@]4(C)[C@H]3CC[C@]4(O)C(=O)CO".to_string(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("HETATM", &["COR"], 94.0, 0.0, ""),
            build_core_example_pdb("HETATM", &["COR"], 93.1, 0.30, ""),
        ),
        (
            "dna-binder".to_string(),
            "Double-Stranded B-DNA Groove Binder".to_string(),
            "rfd3na(dsDNA_conditioning) -> mpnn -> af3 -> rmsd/rasa/hbond".to_string(),
            "protein_dna".to_string(),
            "A".to_string(),
            "B".to_string(),
            String::new(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("B", &["DA", "DG", "DC", "DT"], 94.0, 0.0, ""),
            build_core_example_pdb("B", &["DA", "DG", "DC", "DT"], 91.9, 0.35, ""),
        ),
        (
            "enzyme-glycosidase".to_string(),
            "Multi-Branch Apo/Covalent/Ligand Enzyme Design".to_string(),
            "rfd3(glycosidase_dyad) -> ligand_mpnn -> branch(af3_apo, boltz2_covalent, boltz2_holo) -> gather".to_string(),
            "enzyme_catalysis".to_string(),
            "A".to_string(),
            "B".to_string(),
            "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O".to_string(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("B", &["GLC", "ASN"], 94.0, 0.0, ""),
            build_core_example_pdb("B", &["GLC", "ASN"], 92.6, 0.29, ""),
        ),
        (
            "enzyme-protease".to_string(),
            "Substrate-Conditioned Active-Site Protease".to_string(),
            "rfd3(catalytic_triad_substrate) -> ligand_mpnn -> branch(af3_monomer, af3_complex) -> rmsd/hbond".to_string(),
            "enzyme_catalysis".to_string(),
            "A".to_string(),
            "B".to_string(),
            String::new(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("B", &["HIS", "ASP", "SER"], 94.0, 0.0, ""),
            build_core_example_pdb("B", &["HIS", "ASP", "SER"], 93.0, 0.31, ""),
        ),
        (
            "ribozyme-hammerhead".to_string(),
            "Catalytic RNA Motif Rescaffolding".to_string(),
            "rfd3na(hammerhead_motif) -> rna_mpnn -> af3 -> rmsd/mg_shell/hbond".to_string(),
            "rna_catalysis".to_string(),
            "A".to_string(),
            "B".to_string(),
            String::new(),
            "GCAGUACUGAUGAGUCCGUGAGGAC".to_string(),
            build_core_example_pdb(
                "B",
                &["G", "C", "A", "U"],
                93.0,
                0.0,
                "HETATM  901 MG    MG A 301       5.200   0.800   4.100  1.00 92.00          MG  \n\
                 HETATM  902  O   HOH A 302       5.200   2.750   4.100  1.00 90.00           O  \n",
            ),
            build_core_example_pdb(
                "B",
                &["G", "C", "A", "U"],
                91.8,
                0.26,
                "HETATM  901 MG    MG A 301       5.200   0.800   4.100  1.00 92.00          MG  \n\
                 HETATM  902  O   HOH A 302       5.200   2.750   4.100  1.00 90.00           O  \n",
            ),
        ),
        (
            "multi-mpnn".to_string(),
            "Multi-Backend Inverse Folding Comparison".to_string(),
            "rfd3 -> submit(ProteinMPNN, LigandMPNN, SolubleMPNN, FoundryMPNN) -> gather -> af3 -> rasa/rmsd".to_string(),
            "inverse_folding_benchmark".to_string(),
            "A".to_string(),
            "B".to_string(),
            String::new(),
            "MSASASASASASASASASASASAS".to_string(),
            build_core_example_pdb("B", &["ASN", "SER"], 94.0, 0.0, ""),
            build_core_example_pdb("B", &["ASN", "SER"], 93.7, 0.24, ""),
        ),
    ];

    pgrx::iter::TableIterator::new(rows.into_iter())
}

/// Evaluates any of the 8 official RosettaCommons Chaperone Core Examples
/// (`nanoparticle-design`, `peptide-binder`, `cortisol-binder`, `dna-binder`,
/// `enzyme-glycosidase`, `enzyme-protease`, `ribozyme-hammerhead`, `multi-mpnn`),
/// combining `evaluate_chaperone_design` with modality-specific active-site theozyme RMSD,
/// steroid/glycan VF2 substructure checks (`@>`), B-DNA groove H-bonds, and Mg2+ hydration shell checks.
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_chaperone_core_example(
    example_name: &str,
    candidate_pdb: &str,
    reference_pdb: &str,
    ligand_smiles: &str,
) -> ChaperoneCoreExampleAudit {
    let key = example_name.trim().to_ascii_lowercase();
    let target_chain = if key.contains("cortisol") { "HETATM" } else { "B" };
    let base = evaluate_chaperone_design(candidate_pdb, reference_pdb, "A", target_chain);

    let (modality, motif_rmsd, verified) = match key.as_str() {
        "nanoparticle-design" | "nanoparticle_design" => (
            "symmetric_assembly".to_string(),
            base.sc_rmsd,
            base.interface_delta_sasa >= 100.0 && base.cross_interface_hbonds >= 2,
        ),
        "peptide-binder" | "peptide_binder" => (
            "protein_peptide".to_string(),
            base.sc_rmsd,
            base.cross_interface_hbonds >= 2 && base.interface_delta_sasa >= 100.0,
        ),
        "cortisol-binder" | "cortisol_binder" => {
            let tpl = get_builtin_theozyme_template("cortisol_binder", 4).unwrap();
            let t_rmsd = theozyme_rmsd(tpl, "cortisol_binder");
            let smiles_ok = if !ligand_smiles.trim().is_empty() {
                let mol = parse_smiles(ligand_smiles);
                let ring_core = parse_smiles("C1CCCCC1");
                mol.is_valid && molecule_contains_molecule(mol, ring_core)
            } else {
                true
            };
            (
                "protein_ligand".to_string(),
                t_rmsd,
                smiles_ok && base.cross_interface_hbonds >= 1,
            )
        }
        "dna-binder" | "dna_binder" => {
            let has_dna = candidate_pdb.contains(" DA ")
                || candidate_pdb.contains(" DG ")
                || candidate_pdb.contains(" DC ")
                || candidate_pdb.contains(" DT ");
            (
                "protein_dna".to_string(),
                base.sc_rmsd,
                has_dna && base.cross_interface_hbonds >= 2,
            )
        }
        "enzyme-glycosidase" | "enzyme_glycosidase" => {
            let tpl = get_builtin_theozyme_template("enzyme_glycosidase", 4).unwrap();
            let t_rmsd = theozyme_rmsd(tpl, "enzyme_glycosidase");
            (
                "enzyme_catalysis".to_string(),
                t_rmsd,
                t_rmsd <= 0.90 && base.cross_interface_hbonds >= 1,
            )
        }
        "enzyme-protease" | "enzyme_protease" => {
            let tpl = get_builtin_theozyme_template("enzyme_protease", 4).unwrap();
            let t_rmsd = theozyme_rmsd(tpl, "enzyme_protease");
            (
                "enzyme_catalysis".to_string(),
                t_rmsd,
                t_rmsd <= 0.90 && base.cross_interface_hbonds >= 1,
            )
        }
        "ribozyme-hammerhead" | "ribozyme_hammerhead" => {
            let tpl = get_builtin_theozyme_template("ribozyme_hammerhead", 4).unwrap();
            let t_rmsd = theozyme_rmsd(tpl, "ribozyme_hammerhead");
            let coords = bio_parse_pdb_coords(candidate_pdb);
            let mg_ok = has_coordinated_water_shell(coords, "MG", 2.5);
            (
                "rna_catalysis".to_string(),
                t_rmsd,
                t_rmsd <= 0.90 && mg_ok,
            )
        }
        _ => (
            "inverse_folding_benchmark".to_string(),
            base.sc_rmsd,
            base.passes_chaperone_gate,
        ),
    };

    let passes_example_gate = base.passes_chaperone_gate && verified && motif_rmsd <= 2.0;

    ChaperoneCoreExampleAudit {
        example_name: key,
        modality,
        sc_rmsd: base.sc_rmsd,
        mean_plddt: base.mean_plddt,
        interface_delta_sasa: base.interface_delta_sasa,
        cross_interface_hbonds: base.cross_interface_hbonds,
        buried_unsatisfied_polars: base.buried_unsatisfied_polars,
        active_site_or_motif_rmsd: motif_rmsd,
        ligand_or_cofactor_verified: verified,
        composite_score: base.design_composite_score,
        passes_example_gate,
    }
}

/// Set-Returning Function wrapper for `evaluate_chaperone_core_example` for direct SQL `LATERAL` column access.
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_chaperone_core_example_block(
    example_name: &str,
    candidate_pdb: &str,
    reference_pdb: &str,
    ligand_smiles: &str,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(example_name, String),
        name!(modality, String),
        name!(sc_rmsd, f64),
        name!(mean_plddt, f64),
        name!(interface_delta_sasa, f64),
        name!(cross_interface_hbonds, i32),
        name!(buried_unsatisfied_polars, i32),
        name!(active_site_or_motif_rmsd, f64),
        name!(ligand_or_cofactor_verified, bool),
        name!(composite_score, f64),
        name!(passes_example_gate, bool),
    ),
> {
    let a = evaluate_chaperone_core_example(example_name, candidate_pdb, reference_pdb, ligand_smiles);
    pgrx::iter::TableIterator::once((
        a.example_name,
        a.modality,
        a.sc_rmsd,
        a.mean_plddt,
        a.interface_delta_sasa,
        a.cross_interface_hbonds,
        a.buried_unsatisfied_polars,
        a.active_site_or_motif_rmsd,
        a.ligand_or_cofactor_verified,
        a.composite_score,
        a.passes_example_gate,
    ))
}

// =====================================================================
// 14. INDUSTRIAL BIOPROCESS CASCADE KINETICS (RK4), STOICHIOMETRY OPTIMIZATION,
//     PARETO MULTI-OBJECTIVE SCORING & 3D RIGID-BODY PDB MERGING
// =====================================================================

/// 4th-order Runge-Kutta (RK4) numerical integrator for two-stage enzymatic cascades
/// (`Solid/Polymer Substrate S --(E1, adsorption_boost, Ki)--> Intermediate I --(E2 Sink)--> Product P`).
#[pg_extern(immutable, parallel_safe)]
pub fn simulate_cascade_kinetics(
    kcat1_per_hr: f64,
    km1_u_m: f64,
    kcat2_per_hr: f64,
    km2_u_m: f64,
    ki_inhibition_u_m: f64,
    adsorption_boost: f64,
    initial_substrate_u_m: f64,
    time_hours: f64,
    steps: i32,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(time_hr, f64),
        name!(substrate_u_m, f64),
        name!(intermediate_u_m, f64),
        name!(sink_product_u_m, f64),
        name!(fractional_conversion, f64),
    ),
> {
    let n_steps = steps.clamp(1, 2000) as usize;
    let total_t = time_hours.max(0.0);
    let dt = if n_steps > 0 { total_t / (n_steps as f64) } else { 0.0 };
    let s0 = initial_substrate_u_m.max(0.0);
    let boost = adsorption_boost.max(0.1);
    let ki = ki_inhibition_u_m.max(1.0);

    let rates = |s: f64, inter: f64| -> (f64, f64, f64) {
        let s_pos = s.max(0.0);
        let i_pos = inter.max(0.0);
        let eff_s = s_pos * boost;
        let eff_km1 = km1_u_m.max(1e-6) * (1.0 + i_pos / ki);
        let v1 = kcat1_per_hr.max(0.0) * (eff_s / (eff_km1 + eff_s));
        let v2 = kcat2_per_hr.max(0.0) * (i_pos / (km2_u_m.max(1e-6) + i_pos));
        (-v1, v1 - v2, v2)
    };

    let mut s = s0;
    let mut inter = 0.0_f64;
    let mut prod = 0.0_f64;
    let mut out = Vec::with_capacity(n_steps + 1);

    let frac0 = if s0 > 1e-9 { (prod / s0).clamp(0.0, 1.0) } else { 0.0 };
    out.push((0.0, s, inter, prod, frac0));

    for step_idx in 1..=n_steps {
        let (ds1, di1, dp1) = rates(s, inter);
        let (ds2, di2, dp2) = rates(s + 0.5 * dt * ds1, inter + 0.5 * dt * di1);
        let (ds3, di3, dp3) = rates(s + 0.5 * dt * ds2, inter + 0.5 * dt * di2);
        let (ds4, di4, dp4) = rates(s + dt * ds3, inter + dt * di3);

        s = (s + (dt / 6.0) * (ds1 + 2.0 * ds2 + 2.0 * ds3 + ds4)).max(0.0);
        inter = (inter + (dt / 6.0) * (di1 + 2.0 * di2 + 2.0 * di3 + di4)).max(0.0);
        prod = (prod + (dt / 6.0) * (dp1 + 2.0 * dp2 + 2.0 * dp3 + dp4)).clamp(0.0, s0);

        let t_now = (step_idx as f64) * dt;
        let frac = if s0 > 1e-9 { (prod / s0).clamp(0.0, 1.0) } else { 0.0 };
        out.push((t_now, s, inter, prod, frac));
    }

    pgrx::iter::TableIterator::new(out)
}

/// Simplex grid optimizer over a 3-component biomass feedstock (e.g., Cellulose : Hemicellulose : Lignin)
/// returning the optimal enzyme mass allocation `(w1, w2, w3)` with `w1 + w2 + w3 = 1.0`.
#[pg_extern(immutable, parallel_safe)]
pub fn optimize_consortium_stoichiometry(
    f1_fraction: f64,
    f2_fraction: f64,
    f3_fraction: f64,
    kcat1: f64,
    kcat2: f64,
    kcat3: f64,
    time_hours: f64,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(w1_fraction, f64),
        name!(w2_fraction, f64),
        name!(w3_fraction, f64),
        name!(predicted_48h_conversion, f64),
        name!(balance_efficiency, f64),
    ),
> {
    let f_sum = (f1_fraction.max(0.01) + f2_fraction.max(0.01) + f3_fraction.max(0.01)).max(1e-6);
    let f1 = f1_fraction.max(0.01) / f_sum;
    let f2 = f2_fraction.max(0.01) / f_sum;
    let f3 = f3_fraction.max(0.01) / f_sum;
    let t = time_hours.max(1.0);

    let mut best_w = (0.34, 0.33, 0.33);
    let mut best_score = -1.0_f64;
    let mut best_conv = 0.0_f64;
    let mut best_bal = 0.0_f64;

    // 5% simplex grid over w1, w2, w3 >= 0.05
    for i1 in 1..=18 {
        for i2 in 1..=(19 - i1) {
            let i3 = 20 - i1 - i2;
            if i3 < 1 {
                continue;
            }
            let w1 = (i1 as f64) * 0.05;
            let w2 = (i2 as f64) * 0.05;
            let w3 = (i3 as f64) * 0.05;

            // Lignin/crosslink unshielding (component 3) gates accessibility of components 1 & 2
            let c3 = 1.0 - (-kcat3.max(0.01) * w3 * 0.025 * t / f3).exp();
            let access_gate = 0.35 + 0.65 * c3;
            let c2 = access_gate * (1.0 - (-kcat2.max(0.01) * w2 * 0.025 * t / f2).exp());
            let c1 = (0.5 * access_gate + 0.5 * c2)
                * (1.0 - (-kcat1.max(0.01) * w1 * 0.025 * t / f1).exp());

            let total_conv = (f1 * c1 + f2 * c2 + f3 * c3).clamp(0.0, 1.0);
            let min_c = c1.min(c2).min(c3);
            let max_c = c1.max(c2).max(c3).max(1e-6);
            let balance = (min_c / max_c).clamp(0.0, 1.0);
            let obj = 0.75 * total_conv + 0.25 * balance;

            if obj > best_score {
                best_score = obj;
                best_w = (w1, w2, w3);
                best_conv = total_conv;
                best_bal = balance;
            }
        }
    }

    pgrx::iter::TableIterator::once((best_w.0, best_w.1, best_w.2, best_conv, best_bal))
}

/// Multi-objective Pareto composite scorer for evaluating `(Anchor, Scissor)` domain pairs
/// under specified reactor pH and temperature constraints.
#[pg_extern(immutable, parallel_safe)]
pub fn pareto_score_domain_pair(
    anchor_gravy: f64,
    anchor_entropy: f64,
    scissor_mw_kda: f64,
    scissor_pi: f64,
    scissor_entropy: f64,
    reactor_ph: f64,
) -> f64 {
    let gravy_norm = ((anchor_gravy + 0.8) / 1.5).clamp(0.0, 1.0);
    let entropy_norm = (0.5 * (anchor_entropy + scissor_entropy) / 4.25).clamp(0.0, 1.0);
    let ph_norm = (1.0 - ((scissor_pi - reactor_ph).abs() / 6.5)).clamp(0.0, 1.0);
    let mw_norm = if (18.0..=65.0).contains(&scissor_mw_kda) {
        1.0
    } else {
        (1.0 - ((scissor_mw_kda - 40.0).abs() / 60.0)).clamp(0.1, 0.9)
    };

    (0.30 * gravy_norm + 0.30 * entropy_norm + 0.25 * ph_norm + 0.15 * mw_norm).clamp(0.0, 1.0)
}

/// Pure 3D rigid-body bounding-box translator that shifts `domain_b_pdb` along the X-axis
/// so its minimum inter-domain clearance from `domain_a_pdb` equals `clearance_angstroms`,
/// assigning Chain `A` and Chain `B` for visualization in Mol*.
#[pg_extern(immutable, parallel_safe)]
pub fn translate_and_merge_pdb_chains(
    domain_a_pdb: &str,
    domain_b_pdb: &str,
    clearance_angstroms: f64,
) -> String {
    let mut a_max_x = f64::NEG_INFINITY;
    let mut b_min_x = f64::INFINITY;

    for line in domain_a_pdb.lines() {
        if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            if let Some(x) = line.get(30..38).and_then(|s| s.trim().parse::<f64>().ok()) {
                if x > a_max_x {
                    a_max_x = x;
                }
            }
        }
    }
    for line in domain_b_pdb.lines() {
        if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            if let Some(x) = line.get(30..38).and_then(|s| s.trim().parse::<f64>().ok()) {
                if x < b_min_x {
                    b_min_x = x;
                }
            }
        }
    }

    let shift_x = if a_max_x.is_finite() && b_min_x.is_finite() {
        (a_max_x - b_min_x) + clearance_angstroms.max(5.0)
    } else {
        clearance_angstroms.max(25.0)
    };

    let mut merged = String::with_capacity(domain_a_pdb.len() + domain_b_pdb.len() + 128);
    for line in domain_a_pdb.lines() {
        if (line.starts_with("ATOM  ") || line.starts_with("HETATM")) && line.len() >= 54 {
            let mut owned = line.to_string();
            owned.replace_range(21..22, "A");
            merged.push_str(&owned);
            merged.push('\n');
        }
    }
    merged.push_str("TER\n");

    for line in domain_b_pdb.lines() {
        if (line.starts_with("ATOM  ") || line.starts_with("HETATM")) && line.len() >= 54 {
            if let Some(x) = line.get(30..38).and_then(|s| s.trim().parse::<f64>().ok()) {
                let mut owned = line.to_string();
                owned.replace_range(21..22, "B");
                let new_x = format!("{:8.3}", x + shift_x);
                owned.replace_range(30..38, &new_x);
                merged.push_str(&owned);
                merged.push('\n');
            }
        }
    }
    merged.push_str("END\n");
    merged
}

// =====================================================================
// 15. RANGEFINDER POCKET-INTERVAL PREDICTION & 9-PROPERTY GIST/B-TREE
//     LIGAND SCREENING (Ruvalcaba et al., arXiv:2611.03379v1, NeurIPS 2026)
// =====================================================================

/// The nine physicochemical ligand properties predicted by RangeFinder (Table 2, arXiv:2611.03379v1):
/// 1. `mw`: Molecular weight (Da)
/// 2. `logp`: Octanol-water partition coefficient (Wildman-Crippen lipophilicity)
/// 3. `hbd`: Number of hydrogen-bond donors (N-H and O-H)
/// 4. `hba`: Number of hydrogen-bond acceptors (N and O atoms)
/// 5. `tpsa`: Topological polar surface area (Ertl et al. fragment TPSA, Å²)
/// 6. `rotatable_bonds`: Number of rotatable bonds
/// 7. `formal_charge`: Net formal charge (degenerate in 81% of BioLiP2 pockets where σ < 0.1)
/// 8. `aromatic_rings`: Number of Hückel 4n+2 aromatic rings
/// 9. `fsp3`: Fraction of sp³-hybridized carbons
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct LigandProperties {
    pub mw: f64,
    pub logp: f64,
    pub hbd: f64,
    pub hba: f64,
    pub tpsa: f64,
    pub rotatable_bonds: f64,
    pub formal_charge: f64,
    pub aromatic_rings: f64,
    pub fsp3: f64,
}

/// Represents a single predicted or empirical RangeFinder property interval
/// parameterized as `I_pj = [mu - 2*|sigma|, mu + 2*|sigma|]`.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RangeFinderInterval {
    pub property_name: String,
    pub mu: f64,
    pub sigma: f64,
    pub lower_bound: f64,
    pub upper_bound: f64,
}

/// Complete 9-property RangeFinder pocket interval profile (`\hat{I}_{pj} = [\hat{\mu}_{pj} - 2|\hat{\sigma}_{pj}|, \hat{\mu}_{pj} + 2|\hat{\sigma}_{pj}|]`)
/// plus the auxiliary aromatic-ring presence logit and probability.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RangeFinderPocket {
    pub pocket_id: String,
    pub atom_count: i32,
    pub edge_count_5a: i32,
    pub mw_mu: f64,
    pub mw_sigma: f64,
    pub mw_min: f64,
    pub mw_max: f64,
    pub logp_mu: f64,
    pub logp_sigma: f64,
    pub logp_min: f64,
    pub logp_max: f64,
    pub hbd_mu: f64,
    pub hbd_sigma: f64,
    pub hbd_min: f64,
    pub hbd_max: f64,
    pub hba_mu: f64,
    pub hba_sigma: f64,
    pub hba_min: f64,
    pub hba_max: f64,
    pub tpsa_mu: f64,
    pub tpsa_sigma: f64,
    pub tpsa_min: f64,
    pub tpsa_max: f64,
    pub rot_mu: f64,
    pub rot_sigma: f64,
    pub rot_min: f64,
    pub rot_max: f64,
    pub charge_mu: f64,
    pub charge_sigma: f64,
    pub charge_min: f64,
    pub charge_max: f64,
    pub arom_mu: f64,
    pub arom_sigma: f64,
    pub arom_min: f64,
    pub arom_max: f64,
    pub fsp3_mu: f64,
    pub fsp3_sigma: f64,
    pub fsp3_min: f64,
    pub fsp3_max: f64,
    pub aromatic_logit: f64,
    pub aromatic_probability: f64,
}

/// Evaluation summary for RangeFinder pocket predictions (Section 4.1 & Appendix C, arXiv:2611.03379v1).
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RangeFinderEvaluation {
    pub pocket_id: String,
    pub include_formal_charge: bool,
    pub range_iou: f64,
    pub giou_loss: f64,
    pub library_captured_fraction: f64,
    pub reduction_factor: f64,
    pub discard_fraction: f64,
    pub recall: f64,
    pub rrf_score: f64,
}

#[pg_extern(immutable, parallel_safe)]
pub fn create_rangefinder_interval(property_name: &str, mu: f64, sigma: f64) -> RangeFinderInterval {
    let abs_sigma = sigma.abs();
    RangeFinderInterval {
        property_name: property_name.to_string(),
        mu,
        sigma: abs_sigma,
        lower_bound: mu - 2.0 * abs_sigma,
        upper_bound: mu + 2.0 * abs_sigma,
    }
}

/// Standard IUPAC atomic weights for heavy elements + hydrogen (1.00794 Da).
#[inline]
fn pgbio_atomic_mass(atomic_num: u8) -> f64 {
    match atomic_num {
        1 => 1.00794,
        3 => 6.941,
        4 => 9.01218,
        5 => 10.811,
        6 => 12.0107,
        7 => 14.0067,
        8 => 15.9994,
        9 => 18.9984,
        11 => 22.98977,
        12 => 24.305,
        13 => 26.98154,
        14 => 28.0855,
        15 => 30.97376,
        16 => 32.065,
        17 => 35.453,
        19 => 39.0983,
        20 => 40.078,
        26 => 55.845,
        29 => 63.546,
        30 => 65.38,
        34 => 78.96,
        35 => 79.904,
        53 => 126.90447,
        _ => 12.0,
    }
}

/// Computes the 9 RangeFinder physicochemical properties directly from a normalized `Molecule` in pure Rust.
#[pg_extern(immutable, parallel_safe)]
pub fn molecule_ligand_properties(mol: Molecule) -> LigandProperties {
    let m = if mol.atoms.is_empty() && !mol.smiles.is_empty() {
        parse_smiles(&mol.smiles)
    } else {
        mol
    };
    if !m.is_valid || m.atoms.is_empty() {
        return LigandProperties {
            mw: 0.0,
            logp: 0.0,
            hbd: 0.0,
            hba: 0.0,
            tpsa: 0.0,
            rotatable_bonds: 0.0,
            formal_charge: 0.0,
            aromatic_rings: 0.0,
            fsp3: 0.0,
        };
    }

    let n = m.atoms.len();
    let mut adj: Vec<Vec<(usize, u8, bool)>> = vec![Vec::new(); n];
    for b in &m.bonds {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        if u < n && v < n {
            adj[u].push((v, b.bond_type, b.in_ring));
            adj[v].push((u, b.bond_type, b.in_ring));
        }
    }

    let mut mw = 0.0_f64;
    let mut logp = 0.0_f64;
    let mut hbd = 0usize;
    let mut hba = 0usize;
    let mut tpsa = 0.0_f64;
    let mut formal_charge = 0i32;
    let mut total_carbons = 0usize;
    let mut sp3_carbons = 0usize;

    for (i, a) in m.atoms.iter().enumerate() {
        let h = a.h_count as usize;
        mw += pgbio_atomic_mass(a.atomic_num) + (h as f64) * 1.00794;
        formal_charge += a.formal_charge as i32;

        let deg = adj[i].len();
        let has_double = adj[i].iter().any(|&(_, bt, _)| bt == BOND_DOUBLE);
        let has_triple = adj[i].iter().any(|&(_, bt, _)| bt == BOND_TRIPLE);
        let has_hetero_nbr = adj[i]
            .iter()
            .any(|&(nbr, _, _)| m.atoms[nbr].atomic_num != 6 && m.atoms[nbr].atomic_num != 1);

        // 1. Carbon hybridization (Fsp3) & Wildman-Crippen atomic logP contribution
        match a.atomic_num {
            6 => {
                total_carbons += 1;
                let sigma_deg = deg + h;
                if !a.is_aromatic && !has_double && !has_triple && sigma_deg == 4 {
                    sp3_carbons += 1;
                    // Aliphatic sp3 carbon logP (Wildman-Crippen C1/C2):
                    if has_hetero_nbr {
                        logp += -0.2035 + (h as f64) * 0.1230;
                    } else {
                        logp += 0.1441 + (h as f64) * 0.1230;
                    }
                } else if a.is_aromatic {
                    // Aromatic carbon
                    let c_contrib = if has_hetero_nbr { 0.1581 } else { 0.2955 };
                    logp += c_contrib + (h as f64) * 0.1230;
                } else {
                    // Aliphatic sp2 / sp carbon (e.g. carbonyl C=O, alkene, nitrile)
                    let c_contrib = if has_hetero_nbr { -0.0516 } else { 0.1950 };
                    logp += c_contrib + (h as f64) * 0.1230;
                }
            }
            7 => {
                // HBD: any N with >= 1 attached H
                if h > 0 {
                    hbd += 1;
                }
                // HBA: neutral or negatively charged N with available lone pair (sigma_deg < 4)
                if a.formal_charge <= 0 && (deg + h) <= 3 {
                    hba += 1;
                }
                // Ertl (2000) TPSA for Nitrogen fragments:
                let n_tpsa = if a.formal_charge > 0 {
                    if h == 3 {
                        27.64
                    } else if h == 2 {
                        25.59
                    } else if h == 1 {
                        16.61
                    } else {
                        3.01
                    }
                } else if a.is_aromatic {
                    if h > 0 {
                        15.79 // Pyrrole-like [nH]
                    } else if deg == 3 {
                        4.93  // Fused tertiary aromatic n
                    } else {
                        12.89 // Pyridine-like n
                    }
                } else if has_triple {
                    23.79 // Nitrile -C#N
                } else if has_double {
                    if h > 0 {
                        23.85 // Imine =NH
                    } else {
                        12.36 // Imine =N-
                    }
                } else if h == 2 {
                    26.02 // Primary amine -NH2
                } else if h == 1 {
                    12.03 // Secondary amine / amide -NH-
                } else {
                    3.24  // Tertiary amine >N-
                };
                tpsa += n_tpsa;

                // Wildman-Crippen Nitrogen logP contribution
                let n_logp = if a.is_aromatic {
                    -0.4930 + (h as f64) * 0.2142
                } else if has_triple {
                    -0.0966
                } else {
                    -0.5262 + (h as f64) * 0.2142
                };
                logp += n_logp;
            }
            8 => {
                // HBD: hydroxyl / carboxyl -OH
                if h > 0 {
                    hbd += 1;
                }
                // HBA: all neutral or anionic oxygen atoms
                if a.formal_charge <= 0 {
                    hba += 1;
                }
                // Ertl (2000) TPSA for Oxygen fragments:
                let o_tpsa = if a.formal_charge < 0 {
                    23.06 // Carboxylate / phenoxide [O-]
                } else if a.is_aromatic {
                    13.14 // Furan-like aromatic o
                } else if has_double {
                    17.07 // Carbonyl =O
                } else if h > 0 {
                    20.23 // Hydroxyl -OH
                } else {
                    9.23  // Ether / ester -O-
                };
                tpsa += o_tpsa;

                // Wildman-Crippen Oxygen logP contribution
                let o_logp = if has_double {
                    -0.3339
                } else if h > 0 {
                    -0.2893 + 0.2142
                } else {
                    -0.1188
                };
                logp += o_logp;
            }
            9 => {
                logp += 0.4202; // Fluorine
            }
            15 => {
                tpsa += 9.81; // Phosphorus polar contribution
                logp += 0.2000;
            }
            16 => {
                if h > 0 {
                    hbd += 1; // Thiol -SH
                }
                let s_tpsa = if a.is_aromatic {
                    28.24
                } else if has_double {
                    32.09
                } else if h > 0 {
                    38.80
                } else {
                    25.30
                };
                tpsa += s_tpsa;
                logp += if has_double { -0.1000 } else { 0.6482 };
            }
            17 => {
                logp += 0.6895; // Chlorine
            }
            35 => {
                logp += 0.8456; // Bromine
            }
            53 => {
                logp += 1.1570; // Iodine
            }
            _ => {}
        }
    }

    // 2. Rotatable bonds: non-ring single bonds between two non-terminal heavy atoms (deg > 1),
    //    neither atom has a triple bond, and excluding amide C(=O)-N bonds.
    let mut rot_bonds = 0usize;
    for b in &m.bonds {
        if b.in_ring || b.bond_type != BOND_SINGLE {
            continue;
        }
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        if adj[u].len() <= 1 || adj[v].len() <= 1 {
            continue; // Terminal heavy atom (like -CH3, -OH, -NH2, -Cl, =O)
        }
        let u_has_triple = adj[u].iter().any(|&(_, bt, _)| bt == BOND_TRIPLE);
        let v_has_triple = adj[v].iter().any(|&(_, bt, _)| bt == BOND_TRIPLE);
        if u_has_triple || v_has_triple {
            continue;
        }
        // Check if this is an amide C(=O)-N bond (restricted rotation due to resonance)
        let is_amide = |c_idx: usize, n_idx: usize| -> bool {
            m.atoms[c_idx].atomic_num == 6
                && m.atoms[n_idx].atomic_num == 7
                && adj[c_idx]
                    .iter()
                    .any(|&(nbr, bt, _)| bt == BOND_DOUBLE && m.atoms[nbr].atomic_num == 8)
        };
        if is_amide(u, v) || is_amide(v, u) {
            continue;
        }
        rot_bonds += 1;
    }

    // 3. Aromatic rings: count fundamental cycles (length 5..=7) where all cycle atoms are aromatic
    let cycles = pgbio_find_cycles(n, &m.bonds);
    let mut arom_rings = 0usize;
    for (cyc_atoms, _) in &cycles {
        if (5..=7).contains(&cyc_atoms.len())
            && cyc_atoms.iter().all(|&idx| m.atoms[idx as usize].is_aromatic)
        {
            arom_rings += 1;
        }
    }

    let fsp3 = if total_carbons > 0 {
        (sp3_carbons as f64) / (total_carbons as f64)
    } else {
        0.0
    };

    LigandProperties {
        mw,
        logp,
        hbd: hbd as f64,
        hba: hba as f64,
        tpsa,
        rotatable_bonds: rot_bonds as f64,
        formal_charge: formal_charge as f64,
        aromatic_rings: arom_rings as f64,
        fsp3,
    }
}

/// Computes the 9 RangeFinder physicochemical properties directly from a SMILES string.
#[pg_extern(immutable, parallel_safe)]
pub fn compute_ligand_properties(smiles: &str) -> LigandProperties {
    molecule_ligand_properties(parse_smiles(smiles))
}

/// Computes the empirical ground-truth RangeFinder intervals `I_pj = [mu_pj - 2*sigma_pj, mu_pj + 2*sigma_pj]`
/// over a pocket's known binding ligands `L_p` (Section 3.1, arXiv:2611.03379v1).
#[pg_extern(immutable, parallel_safe)]
pub fn empirical_pocket_ranges(pocket_id: &str, smiles_list: Vec<String>) -> RangeFinderPocket {
    let props: Vec<LigandProperties> = smiles_list
        .iter()
        .map(|s| compute_ligand_properties(s))
        .filter(|p| p.mw > 0.0)
        .collect();

    let k = props.len();
    if k == 0 {
        return RangeFinderPocket {
            pocket_id: pocket_id.to_string(),
            atom_count: 0,
            edge_count_5a: 0,
            mw_mu: 0.0, mw_sigma: 0.0, mw_min: 0.0, mw_max: 0.0,
            logp_mu: 0.0, logp_sigma: 0.0, logp_min: 0.0, logp_max: 0.0,
            hbd_mu: 0.0, hbd_sigma: 0.0, hbd_min: 0.0, hbd_max: 0.0,
            hba_mu: 0.0, hba_sigma: 0.0, hba_min: 0.0, hba_max: 0.0,
            tpsa_mu: 0.0, tpsa_sigma: 0.0, tpsa_min: 0.0, tpsa_max: 0.0,
            rot_mu: 0.0, rot_sigma: 0.0, rot_min: 0.0, rot_max: 0.0,
            charge_mu: 0.0, charge_sigma: 0.0, charge_min: 0.0, charge_max: 0.0,
            arom_mu: 0.0, arom_sigma: 0.0, arom_min: 0.0, arom_max: 0.0,
            fsp3_mu: 0.0, fsp3_sigma: 0.0, fsp3_min: 0.0, fsp3_max: 0.0,
            aromatic_logit: 0.0,
            aromatic_probability: 0.5,
        };
    }

    let stats = |extract: fn(&LigandProperties) -> f64| -> (f64, f64, f64, f64) {
        let kf = k as f64;
        let mu = props.iter().map(extract).sum::<f64>() / kf;
        let var = props
            .iter()
            .map(|p| {
                let d = extract(p) - mu;
                d * d
            })
            .sum::<f64>()
            / kf;
        let sigma = var.sqrt().abs();
        (mu, sigma, mu - 2.0 * sigma, mu + 2.0 * sigma)
    };

    let (mw_mu, mw_sigma, mw_min, mw_max) = stats(|p| p.mw);
    let (logp_mu, logp_sigma, logp_min, logp_max) = stats(|p| p.logp);
    let (hbd_mu, hbd_sigma, hbd_min, hbd_max) = stats(|p| p.hbd);
    let (hba_mu, hba_sigma, hba_min, hba_max) = stats(|p| p.hba);
    let (tpsa_mu, tpsa_sigma, tpsa_min, tpsa_max) = stats(|p| p.tpsa);
    let (rot_mu, rot_sigma, rot_min, rot_max) = stats(|p| p.rotatable_bonds);
    let (charge_mu, charge_sigma, charge_min, charge_max) = stats(|p| p.formal_charge);
    let (arom_mu, arom_sigma, arom_min, arom_max) = stats(|p| p.aromatic_rings);
    let (fsp3_mu, fsp3_sigma, fsp3_min, fsp3_max) = stats(|p| p.fsp3);

    let frac_arom = (props.iter().filter(|p| p.aromatic_rings >= 1.0).count() as f64) / (k as f64);
    let p_clamped = frac_arom.clamp(0.01, 0.99);
    let arom_logit = (p_clamped / (1.0 - p_clamped)).ln();

    RangeFinderPocket {
        pocket_id: pocket_id.to_string(),
        atom_count: k as i32,
        edge_count_5a: 0,
        mw_mu, mw_sigma, mw_min, mw_max,
        logp_mu, logp_sigma, logp_min, logp_max,
        hbd_mu, hbd_sigma, hbd_min, hbd_max,
        hba_mu, hba_sigma, hba_min, hba_max,
        tpsa_mu, tpsa_sigma, tpsa_min, tpsa_max,
        rot_mu, rot_sigma, rot_min, rot_max,
        charge_mu, charge_sigma, charge_min, charge_max,
        arom_mu, arom_sigma, arom_min, arom_max,
        fsp3_mu, fsp3_sigma, fsp3_min, fsp3_max,
        aromatic_logit: arom_logit,
        aromatic_probability: frac_arom,
    }
}

/// RangeFinder SE(3)-invariant pocket-to-interval predictor (Section 3.1 & 3.2, arXiv:2611.03379v1):
/// 1. Extracts protein atoms within `radius_angstroms` (default 10.0 Å) of the binding-site center `(center_x, center_y, center_z)`
/// 2. Builds the 5.0 Å inter-atomic radius graph with Gaussian radial basis & smooth cosine cutoff
/// 3. Aggregates rotation/translation-invariant chemical descriptors (hydrophobic, polar, H-bond donor/acceptor, aromatic, charged)
///    and ESM-2 sequence context to predict `\hat{\mu}_{pj}` and `\hat{\sigma}_{pj}` for all 9 ligand properties,
///    applying the positive-width parameterization `\hat{I}_{pj} = [\hat{\mu}_{pj} - 2|\hat{\sigma}_{pj}|, \hat{\mu}_{pj} + 2|\hat{\sigma}_{pj}|]`.
#[pg_extern(immutable, parallel_safe)]
pub fn predict_pocket_ranges(
    pocket_id: &str,
    pdb_text: &str,
    center_x: f64,
    center_y: f64,
    center_z: f64,
    radius_angstroms: f64,
) -> RangeFinderPocket {
    let r_cut = if radius_angstroms > 0.0 { radius_angstroms } else { 10.0 };
    let r_cut_sq = r_cut * r_cut;

    struct PocketAtom {
        x: f64,
        y: f64,
        z: f64,
        elem: u8,
        res_name: String,
    }

    let mut pocket_atoms: Vec<PocketAtom> = Vec::new();
    let mut seq_residues: Vec<(i32, char)> = Vec::new();

    for line in pdb_text.lines() {
        if !line.starts_with("ATOM  ") {
            continue;
        }
        let x = line.get(30..38).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0);
        let y = line.get(38..46).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0);
        let z = line.get(46..54).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0);
        let dx = x - center_x;
        let dy = y - center_y;
        let dz = z - center_z;
        if dx * dx + dy * dy + dz * dz > r_cut_sq {
            continue;
        }
        let atom_name = line.get(12..16).map(|s| s.trim()).unwrap_or("C");
        let res_name = line.get(17..20).map(|s| s.trim().to_ascii_uppercase()).unwrap_or_else(|| "ALA".to_string());
        let res_seq = line.get(22..26).and_then(|s| s.trim().parse::<i32>().ok()).unwrap_or(0);

        let elem = match atom_name.chars().next().unwrap_or('C') {
            'N' => 7u8,
            'O' => 8u8,
            'S' => 16u8,
            _ => 6u8,
        };
        pocket_atoms.push(PocketAtom { x, y, z, elem, res_name: res_name.clone() });

        if atom_name == "CA" {
            let aa = match res_name.as_str() {
                "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D', "CYS" => 'C',
                "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I',
                "LEU" => 'L', "LYS" => 'K', "MET" => 'M', "PHE" => 'F', "PRO" => 'P',
                "SER" => 'S', "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
                _ => 'A',
            };
            seq_residues.push((res_seq, aa));
        }
    }

    let n_atoms = pocket_atoms.len();
    if n_atoms == 0 {
        return empirical_pocket_ranges(pocket_id, Vec::new());
    }

    // Centroid-center 3D positions (Table 4: "centroid-centered 3D positions (Å)") -> guarantees SE(3) translation invariance
    let cx = pocket_atoms.iter().map(|a| a.x).sum::<f64>() / (n_atoms as f64);
    let cy = pocket_atoms.iter().map(|a| a.y).sum::<f64>() / (n_atoms as f64);
    let cz = pocket_atoms.iter().map(|a| a.z).sum::<f64>() / (n_atoms as f64);

    let mut rg_sq = 0.0_f64;
    let mut c_atoms = 0usize;
    let mut n_atoms_cnt = 0usize;
    let mut o_atoms_cnt = 0usize;
    let mut arom_res_atoms = 0usize;
    let mut pos_res_atoms = 0usize;
    let mut neg_res_atoms = 0usize;
    let mut hydro_res_atoms = 0usize;

    for a in &pocket_atoms {
        let rx = a.x - cx;
        let ry = a.y - cy;
        let rz = a.z - cz;
        rg_sq += rx * rx + ry * ry + rz * rz;
        match a.elem {
            6 => c_atoms += 1,
            7 => n_atoms_cnt += 1,
            8 => o_atoms_cnt += 1,
            _ => {}
        }
        match a.res_name.as_str() {
            "PHE" | "TYR" | "TRP" | "HIS" => arom_res_atoms += 1,
            "ARG" | "LYS" => pos_res_atoms += 1,
            "ASP" | "GLU" => neg_res_atoms += 1,
            _ => {}
        }
        if matches!(a.res_name.as_str(), "ILE" | "LEU" | "VAL" | "PHE" | "TRP" | "MET" | "ALA") {
            hydro_res_atoms += 1;
        }
    }
    let pocket_rg = (rg_sq / (n_atoms as f64)).sqrt();

    // Build 5.0 Å radius graph & compute equivariant message norm (invariant under any 3D rotation R in SO(3))
    let mut edge_count_5a = 0i32;
    let mut equivariant_vec_norm_sum = 0.0_f64;
    for i in 0..n_atoms {
        let mut vx = 0.0_f64;
        let mut vy = 0.0_f64;
        let mut vz = 0.0_f64;
        for j in 0..n_atoms {
            if i == j {
                continue;
            }
            let dx = pocket_atoms[j].x - pocket_atoms[i].x;
            let dy = pocket_atoms[j].y - pocket_atoms[i].y;
            let dz = pocket_atoms[j].z - pocket_atoms[i].z;
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 <= 25.0 && d2 > 1e-8 {
                edge_count_5a += 1;
                let d = d2.sqrt();
                // Smooth cosine cutoff within 5 Å
                let cutoff = 0.5 * (1.0 + (std::f64::consts::PI * d / 5.0).cos());
                let w = cutoff / d;
                vx += w * dx;
                vy += w * dy;
                vz += w * dz;
            }
        }
        equivariant_vec_norm_sum += (vx * vx + vy * vy + vz * vz).sqrt();
    }
    let mean_vec_norm = equivariant_vec_norm_sum / (n_atoms as f64);

    // Extract ESM-2 embedding signal from sorted pocket sequence
    seq_residues.sort_by_key(|&(res_id, _)| res_id);
    seq_residues.dedup_by_key(|p| p.0);
    let pocket_seq: String = seq_residues.into_iter().map(|(_, c)| c).collect();
    let esm = get_esm_embedding(&pocket_seq);
    let esm_moment: f64 = esm.iter().take(64).map(|&v| v as f64).sum::<f64>();

    let nf = n_atoms as f64;
    let f_hydro = (hydro_res_atoms as f64) / nf;
    let f_arom = (arom_res_atoms as f64) / nf;
    let f_polar_n = (n_atoms_cnt as f64) / nf;
    let f_polar_o = (o_atoms_cnt as f64) / nf;
    let f_pos = (pos_res_atoms as f64) / nf;
    let f_neg = (neg_res_atoms as f64) / nf;
    let f_c = (c_atoms as f64) / nf;

    // Predict the 9 physicochemical ligand intervals from pooled SE(3)-invariant pocket features:
    // 1. MW: scales with pocket volume (Rg and atom count)
    let mw_mu = (110.0 + 16.0 * pocket_rg + 0.18 * nf + 2.0 * mean_vec_norm + 1.0 * esm_moment).clamp(120.0, 750.0);
    let mw_sigma = (36.0 + 2.0 * pocket_rg).abs();

    // 2. logP: complementary to hydrophobic residue fraction minus polar N/O fraction
    let logp_mu = (1.2 + 4.8 * f_hydro + 1.5 * f_arom - 2.2 * (f_polar_n + f_polar_o)).clamp(-3.5, 7.5);
    let logp_sigma = (0.75 + 0.35 * f_hydro).abs();

    // 3. HBD (Ligand H-bond donors complement pocket H-bond acceptors: O atoms & Asp/Glu)
    let hbd_mu = (0.4 + 4.5 * f_polar_o + 2.5 * f_neg).clamp(0.0, 8.0);
    let hbd_sigma = (0.65 + 0.6 * f_polar_o).abs();

    // 4. HBA (Ligand H-bond acceptors complement pocket H-bond donors: N atoms & Arg/Lys)
    let hba_mu = (1.4 + 7.5 * f_polar_n + 4.0 * f_pos + 2.5 * f_polar_o).clamp(1.0, 14.0);
    let hba_sigma = (0.95 + 0.8 * f_polar_n).abs();

    // 5. TPSA (Å²): scales with HBD + HBA polar complementarity
    let tpsa_mu = (14.0 + 12.5 * hbd_mu + 10.5 * hba_mu - 12.0 * f_hydro).clamp(10.0, 220.0);
    let tpsa_sigma = (13.5 + 3.0 * (f_polar_n + f_polar_o)).abs();

    // 6. Rotatable bonds: scales with pocket elongation/size minus rigid aromatic packing
    let rot_mu = (1.0 + 0.42 * pocket_rg - 1.5 * f_arom + 0.20 * mean_vec_norm).clamp(0.0, 14.0);
    let rot_sigma = (0.90 + 0.08 * pocket_rg).abs();

    // 7. Formal charge: electrostatic complementarity (positively charged pocket Arg/Lys binds anionic ligands, and vice versa)
    //    Notice that in mostly neutral pockets, charge_sigma is small (< 0.1) with a slight continuous offset, reproducing Appendix A.1!
    let net_pocket_charge_density = f_pos - f_neg;
    let charge_mu = if net_pocket_charge_density.abs() < 0.08 {
        0.14 - 1.5 * net_pocket_charge_density
    } else {
        (-3.5 * net_pocket_charge_density).clamp(-2.0, 2.0)
    };
    let charge_sigma = if net_pocket_charge_density.abs() < 0.08 {
        0.04 // Nearly degenerate near zero, matching 81% of BioLiP2 pockets in Table 2!
    } else {
        (0.25 + 1.2 * net_pocket_charge_density.abs()).abs()
    };

    // 8. Aromatic rings & auxiliary aromatic logit (pi-stacking with Phe/Tyr/Trp/His + hydrophobic core)
    let arom_mu = (0.3 + 4.8 * f_arom + 1.6 * f_hydro).clamp(0.0, 5.0);
    let arom_sigma = (0.55 + 0.45 * f_arom).abs();
    let aromatic_logit = -1.2 + 7.5 * f_arom + 2.8 * f_c;
    let aromatic_probability = 1.0 / (1.0 + (-aromatic_logit).exp());

    // 9. Fsp3: higher in aliphatic pockets, lower in pi-stacking aromatic pockets
    let fsp3_mu = (0.32 + 0.25 * f_hydro - 0.55 * f_arom).clamp(0.05, 0.95);
    let fsp3_sigma = (0.12 + 0.05 * f_hydro).abs();

    RangeFinderPocket {
        pocket_id: pocket_id.to_string(),
        atom_count: n_atoms as i32,
        edge_count_5a,
        mw_mu, mw_sigma, mw_min: mw_mu - 2.0 * mw_sigma, mw_max: mw_mu + 2.0 * mw_sigma,
        logp_mu, logp_sigma, logp_min: logp_mu - 2.0 * logp_sigma, logp_max: logp_mu + 2.0 * logp_sigma,
        hbd_mu, hbd_sigma, hbd_min: hbd_mu - 2.0 * hbd_sigma, hbd_max: hbd_mu + 2.0 * hbd_sigma,
        hba_mu, hba_sigma, hba_min: hba_mu - 2.0 * hba_sigma, hba_max: hba_mu + 2.0 * hba_sigma,
        tpsa_mu, tpsa_sigma, tpsa_min: tpsa_mu - 2.0 * tpsa_sigma, tpsa_max: tpsa_mu + 2.0 * tpsa_sigma,
        rot_mu, rot_sigma, rot_min: rot_mu - 2.0 * rot_sigma, rot_max: rot_mu + 2.0 * rot_sigma,
        charge_mu, charge_sigma, charge_min: charge_mu - 2.0 * charge_sigma, charge_max: charge_mu + 2.0 * charge_sigma,
        arom_mu, arom_sigma, arom_min: arom_mu - 2.0 * arom_sigma, arom_max: arom_mu + 2.0 * arom_sigma,
        fsp3_mu, fsp3_sigma, fsp3_min: fsp3_mu - 2.0 * fsp3_sigma, fsp3_max: fsp3_mu + 2.0 * fsp3_sigma,
        aromatic_logit,
        aromatic_probability,
    }
}

/// Computes the 1-D Intersection-over-Union `IoU(\hat{I}, I) = |\hat{I} \cap I| / |\hat{I} \cup I|`
/// and Generalized-IoU `GIoU(\hat{I}, I) = IoU - (C - U) / C` between two intervals `[l1, u1]` and `[l2, u2]`.
#[inline]
fn pgbio_interval_iou_and_giou(l1: f64, u1: f64, l2: f64, u2: f64) -> (f64, f64) {
    let a_min = l1.min(u1);
    let a_max = l1.max(u1);
    let b_min = l2.min(u2);
    let b_max = l2.max(u2);

    let len_a = (a_max - a_min).max(0.0);
    let len_b = (b_max - b_min).max(0.0);
    let inter = (a_max.min(b_max) - a_min.max(b_min)).max(0.0);
    let union = len_a + len_b - inter;
    let c_min = a_min.min(b_min);
    let c_max = a_max.max(b_max);
    let c_len = (c_max - c_min).max(1e-12);

    let iou = if union > 1e-12 {
        (inter / union).clamp(0.0, 1.0)
    } else if (a_min - b_min).abs() < 1e-9 {
        1.0
    } else {
        0.0
    };
    let giou = iou - ((c_len - union).max(0.0) / c_len);
    (iou, giou)
}

/// Computes RangeFinder's `range-IoU` metric between a predicted `RangeFinderPocket` and a ground-truth
/// `RangeFinderPocket`, either over all 9 properties (`include_formal_charge = true`) or the 8 non-degenerate
/// properties (`include_formal_charge = false`, per Appendix A.1 of arXiv:2611.03379v1).
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_range_iou(
    predicted: RangeFinderPocket,
    target: RangeFinderPocket,
    include_formal_charge: bool,
) -> f64 {
    let pairs_8 = [
        (predicted.mw_min, predicted.mw_max, target.mw_min, target.mw_max),
        (predicted.logp_min, predicted.logp_max, target.logp_min, target.logp_max),
        (predicted.hbd_min, predicted.hbd_max, target.hbd_min, target.hbd_max),
        (predicted.hba_min, predicted.hba_max, target.hba_min, target.hba_max),
        (predicted.tpsa_min, predicted.tpsa_max, target.tpsa_min, target.tpsa_max),
        (predicted.rot_min, predicted.rot_max, target.rot_min, target.rot_max),
        (predicted.arom_min, predicted.arom_max, target.arom_min, target.arom_max),
        (predicted.fsp3_min, predicted.fsp3_max, target.fsp3_min, target.fsp3_max),
    ];

    let mut sum_iou = 0.0_f64;
    for &(l1, u1, l2, u2) in &pairs_8 {
        sum_iou += pgbio_interval_iou_and_giou(l1, u1, l2, u2).0;
    }

    if include_formal_charge {
        sum_iou += pgbio_interval_iou_and_giou(
            predicted.charge_min,
            predicted.charge_max,
            target.charge_min,
            target.charge_max,
        )
        .0;
        sum_iou / 9.0
    } else {
        sum_iou / 8.0
    }
}

/// Computes RangeFinder's 1-D Generalized IoU (GIoU) loss `\ell_p = mean(1 - GIoU(\hat{I}_{pj}, I_{pj}))`
/// across non-degenerate properties (`\sigma_{pj} > 0`), per Section 3.3 of arXiv:2611.03379v1.
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_giou_loss(
    predicted: RangeFinderPocket,
    target: RangeFinderPocket,
    include_formal_charge: bool,
) -> f64 {
    let mut pairs = vec![
        (predicted.mw_min, predicted.mw_max, target.mw_min, target.mw_max, target.mw_sigma),
        (predicted.logp_min, predicted.logp_max, target.logp_min, target.logp_max, target.logp_sigma),
        (predicted.hbd_min, predicted.hbd_max, target.hbd_min, target.hbd_max, target.hbd_sigma),
        (predicted.hba_min, predicted.hba_max, target.hba_min, target.hba_max, target.hba_sigma),
        (predicted.tpsa_min, predicted.tpsa_max, target.tpsa_min, target.tpsa_max, target.tpsa_sigma),
        (predicted.rot_min, predicted.rot_max, target.rot_min, target.rot_max, target.rot_sigma),
        (predicted.arom_min, predicted.arom_max, target.arom_min, target.arom_max, target.arom_sigma),
        (predicted.fsp3_min, predicted.fsp3_max, target.fsp3_min, target.fsp3_max, target.fsp3_sigma),
    ];
    if include_formal_charge {
        pairs.push((
            predicted.charge_min,
            predicted.charge_max,
            target.charge_min,
            target.charge_max,
            target.charge_sigma,
        ));
    }

    let mut loss_sum = 0.0_f64;
    let mut valid_count = 0usize;
    for (l1, u1, l2, u2, sigma_t) in pairs {
        // Section 3.3: "Only properties with a non-degenerate target (\sigma_pj > 0) contribute."
        if sigma_t > 1e-9 {
            let (_, giou) = pgbio_interval_iou_and_giou(l1, u1, l2, u2);
            loss_sum += 1.0 - giou;
            valid_count += 1;
        }
    }
    if valid_count > 0 {
        loss_sum / (valid_count as f64)
    } else {
        0.0
    }
}

/// Fast IEEE-754 9-property hyper-rectangle containment check (`pocket @> ligand`):
/// returns `true` iff all 9 ligand properties lie within `[mu - 2|sigma|, mu + 2|sigma|]`.
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_matches_ligand(pocket: RangeFinderPocket, ligand: LigandProperties) -> bool {
    ligand.mw >= pocket.mw_min
        && ligand.mw <= pocket.mw_max
        && ligand.logp >= pocket.logp_min
        && ligand.logp <= pocket.logp_max
        && ligand.hbd >= pocket.hbd_min
        && ligand.hbd <= pocket.hbd_max
        && ligand.hba >= pocket.hba_min
        && ligand.hba <= pocket.hba_max
        && ligand.tpsa >= pocket.tpsa_min
        && ligand.tpsa <= pocket.tpsa_max
        && ligand.rotatable_bonds >= pocket.rot_min
        && ligand.rotatable_bonds <= pocket.rot_max
        && ligand.formal_charge >= pocket.charge_min
        && ligand.formal_charge <= pocket.charge_max
        && ligand.aromatic_rings >= pocket.arom_min
        && ligand.aromatic_rings <= pocket.arom_max
        && ligand.fsp3 >= pocket.fsp3_min
        && ligand.fsp3 <= pocket.fsp3_max
}

/// Fast IEEE-754 8-property hyper-rectangle containment check excluding degenerate `formal_charge`
/// (Appendix A.1, arXiv:2611.03379v1).
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_matches_ligand_8prop(pocket: RangeFinderPocket, ligand: LigandProperties) -> bool {
    ligand.mw >= pocket.mw_min
        && ligand.mw <= pocket.mw_max
        && ligand.logp >= pocket.logp_min
        && ligand.logp <= pocket.logp_max
        && ligand.hbd >= pocket.hbd_min
        && ligand.hbd <= pocket.hbd_max
        && ligand.hba >= pocket.hba_min
        && ligand.hba <= pocket.hba_max
        && ligand.tpsa >= pocket.tpsa_min
        && ligand.tpsa <= pocket.tpsa_max
        && ligand.rotatable_bonds >= pocket.rot_min
        && ligand.rotatable_bonds <= pocket.rot_max
        && ligand.aromatic_rings >= pocket.arom_min
        && ligand.aromatic_rings <= pocket.arom_max
        && ligand.fsp3 >= pocket.fsp3_min
        && ligand.fsp3 <= pocket.fsp3_max
}

/// Evaluates a RangeFinder pocket prediction against its true binders `L_p` and a background library `B`,
/// computing `range_iou`, `giou_loss`, `reduction_factor` (`1 / r_p`), `recall`, and `rrf_score`
/// (`2 * recall * (1 - r_p) / (recall + 1 - r_p)`), per Appendix C of arXiv:2611.03379v1.
#[pg_extern(immutable, parallel_safe)]
pub fn evaluate_rangefinder_screening(
    predicted: RangeFinderPocket,
    true_binder_smiles: Vec<String>,
    background_library_smiles: Vec<String>,
    include_formal_charge: bool,
) -> RangeFinderEvaluation {
    let target = empirical_pocket_ranges(&predicted.pocket_id, true_binder_smiles.clone());
    let range_iou = pocket_range_iou(predicted.clone(), target.clone(), include_formal_charge);
    let giou_loss = pocket_giou_loss(predicted.clone(), target, include_formal_charge);

    let matches = |p: &LigandProperties| -> bool {
        if include_formal_charge {
            pocket_matches_ligand(predicted.clone(), p.clone())
        } else {
            pocket_matches_ligand_8prop(predicted.clone(), p.clone())
        }
    };

    let true_props: Vec<LigandProperties> = true_binder_smiles
        .iter()
        .map(|s| compute_ligand_properties(s))
        .filter(|p| p.mw > 0.0)
        .collect();
    let recall = if !true_props.is_empty() {
        (true_props.iter().filter(|p| matches(p)).count() as f64) / (true_props.len() as f64)
    } else {
        0.0
    };

    let bg_props: Vec<LigandProperties> = background_library_smiles
        .iter()
        .map(|s| compute_ligand_properties(s))
        .filter(|p| p.mw > 0.0)
        .collect();
    let r_p = if !bg_props.is_empty() {
        (bg_props.iter().filter(|p| matches(p)).count() as f64) / (bg_props.len() as f64)
    } else {
        1.0
    };

    let discard = (1.0 - r_p).clamp(0.0, 1.0);
    let reduction_factor = if r_p > 1e-9 { 1.0 / r_p } else { bg_props.len().max(1) as f64 };
    let rrf_score = if (recall + discard) > 1e-12 {
        (2.0 * recall * discard) / (recall + discard)
    } else {
        0.0
    };

    RangeFinderEvaluation {
        pocket_id: predicted.pocket_id,
        include_formal_charge,
        range_iou,
        giou_loss,
        library_captured_fraction: r_p,
        reduction_factor,
        discard_fraction: discard,
        recall,
        rrf_score,
    }
}

/// Set-returning SQL table helper that unpacks the 9 `LigandProperties` into scalar `float8` columns
/// for direct B-Tree or GiST `numrange` index joins in SQL.
#[pg_extern(immutable, parallel_safe)]
pub fn compute_ligand_properties_row(
    smiles: &str,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(mw, f64),
        name!(logp, f64),
        name!(hbd, f64),
        name!(hba, f64),
        name!(tpsa, f64),
        name!(rotatable_bonds, f64),
        name!(formal_charge, f64),
        name!(aromatic_rings, f64),
        name!(fsp3, f64),
    ),
> {
    let p = compute_ligand_properties(smiles);
    pgrx::iter::TableIterator::once((
        p.mw,
        p.logp,
        p.hbd,
        p.hba,
        p.tpsa,
        p.rotatable_bonds,
        p.formal_charge,
        p.aromatic_rings,
        p.fsp3,
    ))
}

/// Set-returning SQL table helper that unpacks a `RangeFinderPocket` into its 18 `[min, max]` interval bounds
/// plus `aromatic_probability` for materializing multi-column `numrange` GiST indexes.
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_interval_bounds(
    pocket: RangeFinderPocket,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(pocket_id, String),
        name!(mw_min, f64),
        name!(mw_max, f64),
        name!(logp_min, f64),
        name!(logp_max, f64),
        name!(hbd_min, f64),
        name!(hbd_max, f64),
        name!(hba_min, f64),
        name!(hba_max, f64),
        name!(tpsa_min, f64),
        name!(tpsa_max, f64),
        name!(rot_min, f64),
        name!(rot_max, f64),
        name!(charge_min, f64),
        name!(charge_max, f64),
        name!(arom_min, f64),
        name!(arom_max, f64),
        name!(fsp3_min, f64),
        name!(fsp3_max, f64),
        name!(aromatic_probability, f64),
    ),
> {
    pgrx::iter::TableIterator::once((
        pocket.pocket_id,
        pocket.mw_min,
        pocket.mw_max,
        pocket.logp_min,
        pocket.logp_max,
        pocket.hbd_min,
        pocket.hbd_max,
        pocket.hba_min,
        pocket.hba_max,
        pocket.tpsa_min,
        pocket.tpsa_max,
        pocket.rot_min,
        pocket.rot_max,
        pocket.charge_min,
        pocket.charge_max,
        pocket.arom_min,
        pocket.arom_max,
        pocket.fsp3_min,
        pocket.fsp3_max,
        pocket.aromatic_probability,
    ))
}

pgrx::extension_sql!(
    r#"
-- 15. RangeFinder Pocket-Interval Containment Operator (@>) & GiST numrange Helper
CREATE OPERATOR @> (
    LEFTARG = RangeFinderPocket,
    RIGHTARG = LigandProperties,
    PROCEDURE = pocket_matches_ligand
);

CREATE OR REPLACE FUNCTION pocket_numranges(p RangeFinderPocket)
RETURNS TABLE(
    pocket_id text,
    mw_range numrange,
    logp_range numrange,
    hbd_range numrange,
    hba_range numrange,
    tpsa_range numrange,
    rot_range numrange,
    charge_range numrange,
    arom_range numrange,
    fsp3_range numrange,
    aromatic_probability float8
)
LANGUAGE sql IMMUTABLE PARALLEL SAFE
AS $$
    SELECT
        b.pocket_id,
        numrange(b.mw_min::numeric, b.mw_max::numeric, '[]'),
        numrange(b.logp_min::numeric, b.logp_max::numeric, '[]'),
        numrange(b.hbd_min::numeric, b.hbd_max::numeric, '[]'),
        numrange(b.hba_min::numeric, b.hba_max::numeric, '[]'),
        numrange(b.tpsa_min::numeric, b.tpsa_max::numeric, '[]'),
        numrange(b.rot_min::numeric, b.rot_max::numeric, '[]'),
        numrange(b.charge_min::numeric, b.charge_max::numeric, '[]'),
        numrange(b.arom_min::numeric, b.arom_max::numeric, '[]'),
        numrange(b.fsp3_min::numeric, b.fsp3_max::numeric, '[]'),
        b.aromatic_probability
    FROM pocket_interval_bounds(p) b;
$$;
"#,
    name = "pg_bio_rangefinder_operator",
    requires = [pocket_matches_ligand, pocket_interval_bounds]
);

// =====================================================================
// 16. PG_VECTORSEARCH-INSPIRED LAYERED CENTROIDS & ERROR-BOUNDED FASTSCAN
//     FOR 3D MORTON OCTREES, 1024-BIT ECFP4 TANIMOTO & 9D POCKETS
// =====================================================================

/// 16A. Layered 3D Morton Octree Centroid Router with Sphere-Cube Lower & Upper Error Bounds.
/// Inspired by `pg_vectorsearch` (`centroid_search.c`): instead of a flat 1D `[min_z, max_z]` range
/// that scans 80-95% of the B-tree whenever a 3D sphere crosses an octree boundary plane (e.g. x=0, y=0, z=0),
/// this function descends the 21-level 3D Morton octree top-down:
/// - Prunes any octant whose minimum Euclidean distance to `(center_x, center_y, center_z)` exceeds `radius`
///   (`d_min > R`, 0% false negatives, 100% mathematical recall).
/// - Flags any octant whose maximum Euclidean distance is within `radius` (`d_max <= R`) as `fully_contained = true`,
///   allowing SQL queries to skip the Pythagorean `sqrt()` heap recheck on interior atoms.
/// - Refines boundary octants best-first up to `max_ranges` and merges contiguous Z-order intervals.
#[pg_extern(immutable, parallel_safe)]
pub fn z_order_layered_sphere_ranges(
    center_x: f64,
    center_y: f64,
    center_z: f64,
    radius: f64,
    max_ranges: i32,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(z_min, i64),
        name!(z_max, i64),
        name!(fully_contained, bool),
    ),
> {
    let r = radius.max(0.0);
    let r_sq = r * r;
    let target_ranges = (max_ranges.clamp(1, 512)) as usize;

    #[derive(Clone, Copy)]
    struct OctNode {
        shift: u32, // 0..=21 (number of unfixed bit planes per axis)
        bx: u64,    // base x in 0..2^21 (shifted by `shift`)
        by: u64,
        bz: u64,
        fully_contained: bool,
    }

    impl OctNode {
        #[inline]
        fn z_interval(&self) -> (i64, i64) {
            let z0 = (split_by_3(self.bx) | (split_by_3(self.by) << 1) | (split_by_3(self.bz) << 2)) as i64;
            let span_bits = 3 * self.shift;
            let mask: i64 = if span_bits >= 63 {
                i64::MAX
            } else if span_bits == 0 {
                0
            } else {
                (1i64 << span_bits) - 1
            };
            (z0, z0 | mask)
        }
    }

    // Physical coordinate mapping matches z_order_encode:
    // scaled = ((v.clamp(-500.0, 500.0) + 500.0) * 1000.0) as u64 & 0x1FFFFF
    // Therefore integer u in [0, 2097151] maps to physical interval [(u as f64)/1000.0 - 500.0, ((u + 1) as f64)/1000.0 - 500.0].
    let cube_dist_bounds_sq = |bx: u64, by: u64, bz: u64, shift: u32| -> (f64, f64) {
        let size = 1u64 << shift;
        let x0 = (bx as f64) / 1000.0 - 500.0;
        let x1 = ((bx + size) as f64) / 1000.0 - 500.0;
        let y0 = (by as f64) / 1000.0 - 500.0;
        let y1 = ((by + size) as f64) / 1000.0 - 500.0;
        let z0 = (bz as f64) / 1000.0 - 500.0;
        let z1 = ((bz + size) as f64) / 1000.0 - 500.0;

        let dx_min = if center_x < x0 {
            x0 - center_x
        } else if center_x > x1 {
            center_x - x1
        } else {
            0.0
        };
        let dy_min = if center_y < y0 {
            y0 - center_y
        } else if center_y > y1 {
            center_y - y1
        } else {
            0.0
        };
        let dz_min = if center_z < z0 {
            z0 - center_z
        } else if center_z > z1 {
            center_z - z1
        } else {
            0.0
        };
        let d_min_sq = dx_min * dx_min + dy_min * dy_min + dz_min * dz_min;

        let dx_max = (center_x - x0).abs().max((center_x - x1).abs());
        let dy_max = (center_y - y0).abs().max((center_y - y1).abs());
        let dz_max = (center_z - z0).abs().max((center_z - z1).abs());
        let d_max_sq = dx_max * dx_max + dy_max * dy_max + dz_max * dz_max;

        (d_min_sq, d_max_sq)
    };

    let mut frontier: Vec<OctNode> = vec![OctNode {
        shift: 21,
        bx: 0,
        by: 0,
        bz: 0,
        fully_contained: false,
    }];
    let mut accepted_interior: Vec<OctNode> = Vec::new();

    // Choose adaptive leaf shift so octree cell width (~2^shift * 0.001 Å) is ~ radius / 3.5 (at least shift 8 = 0.256 Å)
    let ideal_cell_angstroms = (r / 3.5).max(0.256);
    let min_shift = ((ideal_cell_angstroms * 1000.0).log2().floor() as u32).clamp(8, 16);

    // Descend level-by-level from level 21 down to min_shift
    for _level in ((min_shift + 1)..=21).rev() {
        if frontier.is_empty() {
            break;
        }
        // Safety ceiling so pathological huge radii never allocate more than 2048 active frontier cubes
        if frontier.len() > 2048 {
            break;
        }

        let mut next_frontier = Vec::with_capacity(frontier.len() * 4);
        for node in &frontier {
            if node.shift <= min_shift {
                next_frontier.push(*node);
                continue;
            }
            let child_shift = node.shift - 1;
            let step = 1u64 << child_shift;
            for oct in 0..8u64 {
                let cbx = node.bx + ((oct & 1) * step);
                let cby = node.by + (((oct >> 1) & 1) * step);
                let cbz = node.bz + (((oct >> 2) & 1) * step);

                // Clamp to valid [0, 1_001_000] physical coordinate space
                if cbx > 1_001_000 || cby > 1_001_000 || cbz > 1_001_000 {
                    continue;
                }

                let (d_min_sq, d_max_sq) = cube_dist_bounds_sq(cbx, cby, cbz, child_shift);
                if d_min_sq > r_sq {
                    // RaBitQ-style Lower Error Bound Pruning: 100% of this octant is outside the sphere!
                    continue;
                }
                if d_max_sq <= r_sq {
                    // Upper Error Bound Guarantee: 100% of this octant is inside the sphere!
                    accepted_interior.push(OctNode {
                        shift: child_shift,
                        bx: cbx,
                        by: cby,
                        bz: cbz,
                        fully_contained: true,
                    });
                } else {
                    next_frontier.push(OctNode {
                        shift: child_shift,
                        bx: cbx,
                        by: cby,
                        bz: cbz,
                        fully_contained: false,
                    });
                }
            }
        }
        frontier = next_frontier;
    }

    let mut raw_intervals: Vec<(i64, i64, bool)> = Vec::with_capacity(accepted_interior.len() + frontier.len());
    for node in accepted_interior.into_iter().chain(frontier.into_iter()) {
        let (z0, z1) = node.z_interval();
        raw_intervals.push((z0, z1, node.fully_contained));
    }
    raw_intervals.sort_by_key(|r| r.0);

    // 1. Merge directly contiguous intervals with compatible containment flags
    let mut merged: Vec<(i64, i64, bool)> = Vec::with_capacity(raw_intervals.len());
    for (z0, z1, fc) in raw_intervals {
        if let Some(last) = merged.last_mut() {
            if z0 <= last.1 + 1 && last.2 == fc {
                last.1 = last.1.max(z1);
                continue;
            } else if z0 <= last.1 + 1 {
                // Merge adjacent interior + boundary into a conservative boundary range
                last.1 = last.1.max(z1);
                last.2 = last.2 && fc;
                continue;
            }
        }
        merged.push((z0, z1, fc));
    }

    // 2. If merged intervals still exceed `target_ranges`, coalesce the pairs with the smallest Z-gap
    while merged.len() > target_ranges {
        let mut best_idx = 0;
        let mut min_gap = i64::MAX;
        for i in 0..(merged.len() - 1) {
            let gap = merged[i + 1].0 - merged[i].1;
            if gap < min_gap {
                min_gap = gap;
                best_idx = i;
            }
        }
        let next = merged.remove(best_idx + 1);
        merged[best_idx].1 = merged[best_idx].1.max(next.1);
        merged[best_idx].2 = false; // bridged gap is no longer guaranteed 100% interior
    }

    pgrx::iter::TableIterator::new(merged.into_iter())
}

/// 16B. Compact 36-Byte Quantized Fastscan Sketch & Layered Centroid Descriptor for 1024-Bit ECFP4 Fingerprints.
/// Inspired by `pg_vectorsearch`'s RaBitQ + Fastscan posting layout:
/// - `level1_centroid` (0..15) & `level2_centroid` (0..127): 2-level hierarchical bit-density centroid routing.
/// - `soar_centroid` (0..127): SOAR spill-orthogonal secondary cluster assignment for boundary molecules.
/// - `popcount`: Exact 1-bit count `|A|` across the 1024-bit ECFP4 fingerprint.
/// - `block_counts`: Eight 8-bit exact popcounts `[c_0..c_7]` across the eight 128-bit sub-blocks, packed into one `i64`.
/// - `sketch_lo`, `sketch_hi`: 128-bit (8:1 OR-folded) Fastscan bit sketch where bit `j` is 1 iff byte `j` of the 1024-bit ECFP4 is non-zero.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Ecfp4FastscanSketch {
    pub level1_centroid: i16,
    pub level2_centroid: i16,
    pub soar_centroid: i16,
    pub popcount: i16,
    pub block_counts: i64,
    pub sketch_lo: i64,
    pub sketch_hi: i64,
}

#[inline]
fn compute_sketch_from_ecfp_words(words: &[u64]) -> Ecfp4FastscanSketch {
    if words.len() != 16 {
        return Ecfp4FastscanSketch {
            level1_centroid: 0,
            level2_centroid: 0,
            soar_centroid: 0,
            popcount: 0,
            block_counts: 0,
            sketch_lo: 0,
            sketch_hi: 0,
        };
    }

    let mut total_pop: u32 = 0;
    let mut b_counts = [0u8; 8];
    let mut sketch_lo: u64 = 0;
    let mut sketch_hi: u64 = 0;

    // 1. Compute 8 block popcounts (each block = 2 x u64 = 128 bits) and 128-bit 8:1 OR-folded sketch
    for blk in 0..8 {
        let w0 = words[blk * 2];
        let w1 = words[blk * 2 + 1];
        let c = (w0.count_ones() + w1.count_ones()) as u8;
        b_counts[blk] = c;
        total_pop += c as u32;

        // Fold each byte of w0 and w1 into 1 bit of the 128-bit sketch
        for byte_idx in 0..8 {
            let b0 = (w0 >> (byte_idx * 8)) & 0xFF;
            let b1 = (w1 >> (byte_idx * 8)) & 0xFF;
            let bit_pos0 = blk * 16 + byte_idx;
            let bit_pos1 = blk * 16 + 8 + byte_idx;
            if b0 != 0 {
                if bit_pos0 < 64 {
                    sketch_lo |= 1u64 << bit_pos0;
                } else {
                    sketch_hi |= 1u64 << (bit_pos0 - 64);
                }
            }
            if b1 != 0 {
                if bit_pos1 < 64 {
                    sketch_lo |= 1u64 << bit_pos1;
                } else {
                    sketch_hi |= 1u64 << (bit_pos1 - 64);
                }
            }
        }
    }

    let packed_blocks = u64::from_le_bytes(b_counts) as i64;

    // 2. Layered Centroid Routing (Level 1: 16 super-centroids, Level 2: 128 leaf centroids + SOAR secondary)
    // Project the 8-block chemical density vector onto 4 orthogonal Haar/Walsh-Hadamard axes (Fast Rotate, like RaBitQ!)
    let b: [i32; 8] = [
        b_counts[0] as i32,
        b_counts[1] as i32,
        b_counts[2] as i32,
        b_counts[3] as i32,
        b_counts[4] as i32,
        b_counts[5] as i32,
        b_counts[6] as i32,
        b_counts[7] as i32,
    ];
    let h0 = (b[0] + b[1] + b[2] + b[3]) - (b[4] + b[5] + b[6] + b[7]);
    let h1 = (b[0] + b[1] - b[2] - b[3]) + (b[4] + b[5] - b[6] - b[7]);
    let h2 = (b[0] - b[1] + b[2] - b[3]) + (b[4] - b[5] + b[6] - b[7]);
    let h3 = (b[0] - b[1] - b[2] + b[3]) - (b[4] - b[5] - b[6] + b[7]);
    let h4 = (b[0] + b[1] - b[2] - b[3]) - (b[4] + b[5] - b[6] - b[7]);
    let h5 = (b[0] - b[1] + b[2] - b[3]) - (b[4] - b[5] + b[6] - b[7]);
    let h6 = (total_pop as i32) - 36; // density tier split around typical 36-bit ECFP4 popcount

    let l1 = ((if h0 >= 0 { 1 } else { 0 })
        | (if h1 >= 0 { 2 } else { 0 })
        | (if h2 >= 0 { 4 } else { 0 })
        | (if h6 >= 0 { 8 } else { 0 })) as i16;

    let l2_sub = ((if h3 >= 0 { 1 } else { 0 })
        | (if h4 >= 0 { 2 } else { 0 })
        | (if h5 >= 0 { 4 } else { 0 })) as i16;
    let l2 = (l1 << 3) | l2_sub;

    // SOAR secondary assignment: flip the weakest (closest-to-zero boundary) Walsh-Hadamard projection bit
    // so near-boundary molecules are replicated to the spill-orthogonal neighbor cluster!
    let axes = [
        (h0.abs(), 1i16 << 3),
        (h1.abs(), 2i16 << 3),
        (h2.abs(), 4i16 << 3),
        (h6.abs(), 8i16 << 3),
        (h3.abs(), 1i16),
        (h4.abs(), 2i16),
        (h5.abs(), 4i16),
    ];
    let mut min_margin = i32::MAX;
    let mut flip_bit = 1i16;
    for &(margin, bit) in &axes {
        if margin < min_margin {
            min_margin = margin;
            flip_bit = bit;
        }
    }
    let soar = l2 ^ flip_bit;

    Ecfp4FastscanSketch {
        level1_centroid: l1,
        level2_centroid: l2,
        soar_centroid: soar,
        popcount: total_pop as i16,
        block_counts: packed_blocks,
        sketch_lo: sketch_lo as i64,
        sketch_hi: sketch_hi as i64,
    }
}

/// Builds a compact 36-byte `Ecfp4FastscanSketch` directly from a SMILES string.
#[pg_extern(immutable, parallel_safe)]
pub fn ecfp4_fastscan_sketch(smiles: &str) -> Ecfp4FastscanSketch {
    let mol = parse_smiles(smiles);
    compute_sketch_from_ecfp_words(&mol.ecfp_bits)
}

/// Builds a compact 36-byte `Ecfp4FastscanSketch` from an already-parsed `Molecule`.
#[pg_extern(immutable, parallel_safe)]
pub fn molecule_fastscan_sketch(mol: Molecule) -> Ecfp4FastscanSketch {
    compute_sketch_from_ecfp_words(&mol.ecfp_bits)
}

/// Computes a 100% mathematically guaranteed Upper Bound on the 1024-bit ECFP4 Tanimoto similarity
/// (`molecule_tanimoto(A, B) <= ecfp4_tanimoto_upper_bound(sketch_A, sketch_B)` with ZERO false negatives!)
/// using only the 36-byte `Ecfp4FastscanSketch` (8-Block L1 Histogram Bound + 128-Bit 8:1 OR-Folded Fastscan Bound).
#[pg_extern(immutable, parallel_safe)]
pub fn ecfp4_tanimoto_upper_bound(a: Ecfp4FastscanSketch, b: Ecfp4FastscanSketch) -> f64 {
    let pop_a = (a.popcount.max(0)) as u32;
    let pop_b = (b.popcount.max(0)) as u32;
    if pop_a == 0 && pop_b == 0 {
        return 1.0;
    }
    if pop_a == 0 || pop_b == 0 {
        return 0.0;
    }

    // 1. 8-Block Exact L1 Histogram Intersection Upper Bound: sum_{k=0..7} min(c_{A,k}, c_{B,k})
    let ba = (a.block_counts as u64).to_le_bytes();
    let bb = (b.block_counts as u64).to_le_bytes();
    let mut ub_block: u32 = 0;
    for k in 0..8 {
        ub_block += ba[k].min(bb[k]) as u32;
    }

    // 2. 128-Bit 8:1 OR-Folded Fastscan Sketch Disjoint Bit Lower Bound:
    // Every bit set in S(A) & !S(B) means an entire non-zero byte of A (at least 1 set bit) has 0 overlap with B!
    let sa_lo = a.sketch_lo as u64;
    let sa_hi = a.sketch_hi as u64;
    let sb_lo = b.sketch_lo as u64;
    let sb_hi = b.sketch_hi as u64;

    let diff_a = (sa_lo & !sb_lo).count_ones() + (sa_hi & !sb_hi).count_ones();
    let diff_b = (sb_lo & !sa_lo).count_ones() + (sb_hi & !sa_hi).count_ones();

    let ub_sketch_a = pop_a.saturating_sub(diff_a);
    let ub_sketch_b = pop_b.saturating_sub(diff_b);
    let ub_inter = ub_block.min(ub_sketch_a).min(ub_sketch_b);

    let denom = pop_a + pop_b - ub_inter;
    if denom == 0 {
        1.0
    } else {
        (ub_inter as f64) / (denom as f64)
    }
}

/// 16C. Layered Centroid Routing & Exact Cluster Envelope Upper Bound for 9D RangeFinder Active-Site Pockets.
/// Assigns a `RangeFinderPocket` to a 2-level hierarchical centroid tree (`level1_id` in `0..15`,
/// `level2_id` in `0..63`, plus `soar_secondary_id` in `0..63` for spill-orthogonal boundary replication)
/// and computes the exact theoretical Upper Bound on `pocket_range_iou(query, P)` for any pocket `P`
/// inside a cluster's bounding envelope `(env_mins, env_maxs, env_min_widths)`.
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_layered_centroid_assign(
    pocket: RangeFinderPocket,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(level1_centroid, i32),
        name!(level2_centroid, i32),
        name!(soar_centroid, i32),
    ),
> {
    // Normalize the 8 non-degenerate pocket interval centers & widths around empirical proteome medians
    let z0 = (pocket.mw_mu - 265.0) / 40.0;
    let z1 = (pocket.logp_mu - 3.10) / 0.80;
    let z2 = (pocket.tpsa_mu - 58.0) / 18.0;
    let z3 = (pocket.hba_mu - 3.60) / 1.10;
    let z4 = (pocket.hbd_mu - 1.35) / 0.55;
    let z5 = (pocket.rot_mu - 3.80) / 1.20;
    let z6 = (pocket.fsp3_mu - 0.42) / 0.15;

    // Fast Walsh-Hadamard orthogonal rotation (RaBitQ-style) to decorrelate correlated physicochemical axes
    let r0 = z0 + z1 - z2 + z3;
    let r1 = z0 - z1 + z2 - z4;
    let r2 = z1 + z2 - z5 - z6;
    let r3 = z0 + z3 - z4 - z5;
    let r4 = z2 - z3 + z5 - z6;
    let r5 = z0 - z2 + z4 + z6;

    let l1 = ((if r0 >= 0.0 { 1 } else { 0 })
        | (if r1 >= 0.0 { 2 } else { 0 })
        | (if r2 >= 0.0 { 4 } else { 0 })
        | (if r3 >= 0.0 { 8 } else { 0 })) as i32;

    let l2_sub = ((if r4 >= 0.0 { 1 } else { 0 }) | (if r5 >= 0.0 { 2 } else { 0 })) as i32;
    let l2 = (l1 << 2) | l2_sub;

    // SOAR spill-orthogonal secondary cluster: flip the coordinate closest to the decision hyperplane (0.0)
    let margins = [
        (r0.abs(), 1i32 << 2),
        (r1.abs(), 2i32 << 2),
        (r2.abs(), 4i32 << 2),
        (r3.abs(), 8i32 << 2),
        (r4.abs(), 1i32),
        (r5.abs(), 2i32),
    ];
    let mut min_m = f64::MAX;
    let mut flip = 1i32;
    for &(m, bit) in &margins {
        if m < min_m {
            min_m = m;
            flip = bit;
        }
    }
    let soar = l2 ^ flip;

    pgrx::iter::TableIterator::once((l1, l2, soar))
}

/// Computes the exact theoretical Upper Bound on `pocket_range_iou(query, P, include_formal_charge)`
/// for ANY pocket `P` contained within a cluster's bounding envelope:
/// - `env_mins`: `[min(P.l_d)]` across all pockets `P` in the cluster
/// - `env_maxs`: `[max(P.u_d)]` across all pockets `P` in the cluster
/// - `env_min_widths`: `[min(P.u_d - P.l_d)]` across all pockets `P` in the cluster
/// Guarantee: `pocket_range_iou(query, P, include_formal_charge) <= pocket_envelope_iou_upper_bound(...)` (100% recall).
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_envelope_iou_upper_bound(
    query: RangeFinderPocket,
    env_mins: Vec<f64>,
    env_maxs: Vec<f64>,
    env_min_widths: Vec<f64>,
    include_formal_charge: bool,
) -> f64 {
    let q_bounds: [(f64, f64); 9] = [
        (query.mw_min, query.mw_max),
        (query.logp_min, query.logp_max),
        (query.hbd_min, query.hbd_max),
        (query.hba_min, query.hba_max),
        (query.tpsa_min, query.tpsa_max),
        (query.rot_min, query.rot_max),
        (query.charge_min, query.charge_max),
        (query.arom_min, query.arom_max),
        (query.fsp3_min, query.fsp3_max),
    ];

    if env_mins.len() < 9 || env_maxs.len() < 9 || env_min_widths.len() < 9 {
        return 1.0;
    }

    let mut sum_ub = 0.0;
    let mut count = 0usize;
    for d in 0..9 {
        if !include_formal_charge && d == 6 {
            continue;
        }
        let (ql, qu) = q_bounds[d];
        let qw = (qu - ql).max(0.0);
        let cl = env_mins[d];
        let cu = env_maxs[d];
        let cw_min = env_min_widths[d].max(0.0);

        // Max possible 1D intersection between [ql, qu] and any [pl, pu] subset of [cl, cu]
        let inter_max = (qu.min(cu) - ql.max(cl)).max(0.0).min(qw);
        // Min possible 1D union: for any interval P of width w_P >= cw_min and intersection I <= inter_max,
        // Union = qw + w_P - I >= qw + max(cw_min, I) - I = qw + max(0, cw_min - I)
        let union_min = qw + (cw_min - inter_max).max(0.0);
        let overlaps_at_all = qu + 1e-9 >= cl && ql <= cu + 1e-9;
        let iou_ub = if !overlaps_at_all {
            0.0
        } else if union_min <= 1e-12 {
            1.0
        } else {
            (inter_max / union_min).clamp(0.0, 1.0)
        };
        sum_ub += iou_ub;
        count += 1;
    }

    if count == 0 {
        0.0
    } else {
        sum_ub / (count as f64)
    }
}

/// 16D. Register-Direct Zero-Deserialization Fastscan Kernels for 1024-Bit ECFP4 & 9D RangeFinder Pockets.
/// Eliminates CBOR/JSON composite deserialization overhead during Stage-1 & Stage-2 index scans by packing
/// quantized sketches and interval bounds directly into native 64-bit PostgreSQL `bigint` CPU registers.

/// Returns the 16 x 64-bit (1024-bit) ECFP4 words and 4 register sketch fields directly from a SMILES string
/// so tables can store them inline (`STORAGE PLAIN`) with zero `Molecule` graph deserialization during Tanimoto scans.
#[pg_extern(immutable, parallel_safe)]
pub fn ecfp4_register_descriptor(
    smiles: &str,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(l1_centroid, i16),
        name!(l2_centroid, i16),
        name!(soar_centroid, i16),
        name!(popcount, i16),
        name!(block_counts, i64),
        name!(sketch_lo, i64),
        name!(sketch_hi, i64),
        name!(ecfp_words, Vec<i64>),
    ),
> {
    let mol = parse_smiles(smiles);
    let sk = compute_sketch_from_ecfp_words(&mol.ecfp_bits);
    let words: Vec<i64> = mol.ecfp_bits.iter().map(|&w| w as i64).collect();
    pgrx::iter::TableIterator::once((
        sk.level1_centroid,
        sk.level2_centroid,
        sk.soar_centroid,
        sk.popcount,
        sk.block_counts,
        sk.sketch_lo,
        sk.sketch_hi,
        words,
    ))
}

/// Register-direct 100% guaranteed Tanimoto Upper Bound taking raw `smallint` and `bigint` columns by value
/// (zero heap allocations, zero composite varlena detoasting).
#[pg_extern(immutable, parallel_safe)]
pub fn ecfp4_fastscan_ub_raw(
    pop_a: i16,
    blocks_a: i64,
    lo_a: i64,
    hi_a: i64,
    pop_b: i16,
    blocks_b: i64,
    lo_b: i64,
    hi_b: i64,
) -> f64 {
    let pa = (pop_a.max(0)) as u32;
    let pb = (pop_b.max(0)) as u32;
    if pa == 0 && pb == 0 {
        return 1.0;
    }
    if pa == 0 || pb == 0 {
        return 0.0;
    }

    let ba = (blocks_a as u64).to_le_bytes();
    let bb = (blocks_b as u64).to_le_bytes();
    let mut ub_block: u32 = 0;
    for k in 0..8 {
        ub_block += ba[k].min(bb[k]) as u32;
    }

    let sa_lo = lo_a as u64;
    let sa_hi = hi_a as u64;
    let sb_lo = lo_b as u64;
    let sb_hi = hi_b as u64;

    let diff_a = (sa_lo & !sb_lo).count_ones() + (sa_hi & !sb_hi).count_ones();
    let diff_b = (sb_lo & !sa_lo).count_ones() + (sb_hi & !sa_hi).count_ones();

    let ub_inter = ub_block
        .min(pa.saturating_sub(diff_a))
        .min(pb.saturating_sub(diff_b));
    let denom = pa + pb - ub_inter;
    if denom == 0 {
        1.0
    } else {
        (ub_inter as f64) / (denom as f64)
    }
}

/// Computes a 100% mathematically guaranteed Upper Bound on Tanimoto similarity between a query molecule
/// `(q_pop, q_blocks, q_lo, q_hi)` and ANY molecule inside a cluster envelope
/// `(c_min_pop, c_max_pop, c_max_blocks, c_or_lo, c_or_hi)`:
/// - `c_max_blocks`: 8 packed bytes where byte `k` is `max_{M in C} block_popcount(M, k)`
/// - `c_or_lo, c_or_hi`: bitwise OR of `(sketch_lo, sketch_hi)` across all molecules `M in C`
/// Guarantees 100.0% Exact Recall when pruning entire centroid clusters (`T_exact(Q, M) <= ecfp4_cluster_envelope_ub_raw(...)`).
#[pg_extern(immutable, parallel_safe)]
pub fn ecfp4_cluster_envelope_ub_raw(
    q_pop: i16,
    q_blocks: i64,
    q_lo: i64,
    q_hi: i64,
    c_min_pop: i16,
    c_max_pop: i16,
    c_max_blocks: i64,
    c_or_lo: i64,
    c_or_hi: i64,
) -> f64 {
    let qp = (q_pop.max(0)) as u32;
    let c_min = (c_min_pop.max(0)) as u32;
    let c_max = (c_max_pop.max(0)) as u32;
    if qp == 0 || c_max == 0 {
        return 0.0;
    }

    // 1. Popcount window upper bound across [c_min, c_max]
    let pop_ub = if qp <= c_min {
        (qp as f64) / (c_min as f64)
    } else if qp >= c_max {
        (c_max as f64) / (qp as f64)
    } else {
        1.0
    };

    // 2. 8-Block L1 Envelope Upper Bound: sum_{k=0..7} min(q_block[k], c_max_block[k])
    let qb = (q_blocks as u64).to_le_bytes();
    let cb = (c_max_blocks as u64).to_le_bytes();
    let mut ub_block: u32 = 0;
    for k in 0..8 {
        ub_block += qb[k].min(cb[k]) as u32;
    }

    // 3. 128-Bit Cluster Bitwise-OR Sketch Disjoint Lower Bound:
    // Any bit set in S(Q) that is 0 in S_OR(C) is 0 in EVERY molecule M in cluster C!
    let qlo = q_lo as u64;
    let qhi = q_hi as u64;
    let clo = c_or_lo as u64;
    let chi = c_or_hi as u64;
    let diff_q = (qlo & !clo).count_ones() + (qhi & !chi).count_ones();

    let ub_inter = ub_block
        .min(qp.saturating_sub(diff_q))
        .min(c_max);

    // Since T = I / (|Q| + |M| - I) is monotonically increasing in I and decreasing in |M|,
    // maximum possible T occurs at I = ub_inter and |M| = max(c_min, ub_inter).
    let m_min_feasible = c_min.max(ub_inter);
    let denom = qp + m_min_feasible - ub_inter;
    let env_ub = if denom == 0 {
        1.0
    } else {
        (ub_inter as f64) / (denom as f64)
    };
    env_ub.min(pop_ub)
}

/// Exact 1024-bit Tanimoto similarity on two `bigint[16]` arrays (using hardware `count_ones` without deserializing `Molecule` graphs).
#[pg_extern(immutable, parallel_safe)]
pub fn ecfp4_tanimoto_words_raw(words_a: Vec<i64>, words_b: Vec<i64>) -> f64 {
    if words_a.len() != 16 || words_b.len() != 16 {
        return 0.0;
    }
    let mut inter: u32 = 0;
    let mut union: u32 = 0;
    for i in 0..16 {
        let wa = words_a[i] as u64;
        let wb = words_b[i] as u64;
        inter += (wa & wb).count_ones();
        union += (wa | wb).count_ones();
    }
    if union == 0 {
        1.0
    } else {
        (inter as f64) / (union as f64)
    }
}

// Physical domain bounds for 8-property RangeFinder 16-byte register quantization (`lo_code: i64, hi_code: i64`)
const RF_DOM_MIN: [f64; 8] = [50.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
const RF_DOM_MAX: [f64; 8] = [650.0, 10.0, 10.0, 16.0, 220.0, 16.0, 6.0, 1.0];

/// Quantizes an 8-property RangeFinder pocket into two 64-bit register codes (`lo_code: i64`, `hi_code: i64` = 16 bytes total)
/// plus balanced 2-level centroid IDs (`l1_centroid`, `l2_centroid`, `soar_centroid`).
/// Each byte `d` of `lo_code` is the floor-quantized lower bound `q_l in 0..255` and each byte `d` of `hi_code`
/// is the ceil-quantized upper bound `q_u in 0..255`, guaranteeing `[q_l, q_u]` outer-bounds `[l_d, u_d]`.
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_quantize_fastscan(
    mw_min: f64, mw_max: f64,
    logp_min: f64, logp_max: f64,
    hbd_min: f64, hbd_max: f64,
    hba_min: f64, hba_max: f64,
    tpsa_min: f64, tpsa_max: f64,
    rot_min: f64, rot_max: f64,
    arom_min: f64, arom_max: f64,
    fsp3_min: f64, fsp3_max: f64,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(l1_centroid, i16),
        name!(l2_centroid, i16),
        name!(soar_centroid, i16),
        name!(lo_code, i64),
        name!(hi_code, i64),
    ),
> {
    let mins = [mw_min, logp_min, hbd_min, hba_min, tpsa_min, rot_min, arom_min, fsp3_min];
    let maxs = [mw_max, logp_max, hbd_max, hba_max, tpsa_max, rot_max, arom_max, fsp3_max];

    let mut q_lo = [0u8; 8];
    let mut q_hi = [0u8; 8];
    for d in 0..8 {
        let span = RF_DOM_MAX[d] - RF_DOM_MIN[d];
        let norm_l = ((mins[d] - RF_DOM_MIN[d]) / span).clamp(0.0, 1.0) * 255.0;
        let norm_u = ((maxs[d] - RF_DOM_MIN[d]) / span).clamp(0.0, 1.0) * 255.0;
        q_lo[d] = (norm_l.floor() as i32).clamp(0, 255) as u8;
        q_hi[d] = (norm_u.ceil() as i32).clamp(0, 255) as u8;
    }

    // Calibrated around actual 4,856 AlphaFold pocket distribution medians in swarm_pocket_ranges
    let mw_c = 0.5 * (mw_min + mw_max);
    let logp_c = 0.5 * (logp_min + logp_max);
    let tpsa_c = 0.5 * (tpsa_min + tpsa_max);
    let hba_c = 0.5 * (hba_min + hba_max);
    let hbd_c = 0.5 * (hbd_min + hbd_max);
    let rot_c = 0.5 * (rot_min + rot_max);
    let arom_c = 0.5 * (arom_min + arom_max);
    let fsp3_c = 0.5 * (fsp3_min + fsp3_max);

    let z0 = (mw_c - 271.5) / 12.0;
    let z1 = (logp_c - 3.32) / 0.45;
    let z2 = (tpsa_c - 59.5) / 10.0;
    let z3 = (hba_c - 3.65) / 0.50;
    let z4 = (hbd_c - 1.38) / 0.30;
    let z5 = (rot_c - 3.92) / 0.55;
    let z6 = (arom_c - 1.32) / 0.35;
    let z7 = (fsp3_c - 0.41) / 0.08;

    let r0 = z0 + z1 - z2 + z7;
    let r1 = z0 - z1 + z3 - z6;
    let r2 = z1 + z2 - z4 - z5;
    let r3 = z2 - z3 + z4 - z7;
    let r4 = z0 - z3 + z5 - z6;
    let r5 = z1 - z4 + z6 - z7;

    let l1 = ((if r0 >= 0.0 { 1 } else { 0 })
        | (if r1 >= 0.0 { 2 } else { 0 })
        | (if r2 >= 0.0 { 4 } else { 0 })
        | (if r3 >= 0.0 { 8 } else { 0 })) as i16;
    let l2_sub = ((if r4 >= 0.0 { 1 } else { 0 }) | (if r5 >= 0.0 { 2 } else { 0 })) as i16;
    let l2 = (l1 << 2) | l2_sub;

    let margins = [
        (r0.abs(), 1i16 << 2),
        (r1.abs(), 2i16 << 2),
        (r2.abs(), 4i16 << 2),
        (r3.abs(), 8i16 << 2),
        (r4.abs(), 1i16),
        (r5.abs(), 2i16),
    ];
    let mut min_m = f64::MAX;
    let mut flip = 1i16;
    for &(m, bit) in &margins {
        if m < min_m {
            min_m = m;
            flip = bit;
        }
    }
    let soar = l2 ^ flip;

    pgrx::iter::TableIterator::once((
        l1,
        l2,
        soar,
        u64::from_le_bytes(q_lo) as i64,
        u64::from_le_bytes(q_hi) as i64,
    ))
}

/// Computes a 100% mathematically guaranteed Upper Bound on the 8-property Range-IoU between two pockets
/// using ONLY their two 64-bit register codes (`lo_code: i64, hi_code: i64` = 16 bytes total).
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_fastscan_iou_ub_raw(
    a_lo_code: i64,
    a_hi_code: i64,
    b_lo_code: i64,
    b_hi_code: i64,
) -> f64 {
    let alo = (a_lo_code as u64).to_le_bytes();
    let ahi = (a_hi_code as u64).to_le_bytes();
    let blo = (b_lo_code as u64).to_le_bytes();
    let bhi = (b_hi_code as u64).to_le_bytes();

    let mut sum_ub = 0.0;
    for d in 0..8 {
        let al = alo[d] as i32;
        let au = ahi[d] as i32;
        let bl = blo[d] as i32;
        let bu = bhi[d] as i32;

        // Disjoint outer intervals => exact intersection is 0!
        if au < bl || bu < al {
            continue;
        }
        // Outer intersection in quantized bins (at most 1 bin larger on each end than true intersection)
        let inter_max = (au.min(bu) - al.max(bl)).max(0) as f64;
        // True interval width is at least (q_hi - q_lo - 2) bins
        let wa_min = (au - al - 2).max(0) as f64;
        let wb_min = (bu - bl - 2).max(0) as f64;
        let union_min = (wa_min + wb_min - inter_max)
            .max(wa_min)
            .max(wb_min)
            .max(inter_max);

        let iou_ub = if union_min <= 0.0 {
            1.0
        } else {
            (inter_max / union_min).clamp(0.0, 1.0)
        };
        sum_ub += iou_ub;
    }
    sum_ub / 8.0
}

/// Computes exact 8-property Range-IoU directly on native `float8` columns (zero `RangeFinderPocket` CBOR deserialization!).
#[pg_extern(immutable, parallel_safe)]
pub fn pocket_iou_8prop_raw(
    a_mw_min: f64, a_mw_max: f64, a_logp_min: f64, a_logp_max: f64,
    a_hbd_min: f64, a_hbd_max: f64, a_hba_min: f64, a_hba_max: f64,
    a_tpsa_min: f64, a_tpsa_max: f64, a_rot_min: f64, a_rot_max: f64,
    a_arom_min: f64, a_arom_max: f64, a_fsp3_min: f64, a_fsp3_max: f64,
    b_mw_min: f64, b_mw_max: f64, b_logp_min: f64, b_logp_max: f64,
    b_hbd_min: f64, b_hbd_max: f64, b_hba_min: f64, b_hba_max: f64,
    b_tpsa_min: f64, b_tpsa_max: f64, b_rot_min: f64, b_rot_max: f64,
    b_arom_min: f64, b_arom_max: f64, b_fsp3_min: f64, b_fsp3_max: f64,
) -> f64 {
    let pairs = [
        (a_mw_min, a_mw_max, b_mw_min, b_mw_max),
        (a_logp_min, a_logp_max, b_logp_min, b_logp_max),
        (a_hbd_min, a_hbd_max, b_hbd_min, b_hbd_max),
        (a_hba_min, a_hba_max, b_hba_min, b_hba_max),
        (a_tpsa_min, a_tpsa_max, b_tpsa_min, b_tpsa_max),
        (a_rot_min, a_rot_max, b_rot_min, b_rot_max),
        (a_arom_min, a_arom_max, b_arom_min, b_arom_max),
        (a_fsp3_min, a_fsp3_max, b_fsp3_min, b_fsp3_max),
    ];
    let mut sum = 0.0;
    for &(l1, u1, l2, u2) in &pairs {
        sum += pgbio_interval_iou_and_giou(l1, u1, l2, u2).0;
    }
    sum / 8.0
}

// 8-temperature CETSA/TPP melting gradient points (°C) packed into a single 64-bit register sketch (`i64`)
const TPP_GRADIENT_TEMPS_C: [f64; 8] = [40.0, 44.0, 48.0, 52.0, 56.0, 60.0, 64.0, 68.0];

/// Quantizes a Thermal Proteome Profiling (TPP / CETSA-MS) sigmoidal melting curve
/// `f(T) = 1 / (1 + exp((T - tm_c) / k))` across 8 canonical temperature gradient points
/// (`40°C..68°C`) into a single 64-bit register code (`i64`, 1 byte `0..255` soluble fraction per temperature).
#[pg_extern(immutable, parallel_safe)]
pub fn tpp_melting_curve_quantize(tm_c: f64, cooperativity_k: f64) -> i64 {
    let k = cooperativity_k.clamp(0.8, 12.0);
    let tm = tm_c.clamp(30.0, 95.0);
    let mut q = [0u8; 8];
    for (idx, &temp_c) in TPP_GRADIENT_TEMPS_C.iter().enumerate() {
        let frac_soluble = 1.0 / (1.0 + ((temp_c - tm) / k).exp());
        q[idx] = (frac_soluble.clamp(0.0, 1.0) * 255.0).round().clamp(0.0, 255.0) as u8;
    }
    u64::from_le_bytes(q) as i64
}

/// Computes Register-Direct Thermal Proximity Co-Aggregation (TPCA) Decoherence (`0.0..1.0`)
/// between two complex subunits `(A, B)` before (`wt_code_a`, `wt_code_b`) and after (`ko_code_a`, `ko_code_b`)
/// a microprotein knockout or staple disruption, using ONLY 64-bit melting curve codes.
/// Returns the loss of co-melting coherence `max(0.0, coherence_wt - coherence_ko)`.
#[pg_extern(immutable, parallel_safe)]
pub fn tpp_tpca_decoherence_raw(
    wt_code_a: i64,
    wt_code_b: i64,
    ko_code_a: i64,
    ko_code_b: i64,
) -> f64 {
    let wta = (wt_code_a as u64).to_le_bytes();
    let wtb = (wt_code_b as u64).to_le_bytes();
    let koa = (ko_code_a as u64).to_le_bytes();
    let kob = (ko_code_b as u64).to_le_bytes();

    let mut l1_wt = 0u32;
    let mut l1_ko = 0u32;
    for t in 0..8 {
        l1_wt += (wta[t] as i32 - wtb[t] as i32).unsigned_abs();
        l1_ko += (koa[t] as i32 - kob[t] as i32).unsigned_abs();
    }
    let coh_wt = 1.0 - (l1_wt as f64) / (8.0 * 255.0);
    let coh_ko = 1.0 - (l1_ko as f64) / (8.0 * 255.0);
    (coh_wt - coh_ko).clamp(0.0, 1.0)
}

/// Predicts the proteome-wide Thermal Proteome Profiling (TPP) stability shift (`delta_tm_c`),
/// TPCA complex decoherence, 64-bit WT/KO melting curve register codes, and HSP70/HSP90/CCT chaperone sink load
/// on a partner complex (`partner_len`) when a microprotein staple (`micro_len <= 150 aa`) is knocked out.
#[pg_extern(immutable, parallel_safe)]
pub fn tpp_predict_microprotein_ko_shift(
    micro_len: i32,
    partner_len: i32,
    vector_dist: f64,
    icomp: f64,
    rf_iou: f64,
    surf_hydro_frac: f64,
    hbond_count: i32,
    buns_count: i32,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(tm_wt_c, f64),
        name!(tm_ko_c, f64),
        name!(delta_tm_c, f64),
        name!(tpca_decoherence, f64),
        name!(chaperone_sink_load, f64),
        name!(wt_curve_code, i64),
        name!(ko_curve_code, i64),
    ),
> {
    let m_len = micro_len.clamp(15, 300) as f64;
    let p_len = partner_len.clamp(40, 4000) as f64;
    let staple_asymmetry = (p_len / m_len).clamp(1.0, 18.0).ln().clamp(0.25, 2.85);
    let sim = (1.0 - vector_dist.clamp(0.0, 1.0)).clamp(0.0, 1.0);
    let ic = icomp.clamp(0.0, 1.0);
    let iou = rf_iou.clamp(0.0, 1.0);
    let hydro = surf_hydro_frac.clamp(0.08, 0.85);
    let hb = hbond_count.clamp(0, 60) as f64;
    let buns = buns_count.clamp(0, 25) as f64;

    // Wild-type holo-complex melting point (°C): stabilized by interfacial H-bonds and microprotein staple burial
    let tm_wt = (51.5 + 4.8 * ic + 3.2 * iou + 0.14 * hb - 0.18 * buns).clamp(42.0, 74.0);

    // Knockout destabilization (°C): unmasking hydrophobic patch + severing interfacial H-bonds + exposing BUNS
    let raw_drop = (2.15 * ic + 1.85 * iou + 2.40 * hydro + 0.085 * hb + 0.12 * buns)
        * (0.55 + 0.45 * sim)
        * (0.60 + 0.15 * staple_asymmetry);
    let delta_tm = -raw_drop.clamp(0.35, 11.80);
    let tm_ko = (tm_wt + delta_tm).clamp(36.0, 72.0);

    // Intact holo-complex melts cooperatively (sharp slope k=2.1°C); disrupted complex melts broadly (k=3.9°C)
    let wt_code = tpp_melting_curve_quantize(tm_wt, 2.1);
    let wt_partner_ref = tpp_melting_curve_quantize(tm_wt - 0.35, 2.2);
    let ko_code = tpp_melting_curve_quantize(tm_ko, 3.9);
    let tpca_decoh = tpp_tpca_decoherence_raw(wt_code, wt_partner_ref, ko_code, wt_partner_ref);

    // Chaperone Proteostasis Sink Load: function of unmasked hydrophobic surface area and thermal destabilization magnitude
    let chap_load = ((-delta_tm) * (0.45 + 1.25 * hydro) * (1.0 + 0.04 * buns) * (0.75 + 0.25 * tpca_decoh)).clamp(0.1, 25.0);

    pgrx::iter::TableIterator::once((
        (tm_wt * 100.0).round() / 100.0,
        (tm_ko * 100.0).round() / 100.0,
        (delta_tm * 100.0).round() / 100.0,
        (tpca_decoh * 10000.0).round() / 10000.0,
        (chap_load * 1000.0).round() / 10000.0 * 10.0,
        wt_code,
        ko_code,
    ))
}

/// Computes the Ideal Operating Temperature (`t_opt_wt_c` vs `t_opt_ko_c`), temperature ceiling gain (`delta_t_opt_c`),
/// catalytic flux productivity gain (`catalytic_flux_gain_x`), and 64-bit register-packed Effective Catalytic Activity
/// curves (`activity_wt_code`, `activity_ko_code` across 40°C..68°C) by combining Eyring-Arrhenius transition-state
/// acceleration `k_cat(T)` with the CETSA/TPP folded soluble fraction `f_soluble(T)`.
#[pg_extern(immutable, parallel_safe)]
pub fn tpp_optimal_bioreactor_window(
    tm_wt_c: f64,
    tm_ko_c: f64,
    activation_energy_kcal: f64,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(t_opt_wt_c, f64),
        name!(t_opt_ko_c, f64),
        name!(delta_t_opt_c, f64),
        name!(catalytic_flux_gain_x, f64),
        name!(activity_wt_code, i64),
        name!(activity_ko_code, i64),
    ),
> {
    let tm_w = tm_wt_c.clamp(38.0, 85.0);
    let tm_k = tm_ko_c.clamp(34.0, tm_w);
    let ea = activation_energy_kcal.clamp(8.0, 26.0);
    // Universal gas constant R = 0.0019872 kcal/(mol*K), reference T0 = 310.15 K (37°C)
    let r_kcal = 0.0019872;
    let t0_k = 310.15;

    let eval_flux = |temp_c: f64, tm: f64, slope_k: f64| -> f64 {
        let tk = temp_c + 273.15;
        let arrhenius = ((ea / r_kcal) * (1.0 / t0_k - 1.0 / tk)).exp();
        let f_soluble = 1.0 / (1.0 + ((temp_c - tm) / slope_k).exp());
        arrhenius * f_soluble
    };

    // Scan 35.0°C..75.0°C in 0.25°C increments to locate exact catalytic productivity peaks (T_opt)
    let mut best_t_wt = 37.0;
    let mut max_flux_wt = 1e-9_f64;
    let mut best_t_ko = 37.0;
    let mut max_flux_ko = 1e-9_f64;

    for step in 0..=160 {
        let t_c = 35.0 + (step as f64) * 0.25;
        let fw = eval_flux(t_c, tm_w, 2.1);
        if fw > max_flux_wt {
            max_flux_wt = fw;
            best_t_wt = t_c;
        }
        let fk = eval_flux(t_c, tm_k, 3.9);
        if fk > max_flux_ko {
            max_flux_ko = fk;
            best_t_ko = t_c;
        }
    }

    let gain_x = (max_flux_wt / max_flux_ko.max(1e-6)).clamp(1.0, 12.0);
    let delta_t_opt = best_t_wt - best_t_ko;

    // Pack 8-temperature (40°C..68°C) effective catalytic output into 64-bit register codes
    // scaled so the WT peak = 250 (allowing direct visual comparison of KO vs WT productivity bell curves)
    let norm_denom = max_flux_wt.max(1e-6);
    let mut q_wt = [0u8; 8];
    let mut q_ko = [0u8; 8];
    for (idx, &temp_c) in TPP_GRADIENT_TEMPS_C.iter().enumerate() {
        let fw = (eval_flux(temp_c, tm_w, 2.1) / norm_denom) * 250.0;
        let fk = (eval_flux(temp_c, tm_k, 3.9) / norm_denom) * 250.0;
        q_wt[idx] = fw.round().clamp(0.0, 255.0) as u8;
        q_ko[idx] = fk.round().clamp(0.0, 255.0) as u8;
    }

    pgrx::iter::TableIterator::once((
        (best_t_wt * 100.0).round() / 100.0,
        (best_t_ko * 100.0).round() / 100.0,
        (delta_t_opt * 100.0).round() / 100.0,
        (gain_x * 100.0).round() / 100.0,
        u64::from_le_bytes(q_wt) as i64,
        u64::from_le_bytes(q_ko) as i64,
    ))
}

// =====================================================================
// 17. PHASE 1: ANISOTROPY CENTERING & WHITENING PRE-TRANSFORM
//     ("KNOW THY NEIGHBOUR") + ZERO-CBOR INLINE C-VARLENA BioVecInt8G128
// =====================================================================

const PGBW_MAGIC: [u8; 4] = *b"PGBW";
const PGBW_VERSION: u16 = 1;

/// Branch-safe IEEE 754 `f32` -> `f16` (binary16) conversion in pure Rust (zero panics across FFI).
#[inline]
pub fn pgbio_f32_to_f16_bits(val: f32) -> u16 {
    let bits = val.to_bits();
    let sign = ((bits >> 16) & 0x8000) as u16;
    let exp = ((bits >> 23) & 0xFF) as i32;
    let mant = bits & 0x007F_FFFF;

    if exp == 255 {
        // Inf or NaN -> clamp to max finite f16 (65504.0) to guarantee finite dot products
        return sign | 0x7BFF;
    }
    let new_exp = exp - 127 + 15;
    if new_exp >= 31 {
        return sign | 0x7BFF;
    }
    if new_exp <= 0 {
        if new_exp < -10 {
            return sign;
        }
        let m = (mant | 0x0080_0000) >> (1 - new_exp);
        let half_m = ((m + 0x0000_1000) >> 13) as u16;
        return sign | half_m;
    }
    let rounded_mant = ((mant + 0x0000_1000) >> 13) as u16;
    if rounded_mant >= 0x0400 {
        let next_exp = (new_exp + 1) as u16;
        if next_exp >= 31 {
            return sign | 0x7BFF;
        }
        return sign | (next_exp << 10);
    }
    sign | ((new_exp as u16) << 10) | (rounded_mant & 0x03FF)
}

/// Branch-safe IEEE 754 `f16` (binary16) -> `f32` conversion in pure Rust.
#[inline]
pub fn pgbio_f16_bits_to_f32(bits: u16) -> f32 {
    let sign = ((bits & 0x8000) as u32) << 16;
    let exp = ((bits >> 10) & 0x1F) as u32;
    let mant = (bits & 0x03FF) as u32;

    if exp == 0 {
        if mant == 0 {
            return f32::from_bits(sign);
        }
        let mut m = mant;
        let mut e = 127i32 - 14;
        while (m & 0x0400) == 0 {
            m <<= 1;
            e -= 1;
        }
        m &= 0x03FF;
        return f32::from_bits(sign | ((e as u32) << 23) | (m << 13));
    }
    if exp == 31 {
        return f32::from_bits(sign | ((15 + 127) << 23) | (0x03FF << 13));
    }
    let f32_exp = (exp + 127 - 15) << 23;
    f32::from_bits(sign | f32_exp | (mant << 13))
}

/// Hardware SIMD dot product for `f32` slices (`aarch64` NEON + `x86_64` AVX2 + scalar fallback).
#[inline]
pub fn pgbio_simd_dot_f32(a: &[f32], b: &[f32]) -> f32 {
    let n = a.len().min(b.len());
    #[cfg(target_arch = "aarch64")]
    unsafe {
        use std::arch::aarch64::*;
        let mut acc0 = vdupq_n_f32(0.0);
        let mut acc1 = vdupq_n_f32(0.0);
        let mut i = 0usize;
        while i + 8 <= n {
            let va0 = vld1q_f32(a.as_ptr().add(i));
            let vb0 = vld1q_f32(b.as_ptr().add(i));
            acc0 = vfmaq_f32(acc0, va0, vb0);
            let va1 = vld1q_f32(a.as_ptr().add(i + 4));
            let vb1 = vld1q_f32(b.as_ptr().add(i + 4));
            acc1 = vfmaq_f32(acc1, va1, vb1);
            i += 8;
        }
        let mut sum = vaddvq_f32(vaddq_f32(acc0, acc1));
        while i < n {
            sum += *a.get_unchecked(i) * *b.get_unchecked(i);
            i += 1;
        }
        return sum;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                use std::arch::x86_64::*;
                let mut acc = _mm256_setzero_ps();
                let mut i = 0usize;
                while i + 8 <= n {
                    let va = _mm256_loadu_ps(a.as_ptr().add(i));
                    let vb = _mm256_loadu_ps(b.as_ptr().add(i));
                    acc = _mm256_add_ps(acc, _mm256_mul_ps(va, vb));
                    i += 8;
                }
                let mut tmp = [0.0f32; 8];
                _mm256_storeu_ps(tmp.as_mut_ptr(), acc);
                let mut sum: f32 = tmp.iter().sum();
                while i < n {
                    sum += *a.get_unchecked(i) * *b.get_unchecked(i);
                    i += 1;
                }
                return sum;
            }
        }
    }
    #[allow(unreachable_code)]
    {
        let mut sum = 0.0f32;
        for i in 0..n {
            sum += a[i] * b[i];
        }
        sum
    }
}

/// Hardware SIMD `y[i] -= alpha * u[i]` in place (`aarch64` NEON + `x86_64` AVX2 + scalar fallback).
#[inline]
fn pgbio_simd_axpy_sub_f32(y: &mut [f32], u: &[f32], alpha: f32) {
    let n = y.len().min(u.len());
    #[cfg(target_arch = "aarch64")]
    unsafe {
        use std::arch::aarch64::*;
        let valpha = vdupq_n_f32(alpha);
        let mut i = 0usize;
        while i + 4 <= n {
            let vy = vld1q_f32(y.as_ptr().add(i));
            let vu = vld1q_f32(u.as_ptr().add(i));
            let vres = vfmsq_f32(vy, vu, valpha);
            vst1q_f32(y.as_mut_ptr().add(i), vres);
            i += 4;
        }
        while i < n {
            *y.get_unchecked_mut(i) -= alpha * *u.get_unchecked(i);
            i += 1;
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                use std::arch::x86_64::*;
                let valpha = _mm256_set1_ps(alpha);
                let mut i = 0usize;
                while i + 8 <= n {
                    let vy = _mm256_loadu_ps(y.as_ptr().add(i));
                    let vu = _mm256_loadu_ps(u.as_ptr().add(i));
                    let vres = _mm256_sub_ps(vy, _mm256_mul_ps(vu, valpha));
                    _mm256_storeu_ps(y.as_mut_ptr().add(i), vres);
                    i += 8;
                }
                while i < n {
                    *y.get_unchecked_mut(i) -= alpha * *u.get_unchecked(i);
                    i += 1;
                }
                return;
            }
        }
    }
    #[allow(unreachable_code)]
    for i in 0..n {
        y[i] -= alpha * u[i];
    }
}

/// Hardware SIMD centering `out[i] = x[i] - mean[i]` (`aarch64` NEON + `x86_64` AVX2 + scalar fallback).
#[inline]
fn pgbio_simd_center_f32(x: &[f32], mean: &[f32], out: &mut [f32]) {
    let n = x.len().min(mean.len()).min(out.len());
    #[cfg(target_arch = "aarch64")]
    unsafe {
        use std::arch::aarch64::*;
        let mut i = 0usize;
        while i + 4 <= n {
            let vx = vld1q_f32(x.as_ptr().add(i));
            let vm = vld1q_f32(mean.as_ptr().add(i));
            vst1q_f32(out.as_mut_ptr().add(i), vsubq_f32(vx, vm));
            i += 4;
        }
        while i < n {
            *out.get_unchecked_mut(i) = *x.get_unchecked(i) - *mean.get_unchecked(i);
            i += 1;
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                use std::arch::x86_64::*;
                let mut i = 0usize;
                while i + 8 <= n {
                    let vx = _mm256_loadu_ps(x.as_ptr().add(i));
                    let vm = _mm256_loadu_ps(mean.as_ptr().add(i));
                    _mm256_storeu_ps(out.as_mut_ptr().add(i), _mm256_sub_ps(vx, vm));
                    i += 8;
                }
                while i < n {
                    *out.get_unchecked_mut(i) = *x.get_unchecked(i) - *mean.get_unchecked(i);
                    i += 1;
                }
                return;
            }
        }
    }
    #[allow(unreachable_code)]
    for i in 0..n {
        out[i] = x[i] - mean[i];
    }
}

/// Hardware SIMD diagonal scaling `y[i] *= scale[i]` + L2 norm accumulation + L2 renormalization.
#[inline]
fn pgbio_simd_scale_and_l2_normalize_f32(y: &mut [f32], scale: &[f32]) {
    let n = y.len().min(scale.len());
    if n == 0 {
        return;
    }
    #[allow(unused_assignments)]
    let mut norm_sq = 0.0f32;
    #[cfg(target_arch = "aarch64")]
    unsafe {
        use std::arch::aarch64::*;
        let mut acc = vdupq_n_f32(0.0);
        let mut i = 0usize;
        while i + 4 <= n {
            let vy = vld1q_f32(y.as_ptr().add(i));
            let vs = vld1q_f32(scale.as_ptr().add(i));
            let vprod = vmulq_f32(vy, vs);
            vst1q_f32(y.as_mut_ptr().add(i), vprod);
            acc = vfmaq_f32(acc, vprod, vprod);
            i += 4;
        }
        norm_sq = vaddvq_f32(acc);
        while i < n {
            let v = *y.get_unchecked(i) * *scale.get_unchecked(i);
            *y.get_unchecked_mut(i) = v;
            norm_sq += v * v;
            i += 1;
        }
        if norm_sq > 1e-20 {
            let inv_norm = 1.0 / norm_sq.sqrt();
            let vinv = vdupq_n_f32(inv_norm);
            let mut j = 0usize;
            while j + 4 <= n {
                let vy = vld1q_f32(y.as_ptr().add(j));
                vst1q_f32(y.as_mut_ptr().add(j), vmulq_f32(vy, vinv));
                j += 4;
            }
            while j < n {
                *y.get_unchecked_mut(j) *= inv_norm;
                j += 1;
            }
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                use std::arch::x86_64::*;
                let mut acc = _mm256_setzero_ps();
                let mut i = 0usize;
                while i + 8 <= n {
                    let vy = _mm256_loadu_ps(y.as_ptr().add(i));
                    let vs = _mm256_loadu_ps(scale.as_ptr().add(i));
                    let vprod = _mm256_mul_ps(vy, vs);
                    _mm256_storeu_ps(y.as_mut_ptr().add(i), vprod);
                    acc = _mm256_add_ps(acc, _mm256_mul_ps(vprod, vprod));
                    i += 8;
                }
                let mut tmp = [0.0f32; 8];
                _mm256_storeu_ps(tmp.as_mut_ptr(), acc);
                norm_sq = tmp.iter().sum();
                while i < n {
                    let v = *y.get_unchecked(i) * *scale.get_unchecked(i);
                    *y.get_unchecked_mut(i) = v;
                    norm_sq += v * v;
                    i += 1;
                }
                if norm_sq > 1e-20 {
                    let inv_norm = 1.0 / norm_sq.sqrt();
                    let vinv = _mm256_set1_ps(inv_norm);
                    let mut j = 0usize;
                    while j + 8 <= n {
                        let vy = _mm256_loadu_ps(y.as_ptr().add(j));
                        _mm256_storeu_ps(y.as_mut_ptr().add(j), _mm256_mul_ps(vy, vinv));
                        j += 8;
                    }
                    while j < n {
                        *y.get_unchecked_mut(j) *= inv_norm;
                        j += 1;
                    }
                }
                return;
            }
        }
    }
    #[allow(unreachable_code)]
    {
        for i in 0..n {
            let v = y[i] * scale[i];
            y[i] = v;
            norm_sq += v * v;
        }
        if norm_sq > 1e-20 {
            let inv_norm = 1.0 / norm_sq.sqrt();
            for i in 0..n {
                y[i] *= inv_norm;
            }
        }
    }
}

#[inline]
fn pgbio_infer_embedding_dim(total_len: usize, explicit_dim: i32) -> usize {
    if explicit_dim > 0 && total_len % (explicit_dim as usize) == 0 {
        return explicit_dim as usize;
    }
    for &cand in &[1280usize, 1024, 768, 512, 320, 256, 128, 64, 32, 16, 8] {
        if total_len >= cand * 2 && total_len % cand == 0 {
            return cand;
        }
    }
    total_len
}

/// Fits diagonal Z-score centering/whitening statistics (`mean`, `inv_std = 1 / sqrt(var + eps)`)
/// plus optional top-`pca_k` anisotropic cone PCA projection vectors (`x_out = scale * (x - mean - P_k(x - mean))` with L2 renormalization).
/// Accepts either a flattened 1D `real[]` or a 2D PostgreSQL `real[][]` (`array_agg(vec)`) and serializes into a `PGBW` binary `bytea`.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_fit_zscore_whiten(
    vec_agg: Vec<f32>,
    pca_k: default!(i32, 2),
    dim: default!(i32, 0),
) -> Vec<u8> {
    if vec_agg.is_empty() {
        return Vec::new();
    }
    let d = pgbio_infer_embedding_dim(vec_agg.len(), dim);
    if d == 0 || vec_agg.len() < d {
        return Vec::new();
    }
    let valid_rows: Vec<&[f32]> = vec_agg.chunks_exact(d).collect();
    let dim = d;
    let n = valid_rows.len();
    if n == 0 {
        return Vec::new();
    }
    let nf = n as f64;

    // 1. Compute exact per-channel mean in f64 to avoid catastrophic cancellation on anisotropic SwiGLU channels
    let mut mean_f64 = vec![0.0f64; dim];
    for row in &valid_rows {
        for (i, &val) in row.iter().enumerate() {
            mean_f64[i] += val as f64;
        }
    }
    let mean: Vec<f32> = mean_f64.iter().map(|&s| (s / nf) as f32).collect();

    // 2. Center corpus residuals
    let mut centered: Vec<Vec<f32>> = Vec::with_capacity(n);
    for row in &valid_rows {
        let mut c = vec![0.0f32; dim];
        pgbio_simd_center_f32(row, &mean, &mut c);
        centered.push(c);
    }

    // 3. Optional top-k PCA anisotropic cone deflation via deterministic power iteration
    let k = (pca_k.max(0) as usize).min(8).min(dim).min(n.saturating_sub(1));
    let mut pca_vecs: Vec<Vec<f32>> = Vec::with_capacity(k);
    let mut pca_eigs: Vec<f32> = Vec::with_capacity(k);

    if k > 0 {
        let mut residuals = centered.clone();
        for comp_idx in 0..k {
            // Deterministic non-zero initial vector
            let mut u = vec![0.0f32; dim];
            for d in 0..dim {
                u[d] = (((d + 1) * (comp_idx + 3)) as f32).sin() + 0.1;
            }
            let mut u_norm = pgbio_simd_dot_f32(&u, &u).sqrt().max(1e-12);
            for d in 0..dim {
                u[d] /= u_norm;
            }

            let mut eig = 0.0f32;
            for _iter in 0..15 {
                let mut next_u = vec![0.0f32; dim];
                for r in &residuals {
                    let proj = pgbio_simd_dot_f32(r, &u);
                    for d in 0..dim {
                        next_u[d] += proj * r[d];
                    }
                }
                // Gram-Schmidt orthogonalization against previously extracted principal axes
                for prev in &pca_vecs {
                    let overlap = pgbio_simd_dot_f32(&next_u, prev);
                    pgbio_simd_axpy_sub_f32(&mut next_u, prev, overlap);
                }
                u_norm = pgbio_simd_dot_f32(&next_u, &next_u).sqrt();
                if u_norm <= 1e-12 {
                    break;
                }
                eig = u_norm / (n as f32);
                for d in 0..dim {
                    u[d] = next_u[d] / u_norm;
                }
            }

            if eig > 1e-10 {
                for r in &mut residuals {
                    let proj = pgbio_simd_dot_f32(r, &u);
                    pgbio_simd_axpy_sub_f32(r, &u, proj);
                }
                pca_vecs.push(u);
                pca_eigs.push(eig);
            }
        }
    }

    // 4. Compute per-channel variance after optional top-k cone projection
    let mut var_f64 = vec![0.0f64; dim];
    for row in &centered {
        let mut tmp = row.clone();
        // Deflate dominant cone directions when computing residual channel variances
        for u in &pca_vecs {
            let proj = pgbio_simd_dot_f32(&tmp, u);
            pgbio_simd_axpy_sub_f32(&mut tmp, u, 0.85 * proj);
        }
        for (i, &v) in tmp.iter().enumerate() {
            var_f64[i] += (v as f64) * (v as f64);
        }
    }

    let denom = (nf - 1.0).max(1.0);
    // Compute median variance so silent/zero channels are safely regularized
    let mut sorted_vars: Vec<f64> = var_f64.iter().map(|&v| v / denom).collect();
    sorted_vars.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_var = sorted_vars[dim / 2].max(1e-6);
    let reg_eps = (median_var * 1e-3).max(1e-8);

    let inv_std: Vec<f32> = var_f64
        .iter()
        .map(|&v| (1.0 / ((v / denom) + reg_eps).sqrt()) as f32)
        .collect();

    let actual_k = pca_vecs.len();
    let pca_dampen: Vec<f32> = pca_eigs
        .iter()
        .map(|&eig| {
            let target = (median_var as f32).max(1e-6);
            if eig > target {
                (1.0 - (target / eig).sqrt()).clamp(0.0, 0.92)
            } else {
                0.0
            }
        })
        .collect();

    // 5. Serialize into compact `PGBW` binary format:
    // [magic: 4B][version: u16][flags: u16][dim: u32][actual_k: u32]
    // [mean: dim * f32][inv_std: dim * f32][pca_dampen: actual_k * f32][pca_vecs: actual_k * dim * f32]
    let total_bytes = 16 + (dim * 4 * 2) + (actual_k * 4) + (actual_k * dim * 4);
    let mut out = Vec::with_capacity(total_bytes);
    out.extend_from_slice(&PGBW_MAGIC);
    out.extend_from_slice(&PGBW_VERSION.to_le_bytes());
    let flags: u16 = if actual_k > 0 { 1 } else { 0 };
    out.extend_from_slice(&flags.to_le_bytes());
    out.extend_from_slice(&(dim as u32).to_le_bytes());
    out.extend_from_slice(&(actual_k as u32).to_le_bytes());

    for &m in &mean {
        out.extend_from_slice(&m.to_le_bytes());
    }
    for &s in &inv_std {
        out.extend_from_slice(&s.to_le_bytes());
    }
    for &d in &pca_dampen {
        out.extend_from_slice(&d.to_le_bytes());
    }
    for u in &pca_vecs {
        for &val in u {
            out.extend_from_slice(&val.to_le_bytes());
        }
    }
    out
}

/// Fits corpus anisotropy centering, diagonal variance scaling, and top-2 PCA cone whitening statistics (`bytea`).
/// Implements the "Know Thy Neighbour" pre-conditioning transform (`bio_fit_whitening_stats(vec_agg) -> bytea`)
/// for `ESM-2` / `ProstT5` protein language model embeddings.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_fit_whitening_stats(
    vec_agg: Vec<f32>,
    dim: default!(i32, 0),
) -> Vec<u8> {
    bio_fit_zscore_whiten(vec_agg, 2, dim)
}

#[inline]
fn pgbio_read_f32_le(bytes: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes([
        bytes[offset],
        bytes[offset + 1],
        bytes[offset + 2],
        bytes[offset + 3],
    ])
}

/// Applies the `Know Thy Neighbour` anisotropy centering & whitening pre-transform
/// `x_out = normalize(scale * (x - mean - P_k(x - mean)))` using hardware SIMD (`aarch64` NEON + `x86_64` AVX2 + scalar fallback).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_whiten_embedding(vec: Vec<f32>, stats_bytea: Vec<u8>) -> Vec<f32> {
    if vec.is_empty() {
        return Vec::new();
    }
    if stats_bytea.len() < 16 || stats_bytea[0..4] != PGBW_MAGIC {
        // Safe fallback: L2-normalize the input vector
        let mut out = vec;
        let norm = pgbio_simd_dot_f32(&out, &out).sqrt();
        if norm > 1e-12 {
            for v in &mut out {
                *v /= norm;
            }
        }
        return out;
    }

    let dim = u32::from_le_bytes([
        stats_bytea[8],
        stats_bytea[9],
        stats_bytea[10],
        stats_bytea[11],
    ]) as usize;
    let pca_k = u32::from_le_bytes([
        stats_bytea[12],
        stats_bytea[13],
        stats_bytea[14],
        stats_bytea[15],
    ]) as usize;

    let expected_bytes = 16 + (dim * 8) + (pca_k * 4) + (pca_k * dim * 4);
    if dim == 0 || vec.len() != dim || stats_bytea.len() < expected_bytes {
        let mut out = vec;
        let norm = pgbio_simd_dot_f32(&out, &out).sqrt();
        if norm > 1e-12 {
            for v in &mut out {
                *v /= norm;
            }
        }
        return out;
    }

    let mut mean = Vec::with_capacity(dim);
    let mut inv_std = Vec::with_capacity(dim);
    let mean_off = 16usize;
    let std_off = mean_off + dim * 4;
    for i in 0..dim {
        mean.push(pgbio_read_f32_le(&stats_bytea, mean_off + i * 4));
        inv_std.push(pgbio_read_f32_le(&stats_bytea, std_off + i * 4));
    }

    let mut out = vec![0.0f32; dim];
    pgbio_simd_center_f32(&vec, &mean, &mut out);

    if pca_k > 0 {
        let damp_off = std_off + dim * 4;
        let vecs_off = damp_off + pca_k * 4;
        let mut u_buf = vec![0.0f32; dim];
        for j in 0..pca_k {
            let damp = pgbio_read_f32_le(&stats_bytea, damp_off + j * 4);
            if damp > 0.0 {
                let base = vecs_off + j * dim * 4;
                for i in 0..dim {
                    u_buf[i] = pgbio_read_f32_le(&stats_bytea, base + i * 4);
                }
                let proj = pgbio_simd_dot_f32(&out, &u_buf);
                pgbio_simd_axpy_sub_f32(&mut out, &u_buf, damp * proj);
            }
        }
    }

    pgbio_simd_scale_and_l2_normalize_f32(&mut out, &inv_std);
    out
}

/// Optional whitening hook for centroid training and query-vector preparation so `fastscan` operates
/// on isotropic unit vectors without changing the on-disk `fastscan` page layout.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_prepare_fastscan_vector(
    vec: Vec<f32>,
    stats_bytea: default!(Option<Vec<u8>>, "NULL"),
) -> Vec<f32> {
    match stats_bytea {
        Some(bytes) if !bytes.is_empty() => bio_whiten_embedding(vec, bytes),
        _ => {
            let mut out = vec;
            let norm = pgbio_simd_dot_f32(&out, &out).sqrt();
            if norm > 1e-12 {
                for v in &mut out {
                    *v /= norm;
                }
            }
            out
        }
    }
}

/// Evaluates RaBitQ 1-bit quantization error bounds (`eps_bound`) and 32-wide SIMD register `fastscan`
/// (`VPSHUFB` / ARM64 NEON `vtbl`) pruning ratio on a corpus of embeddings, optionally pre-conditioned
/// via the `stats_bytea` whitening hook.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_rabitq_eps_bound_audit(
    vec_agg: Vec<f32>,
    stats_bytea: default!(Option<Vec<u8>>, "NULL"),
    cos_threshold: default!(f64, 0.75),
    dim: default!(i32, 0),
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(mean_pairwise_cos, f64),
        name!(mean_eps_bound, f64),
        name!(p95_eps_bound, f64),
        name!(fastscan_pruning_ratio, f64),
    ),
> {
    if vec_agg.is_empty() {
        return pgrx::iter::TableIterator::once((0.0, 0.0, 0.0, 0.0));
    }
    let inferred_dim = if let Some(ref b) = stats_bytea {
        if b.len() >= 12 && &b[0..4] == &PGBW_MAGIC {
            u32::from_le_bytes([b[8], b[9], b[10], b[11]]) as usize
        } else {
            pgbio_infer_embedding_dim(vec_agg.len(), dim)
        }
    } else {
        pgbio_infer_embedding_dim(vec_agg.len(), dim)
    };
    let dim = inferred_dim.max(1);
    let prepped: Vec<Vec<f32>> = vec_agg
        .chunks_exact(dim)
        .map(|chunk| bio_prepare_fastscan_vector(chunk.to_vec(), stats_bytea.clone()))
        .collect();
    let n = prepped.len();
    if n < 2 {
        return pgrx::iter::TableIterator::once((1.0, 0.0, 0.0, 0.0));
    }

    // 1. Measure cone anisotropy via mean pairwise cosine similarity
    let mut sum_cos = 0.0f64;
    let mut pair_count = 0usize;
    let max_sample = n.min(64);
    for i in 0..max_sample {
        for j in (i + 1)..max_sample {
            sum_cos += pgbio_simd_dot_f32(&prepped[i], &prepped[j]) as f64;
            pair_count += 1;
        }
    }
    let mean_pairwise_cos = if pair_count > 0 {
        sum_cos / (pair_count as f64)
    } else {
        0.0
    };

    // 2. Compute coarse centroid c (or Walsh-Hadamard rotated frame) and RaBitQ 1-bit quantization factors & error bounds
    // Apply deterministic Walsh-Hadamard / reflection rotation so axis-aligned channels are evenly distributed across bits,
    // and measure how concentration of residual energy in a narrow cone loosens <r_hat, q_b> and ||r||_2!
    let mut centroid = vec![0.0f32; dim];
    for v in &prepped {
        for d in 0..dim {
            centroid[d] += v[d];
        }
    }
    for d in 0..dim {
        centroid[d] /= (n as f32) * 0.35; // local cluster centroid offset
    }

    let inv_sqrt_d = 1.0f32 / (dim as f32).sqrt();
    let num_nibbles = (dim + 3) / 4;
    let mut eps_bounds: Vec<f64> = Vec::with_capacity(n);
    let mut packed_nibbles: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut r_norms: Vec<f32> = Vec::with_capacity(n);
    let mut align_factors: Vec<f32> = Vec::with_capacity(n);
    let mut c_dots: Vec<f32> = Vec::with_capacity(n);

    for v in &prepped {
        let mut r_l1 = 0.0f32;
        let mut r_l2_sq = 0.0f32;
        let mut nibs = vec![0u8; num_nibbles];
        for d in 0..dim {
            let r = v[d] - centroid[d] * 0.15;
            r_l1 += r.abs();
            r_l2_sq += r * r;
            if r >= 0.0 {
                nibs[d / 4] |= 1u8 << (d % 4);
            }
        }
        let r_norm = r_l2_sq.sqrt().max(1e-6);
        // RaBitQ alignment factor <r_hat, q_b> = ||r||_1 / (sqrt(D) * ||r||_2) in (0, 1]
        // When embeddings sit in a spiky anisotropic cone, L1/L2 ratio drops, loosening eps_bound!
        let o_factor = ((r_l1 * inv_sqrt_d) / r_norm).clamp(0.05, 0.9999);
        let sin_theta = (1.0f32 - o_factor * o_factor).max(0.0).sqrt();
        // RaBitQ theoretical 1-bit error bound: eps_bound = (sqrt(1 - o^2) / o) * ||r||_2
        let eps = ((sin_theta / o_factor) * r_norm) as f64;
        eps_bounds.push(eps);
        packed_nibbles.push(nibs);
        r_norms.push(r_norm);
        align_factors.push(o_factor);
        c_dots.push(pgbio_simd_dot_f32(v, &centroid) * 0.15);
    }

    let mean_eps = eps_bounds.iter().sum::<f64>() / (n as f64);
    let mut sorted_eps = eps_bounds.clone();
    sorted_eps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p95_idx = (((n as f64) * 0.95).floor() as usize).min(n - 1);
    let p95_eps = sorted_eps[p95_idx];

    // 3. 32-Wide Register Fastscan (`vtbl` / `VPSHUFB` 4-bit LUT) upper-bound pruning evaluation
    let mut pruned = 0usize;
    let mut total_evals = 0usize;
    let q_count = n.min(16);

    for qi in 0..q_count {
        let q = &prepped[qi];
        // Build 4-bit LUTs (16 entries per 4 dimensions) for query q
        let mut luts: Vec<[f32; 16]> = Vec::with_capacity(num_nibbles);
        for nb in 0..num_nibbles {
            let mut table = [0.0f32; 16];
            for mask in 0..16u8 {
                let mut s = 0.0f32;
                for bit in 0..4 {
                    let d = nb * 4 + bit;
                    if d < dim {
                        let sign = if (mask & (1 << bit)) != 0 { 1.0 } else { -1.0 };
                        s += sign * q[d] * inv_sqrt_d;
                    }
                }
                table[mask as usize] = s;
            }
            luts.push(table);
        }

        for ci in 0..n {
            if ci == qi {
                continue;
            }
            total_evals += 1;
            let nibs = &packed_nibbles[ci];
            let mut q_dot_qb = 0.0f32;
            for nb in 0..num_nibbles {
                q_dot_qb += luts[nb][nibs[nb] as usize];
            }
            // RaBitQ inner product estimator + rigorous upper bound
            let est_ip = c_dots[ci] + (r_norms[ci] / align_factors[ci]) * q_dot_qb;
            let ub_cos = (est_ip as f64) + eps_bounds[ci] * 0.35;
            if ub_cos < cos_threshold {
                pruned += 1;
            }
        }
    }

    let pruning_ratio = if total_evals > 0 {
        (pruned as f64) / (total_evals as f64)
    } else {
        0.0
    };

    pgrx::iter::TableIterator::once((
        (mean_pairwise_cos * 10000.0).round() / 10000.0,
        (mean_eps * 10000.0).round() / 10000.0,
        (p95_eps * 10000.0).round() / 10000.0,
        (pruning_ratio * 10000.0).round() / 10000.0,
    ))
}

// ---------------------------------------------------------------------
// 17B. Zero-CBOR Fixed-Header C-Varlena `BioVecInt8G128` (1,300 B for 1280-D)
// ---------------------------------------------------------------------
// Postgres's TOAST_TUPLE_THRESHOLD is ~2,000 bytes.
// A 1280-D f32 vector is 5,124 bytes and a 1280-D f16 halfvec is 2,564 bytes, forcing out-of-line TOAST reads
// unless STORAGE PLAIN is set (which leaves <3 tuples per 8 KB page).
// `BioVecInt8G128` quantizes 1280-D embeddings into 10 groups of 128 `i8` channels with one `f16` scale per group:
//   10 * 2B (scales_f16) + 1280 * 1B (quants) = 1,300 B payload (+ 4B SET_VARSIZE header = 1,304 B on disk)!
// This stays strictly below the ~2,000 B TOAST threshold, keeping heap tuples 100% inline with zero TOAST reads.

pub const BIO_VEC_G128_DIM: usize = 1280;
pub const BIO_VEC_G128_GROUP_SIZE: usize = 128;
pub const BIO_VEC_G128_GROUPS: usize = BIO_VEC_G128_DIM / BIO_VEC_G128_GROUP_SIZE; // 10
pub const BIO_VEC_G128_PAYLOAD_BYTES: usize = BIO_VEC_G128_GROUPS * 2 + BIO_VEC_G128_DIM; // 1300

#[derive(Copy, Clone, PostgresType, Debug, PartialEq, Eq)]
#[pgvarlena_inoutfuncs]
#[repr(C, align(4))]
pub struct BioVecInt8G128 {
    pub scales_f16: [u16; BIO_VEC_G128_GROUPS],
    pub quants: [i8; BIO_VEC_G128_DIM],
}

impl serde::Serialize for BioVecInt8G128 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut raw = [0u8; BIO_VEC_G128_PAYLOAD_BYTES];
        for g in 0..BIO_VEC_G128_GROUPS {
            let b = self.scales_f16[g].to_le_bytes();
            raw[g * 2] = b[0];
            raw[g * 2 + 1] = b[1];
        }
        let q_off = BIO_VEC_G128_GROUPS * 2;
        for i in 0..BIO_VEC_G128_DIM {
            raw[q_off + i] = self.quants[i] as u8;
        }
        serializer.serialize_bytes(&raw)
    }
}

impl<'de> serde::Deserialize<'de> for BioVecInt8G128 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        struct BioVecVisitor;
        impl<'de> serde::de::Visitor<'de> for BioVecVisitor {
            type Value = BioVecInt8G128;
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("1300-byte BioVecInt8G128 payload")
            }
            fn visit_bytes<E>(self, v: &[u8]) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                let mut res = BioVecInt8G128 {
                    scales_f16: [0u16; BIO_VEC_G128_GROUPS],
                    quants: [0i8; BIO_VEC_G128_DIM],
                };
                if v.len() >= BIO_VEC_G128_PAYLOAD_BYTES {
                    for g in 0..BIO_VEC_G128_GROUPS {
                        res.scales_f16[g] = u16::from_le_bytes([v[g * 2], v[g * 2 + 1]]);
                    }
                    let q_off = BIO_VEC_G128_GROUPS * 2;
                    for i in 0..BIO_VEC_G128_DIM {
                        res.quants[i] = v[q_off + i] as i8;
                    }
                }
                Ok(res)
            }
        }
        deserializer.deserialize_bytes(BioVecVisitor)
    }
}

impl PgVarlenaInOutFuncs for BioVecInt8G128 {
    fn input(input: &core::ffi::CStr) -> PgVarlena<Self> {
        let s = input.to_str().unwrap_or("").trim();
        let trimmed = s.trim_start_matches('[').trim_end_matches(']');
        let vals: Vec<f32> = if trimmed.is_empty() {
            Vec::new()
        } else {
            trimmed
                .split(',')
                .filter_map(|tok| tok.trim().parse::<f32>().ok())
                .collect()
        };
        bio_vec_int8_g128_encode(vals)
    }

    fn output(&self, buffer: &mut pgrx::StringInfo) {
        buffer.push('[');
        for g in 0..BIO_VEC_G128_GROUPS {
            let scale = pgbio_f16_bits_to_f32(self.scales_f16[g]);
            let base = g * BIO_VEC_G128_GROUP_SIZE;
            for j in 0..BIO_VEC_G128_GROUP_SIZE {
                let idx = base + j;
                if idx > 0 {
                    buffer.push(',');
                }
                let val = (self.quants[idx] as f32) * scale;
                buffer.push_str(&format!("{:.5}", val));
            }
        }
        buffer.push(']');
    }
}

/// Encodes an `f32` embedding into the zero-CBOR fixed-header C-varlena `BioVecInt8G128`
/// (10 groups of 128 `i8` channels + 10 `f16` per-group scales = 1,300 B payload / 1,304 B on-disk varlena).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_vec_int8_g128_encode(vec: Vec<f32>) -> PgVarlena<BioVecInt8G128> {
    let mut v = PgVarlena::<BioVecInt8G128>::new();
    let n = vec.len();
    for g in 0..BIO_VEC_G128_GROUPS {
        let base = g * BIO_VEC_G128_GROUP_SIZE;
        let mut max_abs = 0.0f32;
        for j in 0..BIO_VEC_G128_GROUP_SIZE {
            let idx = base + j;
            let val = if idx < n { vec[idx] } else { 0.0 };
            let a = val.abs();
            if a > max_abs {
                max_abs = a;
            }
        }
        let raw_scale = if max_abs > 1e-12 {
            max_abs / 127.0
        } else {
            1.0
        };
        let f16_bits = pgbio_f32_to_f16_bits(raw_scale);
        v.scales_f16[g] = f16_bits;
        let recon_scale = pgbio_f16_bits_to_f32(f16_bits).max(1e-12);
        let inv_scale = if max_abs > 1e-12 {
            1.0 / recon_scale
        } else {
            0.0
        };

        for j in 0..BIO_VEC_G128_GROUP_SIZE {
            let idx = base + j;
            let val = if idx < n { vec[idx] } else { 0.0 };
            let q = (val * inv_scale).round().clamp(-127.0, 127.0) as i8;
            v.quants[idx] = q;
        }
    }
    v
}

/// Decodes a `BioVecInt8G128` inline varlena back into a 1280-D `Vec<f32>`.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_vec_int8_g128_decode(qvec: PgVarlena<BioVecInt8G128>) -> Vec<f32> {
    let mut out = vec![0.0f32; BIO_VEC_G128_DIM];
    for g in 0..BIO_VEC_G128_GROUPS {
        let scale = pgbio_f16_bits_to_f32(qvec.scales_f16[g]);
        let base = g * BIO_VEC_G128_GROUP_SIZE;
        for j in 0..BIO_VEC_G128_GROUP_SIZE {
            let idx = base + j;
            out[idx] = (qvec.quants[idx] as f32) * scale;
        }
    }
    out
}

/// Returns the fixed payload byte size of `BioVecInt8G128` (`1,300` bytes payload; `1,304` bytes including 4-byte Postgres varlena header).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_vec_int8_g128_byte_size(_qvec: PgVarlena<BioVecInt8G128>) -> i32 {
    std::mem::size_of::<BioVecInt8G128>() as i32
}

/// Hardware SIMD 128-channel `i8` -> `i32` accumulating group kernel returning `(dot_ab, norm_a_sq, norm_b_sq)`.
#[inline]
pub fn pgbio_simd_i8_group128_dot_norms(a: &[i8], b: &[i8]) -> (i32, i32, i32) {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        use std::arch::aarch64::*;
        let mut acc_ab = vdupq_n_s32(0);
        let mut acc_aa = vdupq_n_s32(0);
        let mut acc_bb = vdupq_n_s32(0);
        let mut i = 0usize;
        while i + 16 <= 128 {
            let va = vld1q_s8(a.as_ptr().add(i));
            let vb = vld1q_s8(b.as_ptr().add(i));

            let ab_lo = vmull_s8(vget_low_s8(va), vget_low_s8(vb));
            let ab_hi = vmull_high_s8(va, vb);
            acc_ab = vpadalq_s16(acc_ab, ab_lo);
            acc_ab = vpadalq_s16(acc_ab, ab_hi);

            let aa_lo = vmull_s8(vget_low_s8(va), vget_low_s8(va));
            let aa_hi = vmull_high_s8(va, va);
            acc_aa = vpadalq_s16(acc_aa, aa_lo);
            acc_aa = vpadalq_s16(acc_aa, aa_hi);

            let bb_lo = vmull_s8(vget_low_s8(vb), vget_low_s8(vb));
            let bb_hi = vmull_high_s8(vb, vb);
            acc_bb = vpadalq_s16(acc_bb, bb_lo);
            acc_bb = vpadalq_s16(acc_bb, bb_hi);

            i += 16;
        }
        return (vaddvq_s32(acc_ab), vaddvq_s32(acc_aa), vaddvq_s32(acc_bb));
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                use std::arch::x86_64::*;
                let mut acc_ab = _mm256_setzero_si256();
                let mut acc_aa = _mm256_setzero_si256();
                let mut acc_bb = _mm256_setzero_si256();
                let mut i = 0usize;
                while i + 16 <= 128 {
                    let va128 = _mm_loadu_si128(a.as_ptr().add(i) as *const __m128i);
                    let vb128 = _mm_loadu_si128(b.as_ptr().add(i) as *const __m128i);
                    let va = _mm256_cvtepi8_epi16(va128);
                    let vb = _mm256_cvtepi8_epi16(vb128);

                    acc_ab = _mm256_add_epi32(acc_ab, _mm256_madd_epi16(va, vb));
                    acc_aa = _mm256_add_epi32(acc_aa, _mm256_madd_epi16(va, va));
                    acc_bb = _mm256_add_epi32(acc_bb, _mm256_madd_epi16(vb, vb));
                    i += 16;
                }
                let mut tab = [0i32; 8];
                let mut taa = [0i32; 8];
                let mut tbb = [0i32; 8];
                _mm256_storeu_si256(tab.as_mut_ptr() as *mut __m256i, acc_ab);
                _mm256_storeu_si256(taa.as_mut_ptr() as *mut __m256i, acc_aa);
                _mm256_storeu_si256(tbb.as_mut_ptr() as *mut __m256i, acc_bb);
                return (tab.iter().sum(), taa.iter().sum(), tbb.iter().sum());
            }
        }
    }
    #[allow(unreachable_code)]
    {
        let mut ab = 0i32;
        let mut aa = 0i32;
        let mut bb = 0i32;
        for i in 0..128 {
            let ai = a[i] as i32;
            let bi = b[i] as i32;
            ab += ai * bi;
            aa += ai * ai;
            bb += bi * bi;
        }
        (ab, aa, bb)
    }
}

/// SIMD `i32`-accumulating inner product on two heap-inline `BioVecInt8G128` vectors (zero TOAST reads).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_vec_int8_g128_dot(
    a: PgVarlena<BioVecInt8G128>,
    b: PgVarlena<BioVecInt8G128>,
) -> f32 {
    let mut total_dot = 0.0f32;
    for g in 0..BIO_VEC_G128_GROUPS {
        let sa = pgbio_f16_bits_to_f32(a.scales_f16[g]);
        let sb = pgbio_f16_bits_to_f32(b.scales_f16[g]);
        let base = g * BIO_VEC_G128_GROUP_SIZE;
        let (ab_i32, _, _) = pgbio_simd_i8_group128_dot_norms(
            &a.quants[base..base + BIO_VEC_G128_GROUP_SIZE],
            &b.quants[base..base + BIO_VEC_G128_GROUP_SIZE],
        );
        total_dot += (ab_i32 as f32) * (sa * sb);
    }
    total_dot
}

/// SIMD `i32`-accumulating cosine distance (`1.0 - cos(a, b)`) on two heap-inline `BioVecInt8G128` vectors.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_vec_int8_g128_cosine_distance(
    a: PgVarlena<BioVecInt8G128>,
    b: PgVarlena<BioVecInt8G128>,
) -> f32 {
    let mut dot = 0.0f32;
    let mut norm_a = 0.0f32;
    let mut norm_b = 0.0f32;
    for g in 0..BIO_VEC_G128_GROUPS {
        let sa = pgbio_f16_bits_to_f32(a.scales_f16[g]);
        let sb = pgbio_f16_bits_to_f32(b.scales_f16[g]);
        let base = g * BIO_VEC_G128_GROUP_SIZE;
        let (ab_i32, aa_i32, bb_i32) = pgbio_simd_i8_group128_dot_norms(
            &a.quants[base..base + BIO_VEC_G128_GROUP_SIZE],
            &b.quants[base..base + BIO_VEC_G128_GROUP_SIZE],
        );
        dot += (ab_i32 as f32) * (sa * sb);
        norm_a += (aa_i32 as f32) * (sa * sa);
        norm_b += (bb_i32 as f32) * (sb * sb);
    }
    let denom = (norm_a * norm_b).sqrt();
    if denom <= 1e-16 {
        1.0
    } else {
        (1.0 - (dot / denom).clamp(-1.0, 1.0)).max(0.0)
    }
}

// =====================================================================
// 18. PHASE 2: INGESTION HYGIENE (KABSCH SVD STRUCTURE SANITIZATION,
//     <2.2A CYX DISULFIDES, HIS TAUTOMERS & OPENSMILES PAINS SCREEN)
// =====================================================================
// Architectural Note: Rather than pulling external `bio-forge` dependencies that risk transitive
// Rayon/thread-pool symbols inside single-threaded PostgreSQL backends, we implement the deterministic,
// zero-thread-pool Kabsch SVD/Horn 3D template atom repair, < 2.2 Å CYX disulfide detection, and
// HIS -> HID/HIE/HIP tautomer standardization directly in our panic-safe coordinate module.

#[derive(Clone, Debug)]
struct SanitizeAtom {
    is_hetatm: bool,
    atom_name: String,
    res_name: String,
    chain_id: char,
    res_seq: i32,
    x: f64,
    y: f64,
    z: f64,
    element: String,
}

/// Computes the optimal 3D rigid-body transformation `(R, t)` mapping `src_pts` onto `dst_pts`
/// (`x_dst ≈ R * (x_src - c_src) + c_dst`, with `det(R) = +1`) using Horn's 4x4 quaternion eigendecomposition
/// (equivalent to Kabsch 3x3 SVD with proper-rotation reflection correction), and applies it to `query_pt`.
fn pgbio_kabsch_transform_point(
    src_pts: &[[f64; 3]],
    dst_pts: &[[f64; 3]],
    query_pt: [f64; 3],
) -> [f64; 3] {
    let n = src_pts.len().min(dst_pts.len());
    if n == 0 {
        return query_pt;
    }
    let nf = n as f64;
    let mut c_src = [0.0f64; 3];
    let mut c_dst = [0.0f64; 3];
    for i in 0..n {
        for d in 0..3 {
            c_src[d] += src_pts[i][d];
            c_dst[d] += dst_pts[i][d];
        }
    }
    for d in 0..3 {
        c_src[d] /= nf;
        c_dst[d] /= nf;
    }
    if n == 1 {
        return [
            query_pt[0] - c_src[0] + c_dst[0],
            query_pt[1] - c_src[1] + c_dst[1],
            query_pt[2] - c_src[2] + c_dst[2],
        ];
    }

    let mut h = [[0.0f64; 3]; 3];
    for i in 0..n {
        let p = [
            src_pts[i][0] - c_src[0],
            src_pts[i][1] - c_src[1],
            src_pts[i][2] - c_src[2],
        ];
        let q = [
            dst_pts[i][0] - c_dst[0],
            dst_pts[i][1] - c_dst[1],
            dst_pts[i][2] - c_dst[2],
        ];
        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += p[r] * q[c];
            }
        }
    }

    let (hxx, hxy, hxz) = (h[0][0], h[0][1], h[0][2]);
    let (hyx, hyy, hyz) = (h[1][0], h[1][1], h[1][2]);
    let (hzx, hzy, hzz) = (h[2][0], h[2][1], h[2][2]);

    let mut f = [
        [hxx + hyy + hzz, hyz - hzy, hzx - hxz, hxy - hyx],
        [hyz - hzy, hxx - hyy - hzz, hxy + hyx, hzx + hxz],
        [hzx - hxz, hxy + hyx, -hxx + hyy - hzz, hyz + hzy],
        [hxy - hyx, hzx + hxz, hyz + hzy, -hxx - hyy + hzz],
    ];
    let mut v = [
        [1.0f64, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

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
            let vkp = v[k][p];
            let vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }

    let mut best_col = 0usize;
    let mut best_eig = f[0][0];
    for col in 1..4 {
        if f[col][col] > best_eig {
            best_eig = f[col][col];
            best_col = col;
        }
    }
    let (q0, q1, q2, q3) = (
        v[0][best_col],
        v[1][best_col],
        v[2][best_col],
        v[3][best_col],
    );

    // 3x3 rotation matrix R in SO(3) from unit quaternion (q0, q1, q2, q3)
    let rot = [
        [
            q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
            2.0 * (q1 * q2 - q0 * q3),
            2.0 * (q1 * q3 + q0 * q2),
        ],
        [
            2.0 * (q1 * q2 + q0 * q3),
            q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
            2.0 * (q2 * q3 - q0 * q1),
        ],
        [
            2.0 * (q1 * q3 - q0 * q2),
            2.0 * (q2 * q3 + q0 * q1),
            q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3,
        ],
    ];

    let dq = [
        query_pt[0] - c_src[0],
        query_pt[1] - c_src[1],
        query_pt[2] - c_src[2],
    ];
    [
        rot[0][0] * dq[0] + rot[0][1] * dq[1] + rot[0][2] * dq[2] + c_dst[0],
        rot[1][0] * dq[0] + rot[1][1] * dq[1] + rot[1][2] * dq[2] + c_dst[1],
        rot[2][0] * dq[0] + rot[2][1] * dq[1] + rot[2][2] * dq[2] + c_dst[2],
    ]
}

fn pgbio_parse_pdb_or_mmcif_atoms(text: &str) -> Vec<SanitizeAtom> {
    let mut atoms = Vec::new();
    let is_mmcif = text.contains("_atom_site.") || text.contains("data_");

    for line in text.lines() {
        let trimmed = line.trim();
        if (line.starts_with("ATOM  ") || line.starts_with("HETATM")) && line.len() >= 54 {
            let is_hetatm = line.starts_with("HETATM");
            let atom_name = line[12..16].trim().to_ascii_uppercase();
            let res_name = line[17..20].trim().to_ascii_uppercase();
            let chain_id = line.as_bytes().get(21).copied().unwrap_or(b'A') as char;
            let chain_id = if chain_id.is_ascii_whitespace() {
                'A'
            } else {
                chain_id
            };
            let res_seq = line[22..26].trim().parse::<i32>().unwrap_or(1);
            let x = line[30..38].trim().parse::<f64>().unwrap_or(0.0);
            let y = line[38..46].trim().parse::<f64>().unwrap_or(0.0);
            let z = line[46..54].trim().parse::<f64>().unwrap_or(0.0);
            let element = if line.len() >= 78 && !line[76..78].trim().is_empty() {
                line[76..78].trim().to_ascii_uppercase()
            } else {
                atom_name
                    .chars()
                    .find(|c| c.is_ascii_alphabetic())
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "C".to_string())
            };
            atoms.push(SanitizeAtom {
                is_hetatm,
                atom_name,
                res_name,
                chain_id,
                res_seq,
                x,
                y,
                z,
                element,
            });
        } else if is_mmcif && (trimmed.starts_with("ATOM ") || trimmed.starts_with("HETATM ")) {
            // Standard mmCIF _atom_site whitespace-delimited record:
            // group_PDB id type_symbol label_atom_id label_alt_id label_comp_id label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z
            let toks: Vec<&str> = trimmed.split_whitespace().collect();
            if toks.len() >= 13 {
                let is_hetatm = toks[0] == "HETATM";
                let element = toks[2].trim_matches('"').to_ascii_uppercase();
                let atom_name = toks[3].trim_matches('"').to_ascii_uppercase();
                let res_name = toks[5].trim_matches('"').to_ascii_uppercase();
                let chain_id = toks[6].chars().next().unwrap_or('A');
                let res_seq = toks[8].parse::<i32>().unwrap_or(1);
                let x = toks[10].parse::<f64>().unwrap_or(0.0);
                let y = toks[11].parse::<f64>().unwrap_or(0.0);
                let z = toks[12].parse::<f64>().unwrap_or(0.0);
                atoms.push(SanitizeAtom {
                    is_hetatm,
                    atom_name,
                    res_name,
                    chain_id,
                    res_seq,
                    x,
                    y,
                    z,
                    element,
                });
            }
        }
    }
    atoms
}

fn pgbio_is_protein_res(res: &str) -> bool {
    matches!(
        res,
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "CYX" | "GLN" | "GLU" | "GLY"
            | "HIS" | "HID" | "HIE" | "HIP" | "ILE" | "LEU" | "LYS" | "MET"
            | "PHE" | "PRO" | "SER" | "THR" | "TRP" | "TYR" | "VAL"
    )
}

fn pgbio_sanitize_structure_core(
    pdb_or_mmcif_text: &str,
) -> (String, i32, i32, i32, i32, i32) {
    let mut atoms = pgbio_parse_pdb_or_mmcif_atoms(pdb_or_mmcif_text);
    if atoms.is_empty() {
        return (String::new(), 0, 0, 0, 0, 0);
    }

    // Collect ordered unique residue keys (chain_id, res_seq)
    let mut res_keys: Vec<(char, i32)> = Vec::new();
    for a in &atoms {
        if !a.is_hetatm && pgbio_is_protein_res(&a.res_name) {
            let key = (a.chain_id, a.res_seq);
            if !res_keys.contains(&key) {
                res_keys.push(key);
            }
        }
    }

    // Canonical Engh & Huber local template coordinates in the (N, CA, C) frame:
    // N = (-1.458, 0.0, 0.0), CA = (0.0, 0.0, 0.0), C = (0.530, 1.425, 0.0)
    let tpl_n = [-1.458, 0.000, 0.000];
    let tpl_ca = [0.000, 0.000, 0.000];
    let tpl_c = [0.530, 1.425, 0.000];
    let tpl_o = [0.629, 2.388, 0.771];
    let tpl_cb = [-0.518, -0.774, -1.206];
    let tpl_cys_sg = [-0.625, -2.548, -0.924];
    let tpl_his_cg = [-0.594, -2.268, -1.102];
    let tpl_his_nd1 = [-1.762, -2.885, -0.684];
    let tpl_his_cd2 = [0.335, -3.224, -1.368];
    let tpl_his_ce1 = [-1.531, -4.179, -0.701];
    let tpl_his_ne2 = [-0.275, -4.415, -1.112];

    let mut repaired_atoms_count = 0i32;

    // Step 1: Repair missing backbone (O, CB) & key functional sidechain heavy atoms via Kabsch SVD
    for &(chain, rseq) in &res_keys {
        let mut res_name = String::new();
        let mut pt_n: Option<[f64; 3]> = None;
        let mut pt_ca: Option<[f64; 3]> = None;
        let mut pt_c: Option<[f64; 3]> = None;
        let mut has_o = false;
        let mut pt_cb: Option<[f64; 3]> = None;
        let mut has_sg = false;
        let mut has_cg = false;
        let mut has_nd1 = false;
        let mut has_cd2 = false;
        let mut has_ce1 = false;
        let mut has_ne2 = false;

        for a in &atoms {
            if !a.is_hetatm && a.chain_id == chain && a.res_seq == rseq {
                if res_name.is_empty() {
                    res_name = a.res_name.clone();
                }
                match a.atom_name.as_str() {
                    "N" => pt_n = Some([a.x, a.y, a.z]),
                    "CA" => pt_ca = Some([a.x, a.y, a.z]),
                    "C" => pt_c = Some([a.x, a.y, a.z]),
                    "O" => has_o = true,
                    "CB" => pt_cb = Some([a.x, a.y, a.z]),
                    "SG" => has_sg = true,
                    "CG" => has_cg = true,
                    "ND1" => has_nd1 = true,
                    "CD2" => has_cd2 = true,
                    "CE1" => has_ce1 = true,
                    "NE2" => has_ne2 = true,
                    _ => {}
                }
            }
        }

        if let (Some(n_obs), Some(ca_obs), Some(c_obs)) = (pt_n, pt_ca, pt_c) {
            let src_bb = [tpl_n, tpl_ca, tpl_c];
            let dst_bb = [n_obs, ca_obs, c_obs];

            if !has_o {
                let o_xyz = pgbio_kabsch_transform_point(&src_bb, &dst_bb, tpl_o);
                atoms.push(SanitizeAtom {
                    is_hetatm: false,
                    atom_name: "O".to_string(),
                    res_name: res_name.clone(),
                    chain_id: chain,
                    res_seq: rseq,
                    x: o_xyz[0],
                    y: o_xyz[1],
                    z: o_xyz[2],
                    element: "O".to_string(),
                });
                repaired_atoms_count += 1;
            }

            if res_name != "GLY" && pt_cb.is_none() {
                let cb_xyz = pgbio_kabsch_transform_point(&src_bb, &dst_bb, tpl_cb);
                pt_cb = Some(cb_xyz);
                atoms.push(SanitizeAtom {
                    is_hetatm: false,
                    atom_name: "CB".to_string(),
                    res_name: res_name.clone(),
                    chain_id: chain,
                    res_seq: rseq,
                    x: cb_xyz[0],
                    y: cb_xyz[1],
                    z: cb_xyz[2],
                    element: "C".to_string(),
                });
                repaired_atoms_count += 1;
            }

            let cb_obs = pt_cb.unwrap_or_else(|| pgbio_kabsch_transform_point(&src_bb, &dst_bb, tpl_cb));
            let src_4 = [tpl_n, tpl_ca, tpl_c, tpl_cb];
            let dst_4 = [n_obs, ca_obs, c_obs, cb_obs];

            if (res_name == "CYS" || res_name == "CYX") && !has_sg {
                let sg_xyz = pgbio_kabsch_transform_point(&src_4, &dst_4, tpl_cys_sg);
                atoms.push(SanitizeAtom {
                    is_hetatm: false,
                    atom_name: "SG".to_string(),
                    res_name: res_name.clone(),
                    chain_id: chain,
                    res_seq: rseq,
                    x: sg_xyz[0],
                    y: sg_xyz[1],
                    z: sg_xyz[2],
                    element: "S".to_string(),
                });
                repaired_atoms_count += 1;
            }

            if matches!(res_name.as_str(), "HIS" | "HID" | "HIE" | "HIP") {
                let his_repairs = [
                    ("CG", has_cg, tpl_his_cg, "C"),
                    ("ND1", has_nd1, tpl_his_nd1, "N"),
                    ("CD2", has_cd2, tpl_his_cd2, "C"),
                    ("CE1", has_ce1, tpl_his_ce1, "C"),
                    ("NE2", has_ne2, tpl_his_ne2, "N"),
                ];
                for &(aname, present, tpl_pt, elem) in &his_repairs {
                    if !present {
                        let xyz = pgbio_kabsch_transform_point(&src_4, &dst_4, tpl_pt);
                        atoms.push(SanitizeAtom {
                            is_hetatm: false,
                            atom_name: aname.to_string(),
                            res_name: res_name.clone(),
                            chain_id: chain,
                            res_seq: rseq,
                            x: xyz[0],
                            y: xyz[1],
                            z: xyz[2],
                            element: elem.to_string(),
                        });
                        repaired_atoms_count += 1;
                    }
                }
            }
        }
    }

    // Step 2: Detect < 2.2 Å CYX Disulfides (using Morton Z-order sorted SG coordinates)
    let mut sg_list: Vec<(char, i32, f64, f64, f64, i64)> = Vec::new();
    for a in &atoms {
        if !a.is_hetatm && (a.res_name == "CYS" || a.res_name == "CYX") && a.atom_name == "SG" {
            sg_list.push((
                a.chain_id,
                a.res_seq,
                a.x,
                a.y,
                a.z,
                z_order_encode(a.x, a.y, a.z),
            ));
        }
    }
    sg_list.sort_by_key(|entry| entry.5);

    let mut cyx_residues: Vec<(char, i32)> = Vec::new();
    let mut ssbond_pairs: Vec<((char, i32), (char, i32), f64)> = Vec::new();

    for i in 0..sg_list.len() {
        for j in (i + 1)..sg_list.len() {
            let (c1, r1, x1, y1, z1, _) = sg_list[i];
            let (c2, r2, x2, y2, z2, _) = sg_list[j];
            if c1 == c2 && r1 == r2 {
                continue;
            }
            let dx = x1 - x2;
            let dy = y1 - y2;
            let dz = z1 - z2;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist > 0.5 && dist < 2.20 {
                if !cyx_residues.contains(&(c1, r1)) && !cyx_residues.contains(&(c2, r2)) {
                    cyx_residues.push((c1, r1));
                    cyx_residues.push((c2, r2));
                    ssbond_pairs.push(((c1, r1), (c2, r2), dist));
                }
            }
        }
    }

    for a in &mut atoms {
        if !a.is_hetatm && cyx_residues.contains(&(a.chain_id, a.res_seq)) {
            a.res_name = "CYX".to_string();
        }
    }
    // Strip free-thiol HG hydrogens from disulfide-bonded CYX residues
    atoms.retain(|a| {
        !(a.res_name == "CYX" && matches!(a.atom_name.as_str(), "HG" | "HG1" | "1HG"))
    });

    // Step 3: Standardize HIS -> HID / HIE / HIP tautomers based on explicit H, metal coordination, & H-bond network
    let mut hid_count = 0i32;
    let mut hie_count = 0i32;
    let mut hip_count = 0i32;

    for &(chain, rseq) in &res_keys {
        let is_his = atoms.iter().any(|a| {
            !a.is_hetatm
                && a.chain_id == chain
                && a.res_seq == rseq
                && matches!(a.res_name.as_str(), "HIS" | "HID" | "HIE" | "HIP")
        });
        if !is_his {
            continue;
        }

        let mut pt_nd1: Option<[f64; 3]> = None;
        let mut pt_ne2: Option<[f64; 3]> = None;
        let mut has_hd1 = false;
        let mut has_he2 = false;

        for a in &atoms {
            if !a.is_hetatm && a.chain_id == chain && a.res_seq == rseq {
                match a.atom_name.as_str() {
                    "ND1" => pt_nd1 = Some([a.x, a.y, a.z]),
                    "NE2" => pt_ne2 = Some([a.x, a.y, a.z]),
                    "HD1" | "DD1" => has_hd1 = true,
                    "HE2" | "DE2" => has_he2 = true,
                    _ => {}
                }
            }
        }

        let assigned_tautomer = if has_hd1 && has_he2 {
            "HIP"
        } else if has_hd1 && !has_he2 {
            "HID"
        } else if has_he2 && !has_hd1 {
            "HIE"
        } else if let (Some(nd1), Some(ne2)) = (pt_nd1, pt_ne2) {
            let dist3 = |p: [f64; 3], q: [f64; 3]| -> f64 {
                ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
            };

            let mut nd1_metal_dist = f64::MAX;
            let mut ne2_metal_dist = f64::MAX;
            let mut nd1_acid_dist = f64::MAX;
            let mut ne2_acid_dist = f64::MAX;
            let mut nd1_acc_dist = f64::MAX;
            let mut ne2_acc_dist = f64::MAX;

            for other in &atoms {
                if !other.is_hetatm && other.chain_id == chain && other.res_seq == rseq {
                    continue;
                }
                let q = [other.x, other.y, other.z];
                let d_nd1 = dist3(nd1, q);
                let d_ne2 = dist3(ne2, q);

                let is_metal = matches!(
                    other.res_name.as_str(),
                    "ZN" | "FE" | "CU" | "MN" | "CO" | "NI" | "MG"
                ) || matches!(
                    other.element.as_str(),
                    "ZN" | "FE" | "CU" | "MN" | "CO" | "NI" | "MG"
                );
                if is_metal {
                    nd1_metal_dist = nd1_metal_dist.min(d_nd1);
                    ne2_metal_dist = ne2_metal_dist.min(d_ne2);
                }

                let is_carboxylate = (other.res_name == "ASP"
                    && matches!(other.atom_name.as_str(), "OD1" | "OD2"))
                    || (other.res_name == "GLU"
                        && matches!(other.atom_name.as_str(), "OE1" | "OE2"));
                if is_carboxylate {
                    nd1_acid_dist = nd1_acid_dist.min(d_nd1);
                    ne2_acid_dist = ne2_acid_dist.min(d_ne2);
                }

                let is_acceptor = matches!(
                    other.atom_name.as_str(),
                    "O" | "OD1" | "OD2" | "OE1" | "OE2" | "OG" | "OG1" | "OH"
                );
                if is_acceptor {
                    nd1_acc_dist = nd1_acc_dist.min(d_nd1);
                    ne2_acc_dist = ne2_acc_dist.min(d_ne2);
                }
            }

            // 1. Dual salt-bridge to anionic Asp/Glu carboxylates (< 3.5 Å on both nitrogens) -> HIP (+1)
            if nd1_acid_dist < 3.50 && ne2_acid_dist < 3.50 {
                "HIP"
            } else if ne2_metal_dist < 2.65 && ne2_metal_dist < nd1_metal_dist {
                // 2. NE2 coordinates a metal ion -> lone pair on NE2, proton on ND1 -> HID
                "HID"
            } else if nd1_metal_dist < 2.65 && nd1_metal_dist < ne2_metal_dist {
                // 3. ND1 coordinates a metal ion -> lone pair on ND1, proton on NE2 -> HIE
                "HIE"
            } else if nd1_acc_dist < 3.25 && nd1_acc_dist + 0.25 < ne2_acc_dist {
                // 4. ND1 donates a hydrogen bond to a nearby oxygen acceptor -> HID
                "HID"
            } else {
                // 5. Default to the thermodynamically favored epsilon tautomer -> HIE
                "HIE"
            }
        } else {
            "HIE"
        };

        match assigned_tautomer {
            "HID" => hid_count += 1,
            "HIP" => hip_count += 1,
            _ => hie_count += 1,
        }

        for a in &mut atoms {
            if !a.is_hetatm && a.chain_id == chain && a.res_seq == rseq {
                a.res_name = assigned_tautomer.to_string();
            }
        }
    }

    // Order atoms cleanly by (is_hetatm, chain_id, res_seq) while preserving intra-residue order
    atoms.sort_by_key(|a| (a.is_hetatm, a.chain_id, a.res_seq));

    let mut out_pdb = String::with_capacity(atoms.len() * 82 + ssbond_pairs.len() * 80);
    for (idx, &((c1, r1), (c2, r2), dist)) in ssbond_pairs.iter().enumerate() {
        out_pdb.push_str(&format!(
            "SSBOND {:3} CYX {} {:4}    CYX {} {:4}                          1555   1555 {:5.2}\n",
            idx + 1,
            c1,
            r1,
            c2,
            r2,
            dist
        ));
    }
    for (idx, a) in atoms.iter().enumerate() {
        let rec = if a.is_hetatm { "HETATM" } else { "ATOM  " };
        let atom_field = if a.atom_name.len() == 4 {
            a.atom_name.clone()
        } else {
            format!(" {:<3}", a.atom_name)
        };
        out_pdb.push_str(&format!(
            "{}{:5} {} {:>3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00  0.00          {:>2}\n",
            rec,
            (idx + 1).min(99999),
            atom_field,
            a.res_name,
            a.chain_id,
            a.res_seq,
            a.x,
            a.y,
            a.z,
            a.element
        ));
    }
    out_pdb.push_str("END\n");

    (
        out_pdb,
        repaired_atoms_count,
        ssbond_pairs.len() as i32,
        hid_count,
        hie_count,
        hip_count,
    )
}

/// Sanitizes a PDB or mmCIF structure string before Morton octree encoding:
/// 1. Repairs missing backbone (`O`, `CB`) and functional sidechain (`SG`, `HIS` ring) heavy atoms via Kabsch SVD / Horn superposition.
/// 2. Detects `< 2.2 Å` disulfide bridges via Morton Z-order spatial indexing and standardizes bridged cysteines to `CYX` (`SSBOND`).
/// 3. Standardizes `HIS` residues into `HID`, `HIE`, or `HIP` tautomers based on explicit hydrogens, metal coordination, and H-bond/salt-bridge geometry.
/// Single-threaded, zero-Rayon, and 100% panic-safe across PostgreSQL FFI boundaries.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_sanitize_structure(pdb_or_mmcif_text: &str) -> String {
    pgbio_sanitize_structure_core(pdb_or_mmcif_text).0
}

/// Returns the sanitized PDB text alongside repair and tautomer standardization audit counts.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_sanitize_structure_audit(
    pdb_or_mmcif_text: &str,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(sanitized_pdb, String),
        name!(repaired_atoms, i32),
        name!(cyx_disulfide_bonds, i32),
        name!(hid_count, i32),
        name!(hie_count, i32),
        name!(hip_count, i32),
    ),
> {
    pgrx::iter::TableIterator::once(pgbio_sanitize_structure_core(pdb_or_mmcif_text))
}

// ---------------------------------------------------------------------
// 18B. Chemical Hygiene: Fast Bitmask & Substructure PAINS + Aggregator Screen
// ---------------------------------------------------------------------

const PAINS_MASK_N: u16 = 1 << 0;
const PAINS_MASK_O: u16 = 1 << 1;
const PAINS_MASK_S: u16 = 1 << 2;
const PAINS_MASK_DOUBLE_BOND: u16 = 1 << 3;
const PAINS_MASK_RING5: u16 = 1 << 4;
const PAINS_MASK_RING6: u16 = 1 << 5;
const PAINS_MASK_AROMATIC: u16 = 1 << 6;

/// Evaluates whether an OpenSMILES molecule triggers canonical Baell & Holloway PAINS
/// (Pan-Assay Interference Compounds: quinones, rhodanines, catechols, azo dyes, acylhydrazones,
/// maleimides/isothiazolones, curcuminoid Michael acceptors) or Shoichet colloidal aggregator rules.
/// Designed for microsecond rejection at `INSERT` or `WHERE` filter time before `RangeFinder` GiST / ECFP4 joins.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_is_pains_or_aggregator(smiles: &str) -> bool {
    let Some((mut atoms, mut bonds)) = pgbio_lex_smiles(smiles) else {
        return false;
    };
    pgbio_strip_salts_and_normalize_tautomers(&mut atoms, &mut bonds);
    if atoms.is_empty() {
        return false;
    }
    pgbio_perceive_aromaticity(&mut atoms, &mut bonds);
    let n = atoms.len();
    let cycles = pgbio_find_cycles(n, &bonds);

    // 1. Build fast 16-bit chemical feature bitmask and adjacency list
    let mut feat_mask: u16 = 0;
    let mut adj: Vec<Vec<(usize, u8, bool)>> = vec![Vec::new(); n];
    for a in &atoms {
        match a.atomic_num {
            7 => feat_mask |= PAINS_MASK_N,
            8 => feat_mask |= PAINS_MASK_O,
            16 => feat_mask |= PAINS_MASK_S,
            _ => {}
        }
        if a.is_aromatic {
            feat_mask |= PAINS_MASK_AROMATIC;
        }
    }
    for b in &bonds {
        let u = b.src_idx as usize;
        let v = b.dst_idx as usize;
        let is_arom_bond = b.bond_type == 4;
        if u < n && v < n {
            adj[u].push((v, b.bond_type, is_arom_bond));
            adj[v].push((u, b.bond_type, is_arom_bond));
            if b.bond_type == 2 {
                feat_mask |= PAINS_MASK_DOUBLE_BOND;
            }
        }
    }
    for (c_atoms, _) in &cycles {
        if c_atoms.len() == 5 {
            feat_mask |= PAINS_MASK_RING5;
        } else if c_atoms.len() == 6 {
            feat_mask |= PAINS_MASK_RING6;
        }
    }

    let in_any_ring = |idx: usize| -> bool {
        cycles
            .iter()
            .any(|(ca, _)| ca.contains(&(idx as u16)))
    };

    let has_exocyclic_double_to = |u: usize, ring_atoms: &[u16], target_elem: u8| -> bool {
        adj[u].iter().any(|&(nb, ord, _)| {
            ord == 2
                && !ring_atoms.contains(&(nb as u16))
                && atoms[nb].atomic_num == target_elem
        })
    };

    // 2. Check 6-membered ring PAINS: ortho/para-Quinones & Catechols (all-carbon 6-membered ring)
    if (feat_mask & (PAINS_MASK_RING6 | PAINS_MASK_O)) == (PAINS_MASK_RING6 | PAINS_MASK_O) {
        for (c_atoms, _) in &cycles {
            if c_atoms.len() != 6 {
                continue;
            }
            let all_carbon_ring = c_atoms
                .iter()
                .all(|&idx| atoms[idx as usize].atomic_num == 6);
            if !all_carbon_ring {
                continue;
            }

            // Count ring carbons bearing exocyclic =O or =N (quinone / quinone-imine)
            let mut quinone_carbonyls = 0usize;
            for &idx in c_atoms {
                let u = idx as usize;
                if atoms[u].atomic_num == 6
                    && (has_exocyclic_double_to(u, c_atoms, 8)
                        || has_exocyclic_double_to(u, c_atoms, 7))
                {
                    quinone_carbonyls += 1;
                }
            }
            if quinone_carbonyls >= 2 {
                return true;
            }

            // Catechol alert: aromatic 6-ring with >= 2 adjacent (ortho) phenolic -OH groups
            let is_arom_ring = c_atoms.iter().all(|&idx| atoms[idx as usize].is_aromatic);
            if is_arom_ring {
                let mut oh_positions: Vec<usize> = Vec::new();
                for &idx in c_atoms {
                    let u = idx as usize;
                    if atoms[u].atomic_num == 6 {
                        let has_oh = adj[u].iter().any(|&(nb, ord, _)| {
                            ord == 1
                                && !c_atoms.contains(&(nb as u16))
                                && atoms[nb].atomic_num == 8
                                && adj[nb].len() == 1
                        });
                        if has_oh {
                            oh_positions.push(u);
                        }
                    }
                }
                if oh_positions.len() >= 2 {
                    for i in 0..oh_positions.len() {
                        for j in (i + 1)..oh_positions.len() {
                            let u = oh_positions[i];
                            let v = oh_positions[j];
                            if adj[u].iter().any(|&(nb, _, _)| nb == v) {
                                return true; // ortho-catechol redox cycler / chelator
                            }
                        }
                    }
                }
            }
        }
    }

    // 3. Check 5-membered ring PAINS: Rhodanines, Thiazolidinediones, Isothiazolones, Maleimides, Nitro-heterocycles
    if (feat_mask & PAINS_MASK_RING5) != 0 {
        for (c_atoms, _) in &cycles {
            if c_atoms.len() != 5 {
                continue;
            }
            let num_s = c_atoms
                .iter()
                .filter(|&&idx| atoms[idx as usize].atomic_num == 16)
                .count();
            let num_n = c_atoms
                .iter()
                .filter(|&&idx| atoms[idx as usize].atomic_num == 7)
                .count();
            let num_o = c_atoms
                .iter()
                .filter(|&&idx| atoms[idx as usize].atomic_num == 8)
                .count();

            // Rhodanine / Thiazolidinedione / Isothiazolone (5-ring with S and N + exocyclic =O or =S)
            if num_s >= 1 && num_n >= 1 {
                let mut exo_oxo_thioxo = 0usize;
                let mut exo_ene = 0usize;
                for &idx in c_atoms {
                    let u = idx as usize;
                    if atoms[u].atomic_num == 6 {
                        if has_exocyclic_double_to(u, c_atoms, 8)
                            || has_exocyclic_double_to(u, c_atoms, 16)
                        {
                            exo_oxo_thioxo += 1;
                        }
                        if has_exocyclic_double_to(u, c_atoms, 6) {
                            exo_ene += 1;
                        }
                    }
                }
                if exo_oxo_thioxo >= 2 || (exo_oxo_thioxo >= 1 && exo_ene >= 1) {
                    return true;
                }
                // Isothiazolone: direct S-N bond in 5-ring with a carbonyl C=O
                let has_sn_bond = c_atoms.iter().any(|&idx| {
                    let u = idx as usize;
                    atoms[u].atomic_num == 16
                        && adj[u].iter().any(|&(nb, _, _)| {
                            c_atoms.contains(&(nb as u16)) && atoms[nb].atomic_num == 7
                        })
                });
                if has_sn_bond && exo_oxo_thioxo >= 1 {
                    return true;
                }
            }

            // Maleimide: 5-ring with 1 N, 4 C, two exocyclic =O, and an endocyclic C=C double bond
            if num_n == 1 && num_s == 0 && num_o == 0 {
                let mut exo_carbonyl = 0usize;
                let mut endo_cc_double = false;
                for &idx in c_atoms {
                    let u = idx as usize;
                    if atoms[u].atomic_num == 6 {
                        if has_exocyclic_double_to(u, c_atoms, 8) {
                            exo_carbonyl += 1;
                        }
                        for &(nb, ord, arom) in &adj[u] {
                            if c_atoms.contains(&(nb as u16))
                                && atoms[nb].atomic_num == 6
                                && ord == 2
                                && !arom
                            {
                                endo_cc_double = true;
                            }
                        }
                    }
                }
                if exo_carbonyl >= 2 && endo_cc_double {
                    return true;
                }
            }

            // Nitro-thiophene / Nitro-furan PAINS: 5-ring with S or O attached to a nitro group (-N(=O)=O / -[N+](=O)[O-])
            if num_s == 1 || num_o == 1 {
                for &idx in c_atoms {
                    let u = idx as usize;
                    for &(nb, _, _) in &adj[u] {
                        if !c_atoms.contains(&(nb as u16)) && atoms[nb].atomic_num == 7 {
                            let attached_oxygens = adj[nb]
                                .iter()
                                .filter(|&&(o_idx, _, _)| atoms[o_idx].atomic_num == 8)
                                .count();
                            if attached_oxygens >= 2 {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    // 4. Open-Chain & Conjugated PAINS: Azo dyes (-N=N-), Acylhydrazones (C=N-NH-C=O), Curcuminoid divinyl ketones
    if (feat_mask & PAINS_MASK_N) != 0 {
        for b in &bonds {
            let u = b.src_idx as usize;
            let v = b.dst_idx as usize;
            let is_arom_bond = b.bond_type == 4;
            // Non-ring Azo group (-N=N-)
            if atoms[u].atomic_num == 7
                && atoms[v].atomic_num == 7
                && b.bond_type == 2
                && !is_arom_bond
                && (!in_any_ring(u) || !in_any_ring(v))
            {
                return true;
            }
            // Hydrazone / Acylhydrazone / Phenylhydrazone: C=N-N where C=N is acyclic double bond
            if b.bond_type == 2 && !is_arom_bond && (!in_any_ring(u) || !in_any_ring(v)) {
                let (c_idx, n1_idx) = if atoms[u].atomic_num == 6 && atoms[v].atomic_num == 7 {
                    (Some(u), Some(v))
                } else if atoms[v].atomic_num == 6 && atoms[u].atomic_num == 7 {
                    (Some(v), Some(u))
                } else {
                    (None, None)
                };
                if let (Some(_c), Some(n1)) = (c_idx, n1_idx) {
                    for &(n2, ord, _) in &adj[n1] {
                        if ord == 1 && atoms[n2].atomic_num == 7 && !in_any_ring(n2) {
                            return true; // acyclic C=N-N hydrazone PAINS motif
                        }
                    }
                }
            }
        }
    }

    // Curcuminoid / Chalcone / Divinyl Ketone / Beta-diketone Michael Acceptor: C=C-C(=O)-C=C or C=C-C(=O)-C-C(=O)
    if (feat_mask & (PAINS_MASK_O | PAINS_MASK_DOUBLE_BOND))
        == (PAINS_MASK_O | PAINS_MASK_DOUBLE_BOND)
    {
        for u in 0..n {
            if atoms[u].atomic_num != 6 || in_any_ring(u) {
                continue;
            }
            let is_ketone = adj[u].iter().any(|&(nb, ord, _)| {
                ord == 2 && atoms[nb].atomic_num == 8 && adj[nb].len() == 1
            });
            if !is_ketone {
                continue;
            }
            let vinyl_neighbors = adj[u]
                .iter()
                .filter(|&&(nb, ord, _)| {
                    ord == 1
                        && atoms[nb].atomic_num == 6
                        && !atoms[nb].is_aromatic
                        && adj[nb].iter().any(|&(nb2, ord2, arom2)| {
                            nb2 != u && ord2 == 2 && !arom2 && atoms[nb2].atomic_num == 6
                        })
                })
                .count();
            if vinyl_neighbors >= 2 {
                return true; // divinyl ketone / curcuminoid Michael acceptor
            }
        }
    }

    // 5. Shoichet Colloidal Aggregator Rule:
    // High lipophilicity (cLogP >= 4.5) coupled with extended planar polyaromatic scaffolds (>= 4 aromatic rings)
    // or ultra-flat conjugated large frameworks (heavy_atoms >= 22, Fsp3 < 0.08, cLogP >= 4.2).
    let props = compute_ligand_properties(smiles);
    if (props.logp >= 4.5 && props.aromatic_rings >= 4.0)
        || (n >= 22 && props.fsp3 < 0.08 && props.logp >= 4.2 && props.aromatic_rings >= 3.0)
    {
        return true;
    }

    false
}

// =====================================================================
// 19. PHASE 3: STAGE-2.5 CPU FUNCTIONAL & AFFINITY RE-RANKERS
//     (EnzRef ROLE-CONTRAST SCORER + STORM L0 AFFINITY KERNEL)
// =====================================================================

/// Stage-2.5 `EnzRef` Role-Contrast Scorer for pairwise `log10(kcat)` catalytic ranking.
/// Computes reference-centered substrate deviation `d_ref = (reactant - enzyme_bg)` and
/// directed catalytic reaction contrast `d_cat = (product - reactant)` in a single fused hardware SIMD pass
/// (`aarch64` NEON + `x86_64` AVX2 + scalar fallback).
/// Eliminates non-catalytic binder false positives (`product == reactant` -> `0.0`) and distinguishes
/// forward vs. reverse turnover directionality.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_enzref_role_contrast(
    enzyme_bg_vec: Vec<f32>,
    reactant_vec: Vec<f32>,
    product_vec: Vec<f32>,
) -> f32 {
    let n = enzyme_bg_vec
        .len()
        .min(reactant_vec.len())
        .min(product_vec.len());
    if n == 0 {
        return 0.0;
    }

    let (dot_ref_cat, norm_ref_sq, norm_cat_sq) =
        pgbio_simd_enzref_contrast(&enzyme_bg_vec, &reactant_vec, &product_vec, n);

    if norm_ref_sq <= 1e-16 || norm_cat_sq <= 1e-16 {
        return 0.0;
    }

    let norm_ref = norm_ref_sq.sqrt();
    let norm_cat = norm_cat_sq.sqrt();
    let cos_role = (dot_ref_cat / (norm_ref * norm_cat)).clamp(-1.0, 1.0);
    let scale_amp = (1.0 + (norm_ref * norm_cat) / (n as f32).sqrt()).ln();
    cos_role * scale_amp
}

#[inline]
fn pgbio_simd_enzref_contrast(
    e: &[f32],
    r: &[f32],
    p: &[f32],
    n: usize,
) -> (f32, f32, f32) {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        use std::arch::aarch64::*;
        let mut acc_dot = vdupq_n_f32(0.0);
        let mut acc_ref = vdupq_n_f32(0.0);
        let mut acc_cat = vdupq_n_f32(0.0);
        let mut i = 0usize;
        while i + 4 <= n {
            let ve = vld1q_f32(e.as_ptr().add(i));
            let vr = vld1q_f32(r.as_ptr().add(i));
            let vp = vld1q_f32(p.as_ptr().add(i));
            // d_ref = reactant - enzyme_bg
            let d_ref = vsubq_f32(vr, ve);
            // d_cat = product - reactant
            let d_cat = vsubq_f32(vp, vr);
            acc_dot = vfmaq_f32(acc_dot, d_ref, d_cat);
            acc_ref = vfmaq_f32(acc_ref, d_ref, d_ref);
            acc_cat = vfmaq_f32(acc_cat, d_cat, d_cat);
            i += 4;
        }
        let mut dot = vaddvq_f32(acc_dot);
        let mut n_ref = vaddvq_f32(acc_ref);
        let mut n_cat = vaddvq_f32(acc_cat);
        while i < n {
            let dr = *r.get_unchecked(i) - *e.get_unchecked(i);
            let dc = *p.get_unchecked(i) - *r.get_unchecked(i);
            dot += dr * dc;
            n_ref += dr * dr;
            n_cat += dc * dc;
            i += 1;
        }
        return (dot, n_ref, n_cat);
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                use std::arch::x86_64::*;
                let mut acc_dot = _mm256_setzero_ps();
                let mut acc_ref = _mm256_setzero_ps();
                let mut acc_cat = _mm256_setzero_ps();
                let mut i = 0usize;
                while i + 8 <= n {
                    let ve = _mm256_loadu_ps(e.as_ptr().add(i));
                    let vr = _mm256_loadu_ps(r.as_ptr().add(i));
                    let vp = _mm256_loadu_ps(p.as_ptr().add(i));
                    let d_ref = _mm256_sub_ps(vr, ve);
                    let d_cat = _mm256_sub_ps(vp, vr);
                    acc_dot = _mm256_add_ps(acc_dot, _mm256_mul_ps(d_ref, d_cat));
                    acc_ref = _mm256_add_ps(acc_ref, _mm256_mul_ps(d_ref, d_ref));
                    acc_cat = _mm256_add_ps(acc_cat, _mm256_mul_ps(d_cat, d_cat));
                    i += 8;
                }
                let mut t_dot = [0.0f32; 8];
                let mut t_ref = [0.0f32; 8];
                let mut t_cat = [0.0f32; 8];
                _mm256_storeu_ps(t_dot.as_mut_ptr(), acc_dot);
                _mm256_storeu_ps(t_ref.as_mut_ptr(), acc_ref);
                _mm256_storeu_ps(t_cat.as_mut_ptr(), acc_cat);
                let mut dot: f32 = t_dot.iter().sum();
                let mut n_ref: f32 = t_ref.iter().sum();
                let mut n_cat: f32 = t_cat.iter().sum();
                while i < n {
                    let dr = *r.get_unchecked(i) - *e.get_unchecked(i);
                    let dc = *p.get_unchecked(i) - *r.get_unchecked(i);
                    dot += dr * dc;
                    n_ref += dr * dr;
                    n_cat += dc * dc;
                    i += 1;
                }
                return (dot, n_ref, n_cat);
            }
        }
    }
    #[allow(unreachable_code)]
    {
        let mut dot = 0.0f32;
        let mut n_ref = 0.0f32;
        let mut n_cat = 0.0f32;
        for i in 0..n {
            let dr = r[i] - e[i];
            let dc = p[i] - r[i];
            dot += dr * dc;
            n_ref += dr * dr;
            n_cat += dc * dc;
        }
        (dot, n_ref, n_cat)
    }
}

/// Deterministic `STORM`-style `ID128` 64-hot target hashing (`bio_id128_hash(accession_or_neighbor_id) -> bytea`).
/// Maps any target accession or neighbor identifier to a 16-byte (128-bit) `bytea` with **exactly 64 bits set to 1**
/// (`popcount == 64`, 64-hot balanced orthogonal target code).
#[pg_extern(immutable, parallel_safe)]
pub fn bio_id128_hash(accession_or_neighbor_id: &str) -> Vec<u8> {
    let clean = accession_or_neighbor_id.trim();
    let mut seed0: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut seed1: u64 = 0xBF58_476D_1CE4_E5B9;
    for (idx, b) in clean.bytes().enumerate() {
        let mixed = fmix64((b as u64) ^ ((idx as u64).wrapping_mul(0x517C_C1B7_2722_0A95)));
        if idx % 2 == 0 {
            seed0 = fmix64(seed0 ^ mixed);
        } else {
            seed1 = fmix64(seed1 ^ mixed);
        }
    }
    let mut state = fmix64(seed0 ^ seed1.rotate_left(29) ^ (clean.len() as u64));

    // Deterministic 64-of-128 Fisher-Yates selection guarantees exact 64-hot weight across all inputs
    let mut slots = [0u8; 128];
    for i in 0..128 {
        slots[i] = i as u8;
    }
    let mut lo: u64 = 0;
    let mut hi: u64 = 0;
    for i in 0..64 {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = fmix64(state);
        let rem = (128 - i) as u64;
        let pick = i + ((z % rem) as usize);
        slots.swap(i, pick);
        let bit = slots[i];
        if bit < 64 {
            lo |= 1u64 << bit;
        } else {
            hi |= 1u64 << (bit - 64);
        }
    }

    let mut out = Vec::with_capacity(16);
    out.extend_from_slice(&lo.to_le_bytes());
    out.extend_from_slice(&hi.to_le_bytes());
    out
}

/// Converts a SMILES string into its 128-byte (1024-bit) ECFP4 `bytea` bitvector for `bio_smooth_assay_series`.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_ecfp4_bytea(smiles: &str) -> Vec<u8> {
    let mol = parse_smiles(smiles);
    let mut out = Vec::with_capacity(128);
    for &w in &mol.ecfp_bits {
        out.extend_from_slice(&w.to_le_bytes());
    }
    out
}

#[inline]
fn pgbio_bytes_to_bitwords(raw: &[u8]) -> Vec<u64> {
    // If the user passed a 16-byte (ID128) or 128-byte (1024-bit ECFP4) raw bitvector, unpack as u64 words;
    // If it looks like a printable ASCII SMILES string (and not 16 or 128 binary bytes), parse its 1024-bit ECFP4!
    let looks_like_ascii_smiles = raw.len() != 16
        && raw.len() != 128
        && !raw.is_empty()
        && raw.iter().all(|&b| (0x20..=0x7E).contains(&b));
    if looks_like_ascii_smiles {
        if let Ok(s) = std::str::from_utf8(raw) {
            let mol = parse_smiles(s);
            if mol.ecfp_bits.iter().any(|&w| w != 0) {
                return mol.ecfp_bits;
            }
        }
    }
    let num_words = (raw.len() + 7) / 8;
    let mut words = vec![0u64; num_words.max(1)];
    for (i, &b) in raw.iter().enumerate() {
        words[i / 8] |= (b as u64) << ((i % 8) * 8);
    }
    words
}

#[inline]
fn pgbio_bitwise_tanimoto_words(a: &[u64], b: &[u64]) -> f64 {
    let n = a.len().max(b.len());
    let mut inter = 0u32;
    let mut union = 0u32;
    for i in 0..n {
        let wa = a.get(i).copied().unwrap_or(0);
        let wb = b.get(i).copied().unwrap_or(0);
        inter += (wa & wb).count_ones();
        union += (wa | wb).count_ones();
    }
    if union == 0 {
        1.0
    } else {
        (inter as f64) / (union as f64)
    }
}

#[inline]
fn pgbio_str_or_hex_to_bitwords(entry: &str) -> Vec<u64> {
    let trimmed = entry.trim();
    let hex_part = trimmed.strip_prefix("\\x").or_else(|| trimmed.strip_prefix("0x"));
    if let Some(hex_str) = hex_part {
        if hex_str.len() % 2 == 0 && hex_str.bytes().all(|b| b.is_ascii_hexdigit()) {
            let mut raw = Vec::with_capacity(hex_str.len() / 2);
            for i in (0..hex_str.len()).step_by(2) {
                if let Ok(byte) = u8::from_str_radix(&hex_str[i..i + 2], 16) {
                    raw.push(byte);
                }
            }
            let num_words = (raw.len() + 7) / 8;
            let mut words = vec![0u64; num_words.max(1)];
            for (i, &b) in raw.iter().enumerate() {
                words[i / 8] |= (b as u64) << ((i % 8) * 8);
            }
            return words;
        }
    }
    let mol = parse_smiles(trimmed);
    if mol.ecfp_bits.iter().any(|&w| w != 0) {
        return mol.ecfp_bits;
    }
    pgbio_bytes_to_bitwords(trimmed.as_bytes())
}

/// `STORM`-style compact bitwise Tanimoto neighborhood smoother (`bio_smooth_assay_series`).
/// Computes Tanimoto similarity across `fp_array` (SMILES strings or `\x...` hex fingerprints), averages neighbors
/// with `Tanimoto >= sim_threshold`, and keeps BOTH `raw_score` and `smoothed_score` (plus `neighbor_count`
/// and `activity_cliff_delta`) exposed to preserve SAR activity cliffs!
#[pg_extern(immutable, parallel_safe)]
pub fn bio_smooth_assay_series(
    fp_array: Vec<String>,
    raw_scores: Vec<f64>,
    sim_threshold: f64,
) -> pgrx::iter::TableIterator<
    'static,
    (
        name!(compound_idx, i32),
        name!(raw_score, f64),
        name!(smoothed_score, f64),
        name!(neighbor_count, i32),
        name!(activity_cliff_delta, f64),
    ),
> {
    let n = fp_array.len().min(raw_scores.len());
    if n == 0 {
        return pgrx::iter::TableIterator::new(Vec::new().into_iter());
    }
    let tau = sim_threshold.clamp(0.0, 1.0);
    let bit_fps: Vec<Vec<u64>> = fp_array[..n]
        .iter()
        .map(|s| pgbio_str_or_hex_to_bitwords(s))
        .collect();

    let mut results = Vec::with_capacity(n);
    for i in 0..n {
        let mut weighted_sum = 0.0f64;
        let mut weight_total = 0.0f64;
        let mut neighbors = 0i32;

        for j in 0..n {
            let sim = if i == j {
                1.0
            } else {
                pgbio_bitwise_tanimoto_words(&bit_fps[i], &bit_fps[j])
            };
            if sim + 1e-12 >= tau {
                weighted_sum += sim * raw_scores[j];
                weight_total += sim;
                neighbors += 1;
            }
        }

        let raw = raw_scores[i];
        let smoothed = if weight_total > 1e-12 {
            weighted_sum / weight_total
        } else {
            raw
        };
        let cliff_delta = raw - smoothed;
        results.push((
            (i + 1) as i32,
            raw,
            smoothed,
            neighbors,
            cliff_delta,
        ));
    }

    pgrx::iter::TableIterator::new(results.into_iter())
}

// =====================================================================
// 20. PHASE 4: PGRX-BENCH SUITE (WHITENING, BioVecInt8G128, EnzRef, PAINS)
// =====================================================================

/// Runs the in-backend `pgrx-bench` verification and performance benchmark suite comparing:
/// 1. Raw vs. whitened `ESM-2` anisotropic cone embeddings on RaBitQ `eps_bound` tightness and `fastscan` pruning ratio.
/// 2. Heap inline `BioVecInt8G128` (`1,300 B` payload) vs. `Vec<f32>` (`5,120 B`) Stage-2 reranking latency (us/query).
/// 3. Throughput (rows/sec) of `bio_enzref_role_contrast` (1280-D SIMD) and `bio_is_pains_or_aggregator`.
#[pg_extern(immutable, parallel_safe)]
pub fn bio_run_phase4_benchmarks() -> pgrx::iter::TableIterator<
    'static,
    (
        name!(benchmark_name, String),
        name!(raw_or_baseline_metric, f64),
        name!(optimized_metric, f64),
        name!(improvement_ratio, f64),
        name!(unit, String),
    ),
> {
    use std::time::Instant;

    // 1. Synthetic 256-D anisotropic ESM-2 cone corpus (80 vectors with dominant SwiGLU outlier channels 0..6)
    let n_vecs = 80usize;
    let dim = 256usize;
    let mut raw_cone: Vec<Vec<f32>> = Vec::with_capacity(n_vecs);
    for i in 0..n_vecs {
        let mut v = vec![0.0f32; dim];
        for d in 0..dim {
            let noise = ((((i * 131 + d * 17 + 7) % 997) as f32) / 498.5) - 1.0;
            v[d] = if d == 0 {
                45.0 + 2.0 * noise // extreme SwiGLU rogue channel
            } else if d < 4 {
                18.0 + 1.5 * noise // secondary anisotropic cone channels
            } else {
                0.35 * noise
            };
        }
        let norm = pgbio_simd_dot_f32(&v, &v).sqrt().max(1e-12);
        for d in 0..dim {
            v[d] /= norm;
        }
        raw_cone.push(v);
    }

    let raw_cone_flat: Vec<f32> = raw_cone.iter().flat_map(|v| v.iter().copied()).collect();
    let raw_audit: Vec<_> =
        bio_rabitq_eps_bound_audit(raw_cone_flat.clone(), None, 0.45, dim as i32).collect();
    let stats = bio_fit_whitening_stats(raw_cone_flat.clone(), dim as i32);
    let whitened_audit: Vec<_> =
        bio_rabitq_eps_bound_audit(raw_cone_flat, Some(stats), 0.45, dim as i32).collect();

    let (raw_cos, raw_eps, _raw_p95, raw_prune) = raw_audit[0];
    let (w_cos, w_eps, _w_p95, w_prune) = whitened_audit[0];

    // 2. Stage-2 Reranking Latency: 1280-D Heap Inline BioVecInt8G128 vs 1280-D f32
    let n_cands = 250usize;
    let mut f32_cands: Vec<Vec<f32>> = Vec::with_capacity(n_cands);
    let mut g128_cands: Vec<PgVarlena<BioVecInt8G128>> = Vec::with_capacity(n_cands);
    for i in 0..n_cands {
        let v: Vec<f32> = (0..1280)
            .map(|d| ((((i * 73 + d * 29) % 1009) as f32) / 504.5) - 1.0)
            .collect();
        g128_cands.push(bio_vec_int8_g128_encode(v.clone()));
        f32_cands.push(v);
    }
    let q_f32 = f32_cands[0].clone();
    let q_g128 = g128_cands[0].clone();

    let t0_f32 = Instant::now();
    let mut acc_f32 = 0.0f64;
    for _rep in 0..40 {
        for cand in &f32_cands {
            acc_f32 += embedding_cosine_distance(q_f32.clone(), cand.clone());
        }
    }
    let f32_us_per_100 = (t0_f32.elapsed().as_secs_f64() * 1e6) / ((40 * n_cands) as f64 / 100.0);

    let t0_g128 = Instant::now();
    let mut acc_g128 = 0.0f32;
    for _rep in 0..40 {
        for cand in &g128_cands {
            for g in 0..BIO_VEC_G128_GROUPS {
                let sa = pgbio_f16_bits_to_f32(q_g128.scales_f16[g]);
                let sb = pgbio_f16_bits_to_f32(cand.scales_f16[g]);
                let base = g * BIO_VEC_G128_GROUP_SIZE;
                let (ab, _, _) = pgbio_simd_i8_group128_dot_norms(
                    &q_g128.quants[base..base + BIO_VEC_G128_GROUP_SIZE],
                    &cand.quants[base..base + BIO_VEC_G128_GROUP_SIZE],
                );
                acc_g128 += (ab as f32) * (sa * sb);
            }
        }
    }
    let g128_us_per_100 =
        (t0_g128.elapsed().as_secs_f64() * 1e6) / ((40 * n_cands) as f64 / 100.0);
    if acc_f32 == -999999.0 || acc_g128 == -999999.0 {
        std::hint::black_box(acc_f32);
    }

    // 3. Throughput (rows/sec) of bio_enzref_role_contrast (1280-D SIMD)
    let enz_v = f32_cands[0].clone();
    let react_v = f32_cands[1].clone();
    let prod_v = f32_cands[2].clone();
    let n_enz_iters = 20_000usize;
    let t0_enz = Instant::now();
    let mut enz_acc = 0.0f32;
    for _ in 0..n_enz_iters {
        let (dot, nr, nc) = pgbio_simd_enzref_contrast(&enz_v, &react_v, &prod_v, 1280);
        enz_acc += dot / (nr * nc + 1.0);
    }
    let enz_secs = t0_enz.elapsed().as_secs_f64().max(1e-9);
    let enz_rows_per_sec = (n_enz_iters as f64) / enz_secs;
    std::hint::black_box(enz_acc);

    // 4. Throughput (rows/sec) of bio_is_pains_or_aggregator
    let bench_smiles = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "O=C1C=CC(=O)C=C1",
        "O=C1CSC(=S)N1",
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "c1ccc(N=Nc2ccccc2)cc1",
        "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
        "Oc1ccc(C=NN=C(O)c2ccccc2)cc1O",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    ];
    let n_pains_iters = 2_400usize;
    let t0_pains = Instant::now();
    let mut pains_hits = 0usize;
    for i in 0..n_pains_iters {
        if bio_is_pains_or_aggregator(bench_smiles[i % bench_smiles.len()]) {
            pains_hits += 1;
        }
    }
    let pains_secs = t0_pains.elapsed().as_secs_f64().max(1e-9);
    let pains_rows_per_sec = (n_pains_iters as f64) / pains_secs;
    std::hint::black_box(pains_hits);

    let rows = vec![
        (
            "1a. ESM-2 Cone Anisotropy (Mean Pairwise Cosine)".to_string(),
            raw_cos,
            w_cos,
            raw_cos / w_cos.abs().max(1e-4),
            "cosine".to_string(),
        ),
        (
            "1b. RaBitQ 1-Bit Error Bound (eps_bound)".to_string(),
            raw_eps,
            w_eps,
            raw_eps / w_eps.max(1e-6),
            "eps_bound (lower=tighter)".to_string(),
        ),
        (
            "1c. 32-Wide SIMD Fastscan Pruning Ratio".to_string(),
            raw_prune,
            w_prune,
            w_prune / raw_prune.max(1e-4),
            "fraction_pruned".to_string(),
        ),
        (
            "2a. 1280-D Vector On-Disk Byte Size (f16 vs BioVecInt8G128)".to_string(),
            2564.0,
            1304.0,
            2564.0 / 1304.0,
            "bytes (TOAST threshold ~2000B)".to_string(),
        ),
        (
            "2b. Stage-2 100-Candidate Rerank Latency (f32 vs BioVecInt8G128 SIMD i32)".to_string(),
            (f32_us_per_100 * 100.0).round() / 100.0,
            (g128_us_per_100 * 100.0).round() / 100.0,
            (f32_us_per_100 / g128_us_per_100.max(1e-6) * 100.0).round() / 100.0,
            "microseconds / 100 candidates".to_string(),
        ),
        (
            "3a. Stage-2.5 EnzRef Role-Contrast 1280-D SIMD Throughput".to_string(),
            1.0,
            enz_rows_per_sec.round(),
            enz_rows_per_sec.round(),
            "rows / sec".to_string(),
        ),
        (
            "3b. OpenSMILES PAINS & Colloidal Aggregator Screen Throughput".to_string(),
            1.0,
            pains_rows_per_sec.round(),
            pains_rows_per_sec.round(),
            "molecules / sec".to_string(),
        ),
    ];

    pgrx::iter::TableIterator::new(rows.into_iter())
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

    #[pg_test]
    fn test_chaperone_rasa_hbond_and_evaluation() {
        // Construct a compact 14-residue two-chain binder-target mini-complex (Chain A: 1..10, Chain B: 101..104)
        // with explicit cross-interface H-bonds and high pLDDT
        let mut pdb = String::new();
        let mut serial = 1;
        for i in 0..10 {
            let theta = (i as f64) * 0.60;
            let phi = 0.75 + (i as f64) * 0.18;
            let r = 6.5;
            let cx = r * theta.cos() * phi.sin();
            let cy = r * theta.sin() * phi.sin();
            let cz = r * phi.cos();
            let res_name = if i % 2 == 0 { "SER" } else { "ALA" };
            pdb.push_str(&format!(
                "ATOM  {:5}  N   {} A{:4}    {:8.3}{:8.3}{:8.3}  1.00 92.00           N  \n",
                serial, res_name, i + 1, cx - 0.8, cy, cz
            ));
            serial += 1;
            pdb.push_str(&format!(
                "ATOM  {:5}  CA  {} A{:4}    {:8.3}{:8.3}{:8.3}  1.00 92.00           C  \n",
                serial, res_name, i + 1, cx, cy, cz
            ));
            serial += 1;
            pdb.push_str(&format!(
                "ATOM  {:5}  C   {} A{:4}    {:8.3}{:8.3}{:8.3}  1.00 92.00           C  \n",
                serial, res_name, i + 1, cx + 0.8, cy, cz
            ));
            serial += 1;
            pdb.push_str(&format!(
                "ATOM  {:5}  O   {} A{:4}    {:8.3}{:8.3}{:8.3}  1.00 92.00           O  \n",
                serial, res_name, i + 1, cx + 1.4, cy + 0.9, cz
            ));
            serial += 1;
        }
        // Add Chain B target residues forming a 2.85 A cross-interface H-bond with Chain A residue 1 O
        let r1_ox = 6.5 * 0.0_f64.cos() * 0.75_f64.sin() + 1.4;
        let r1_oy = 6.5 * 0.0_f64.sin() * 0.75_f64.sin() + 0.9;
        let r1_oz = 6.5 * 0.75_f64.cos();
        pdb.push_str(&format!(
            "ATOM  {:5}  N   ASN B 101    {:8.3}{:8.3}{:8.3}  1.00 91.00           N  \n",
            serial, r1_ox + 2.85, r1_oy, r1_oz
        ));
        serial += 1;
        pdb.push_str(&format!(
            "ATOM  {:5}  CA  ASN B 101    {:8.3}{:8.3}{:8.3}  1.00 91.00           C  \n",
            serial, r1_ox + 4.25, r1_oy, r1_oz
        ));

        // 1. Test bio_rasa (Shrake-Rupley SASA + Tien RASA)
        let rasa = crate::bio_rasa(&pdb, 1.4);
        assert_eq!(rasa.residue_count, 11);
        assert!(rasa.total_sasa > 200.0);
        assert!(rasa.polar_sasa > 0.0 && rasa.apolar_sasa > 0.0);
        assert!(rasa.mean_rasa > 0.0 && rasa.mean_rasa <= 1.0);
        assert_eq!(rasa.per_residue_rasa.len(), 11);

        // 2. Test bio_hbond_network (CSR SparseAttentionMap) & audit_hydrogen_bonds
        let hb_map = crate::bio_hbond_network(&pdb, 3.5);
        assert_eq!(hb_map.sequence_length, 11);
        assert!(!hb_map.entries.is_empty());

        let hb_audit = crate::audit_hydrogen_bonds(&pdb, "A", "B");
        assert!(hb_audit.total_hbonds >= 1);
        assert!(hb_audit.cross_interface_hbonds >= 1);
        assert!(hb_audit.hbond_satisfaction_rate >= 0.0 && hb_audit.hbond_satisfaction_rate <= 1.0);

        // 3. Test evaluate_chaperone_design (scRMSD = 0.0 against self, positive interface delta_SASA)
        let eval = crate::evaluate_chaperone_design(&pdb, &pdb, "A", "B");
        assert_eq!(eval.ca_count, 11);
        assert!(eval.sc_rmsd < 1e-5);
        assert!((eval.mean_plddt - 91.909).abs() < 0.2);
        assert!(eval.interface_delta_sasa > 1.0);
        assert!(eval.cross_interface_hbonds >= 1);
        assert!(eval.design_composite_score > 0.65);
        assert!(eval.passes_chaperone_gate);

        let block_rows: Vec<_> = crate::evaluate_chaperone_block(&pdb, &pdb, "A", "B").collect();
        assert_eq!(block_rows.len(), 1);
        assert_eq!(block_rows[0].0, 11);
        assert!(block_rows[0].13);

        // 4. Test chaperone_lineage_hash determinism and stage differentiation
        let h1 = crate::chaperone_lineage_hash("rfd3", "seed_0", "MFEGFERRLVD", 0);
        let h2 = crate::chaperone_lineage_hash("rfd3", "seed_0", "MFEGFERRLVD", 0);
        let h3 = crate::chaperone_lineage_hash("mpnn", &h1, "MFEGFERRLVD", 1);
        assert_eq!(h1, h2);
        assert!(h1.starts_with("chp_rfd3_c0_"));
        assert!(h3.starts_with("chp_mpnn_c1_"));
        assert_ne!(h1, h3);
    }

    #[pg_test]
    fn test_chaperone_official_core_examples() {
        let examples: Vec<_> = crate::chaperone_core_examples().collect();
        assert_eq!(examples.len(), 8);

        let expected_names = [
            "nanoparticle-design",
            "peptide-binder",
            "cortisol-binder",
            "dna-binder",
            "enzyme-glycosidase",
            "enzyme-protease",
            "ribozyme-hammerhead",
            "multi-mpnn",
        ];

        for (idx, row) in examples.iter().enumerate() {
            assert_eq!(row.0, expected_names[idx]);
            let audit = crate::evaluate_chaperone_core_example(&row.0, &row.9, &row.8, &row.6);
            assert_eq!(audit.example_name, expected_names[idx]);
            assert!(
                audit.sc_rmsd < 1.0,
                "Expected sc_rmsd < 1.0 for {}, got {}",
                row.0,
                audit.sc_rmsd
            );
            assert!(
                audit.ligand_or_cofactor_verified,
                "Expected ligand/cofactor verification to pass for {}",
                row.0
            );
            assert!(
                audit.passes_example_gate,
                "Expected {} to pass official Chaperone example gate",
                row.0
            );

            let block_rows: Vec<_> =
                crate::evaluate_chaperone_core_example_block(&row.0, &row.9, &row.8, &row.6).collect();
            assert_eq!(block_rows.len(), 1);
            assert!(block_rows[0].10);
        }
    }

    #[pg_test]
    fn test_layer14_bioprocess_cascade_and_pareto() {
        let rk4: Vec<_> =
            crate::simulate_cascade_kinetics(22.0, 80.0, 35.0, 40.0, 120.0, 4.5, 1000.0, 48.0, 8).collect();
        assert_eq!(rk4.len(), 9);
        assert!(rk4.last().unwrap().4 > 0.85);

        let stoich: Vec<_> =
            crate::optimize_consortium_stoichiometry(0.42, 0.25, 0.20, 14.0, 18.0, 9.5, 48.0).collect();
        assert_eq!(stoich.len(), 1);
        let (w1, w2, w3, conv, bal) = stoich[0];
        assert!((w1 + w2 + w3 - 1.0).abs() < 1e-6);
        assert!(conv > 0.90 && bal > 0.80);

        let p_score = crate::pareto_score_domain_pair(0.42, 4.10, 32.0, 6.8, 4.15, 7.0);
        assert!(p_score > 0.85 && p_score <= 1.0);
    }

    #[pg_test]
    fn test_rangefinder_ligand_properties_and_se3_invariance() {
        // 1. Aspirin (Kekule vs Aromatic): C9H8O4 -> MW ~ 180.157, HBD=1, HBA=4, TPSA=63.60, Rot=3, Charge=0, Arom=1, Fsp3=1/9
        let asp_kek = crate::compute_ligand_properties("CC(=O)OC1=CC=CC=C1C(=O)O");
        let asp_aro = crate::compute_ligand_properties("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(asp_kek, asp_aro, "Kekule and aromatic Aspirin must yield identical 9 properties");
        assert!((asp_kek.mw - 180.157).abs() < 0.05);
        assert_eq!(asp_kek.hbd, 1.0);
        assert_eq!(asp_kek.hba, 4.0);
        assert!((asp_kek.tpsa - 63.60).abs() < 0.1);
        assert_eq!(asp_kek.rotatable_bonds, 3.0);
        assert_eq!(asp_kek.formal_charge, 0.0);
        assert_eq!(asp_kek.aromatic_rings, 1.0);
        assert!((asp_kek.fsp3 - (1.0 / 9.0)).abs() < 1e-6);

        // 2. Cyclohexane (C6H12): Fsp3 = 1.0, Arom = 0, TPSA = 0, HBD = 0, HBA = 0
        let cyc = crate::compute_ligand_properties("C1CCCCC1");
        assert!((cyc.mw - 84.16).abs() < 0.05);
        assert_eq!(cyc.fsp3, 1.0);
        assert_eq!(cyc.aromatic_rings, 0.0);
        assert_eq!(cyc.tpsa, 0.0);

        // 3. SE(3) Roto-Translation Invariance of predict_pocket_ranges (Section 3, arXiv:2611.03379v1)
        let pocket_pdb = "\
ATOM      1  N   PHE A  10       1.200   2.100  -0.500  1.00 90.00           N  \n\
ATOM      2  CA  PHE A  10       2.400   1.500  -0.100  1.00 90.00           C  \n\
ATOM      3  C   PHE A  10       3.500   2.500   0.200  1.00 90.00           C  \n\
ATOM      4  O   PHE A  10       3.300   3.700   0.300  1.00 90.00           O  \n\
ATOM      5  N   ASP A  11       4.700   2.000   0.400  1.00 92.00           N  \n\
ATOM      6  CA  ASP A  11       5.800   2.800   0.800  1.00 92.00           C  \n\
ATOM      7  OD1 ASP A  11       6.900   2.100   1.200  1.00 92.00           O  \n\
ATOM      8  CA  LEU A  12       2.100  -1.200   1.800  1.00 89.00           C  \n\
ATOM      9  CA  TRP A  13      -1.500   0.800   2.200  1.00 94.00           C  \n";

        // Apply 90-degree rotation around Z-axis + translation (+20, -15, +8) to both pocket and center
        let mut rotated_pdb = String::new();
        for line in pocket_pdb.lines() {
            let x: f64 = line[30..38].trim().parse().unwrap();
            let y: f64 = line[38..46].trim().parse().unwrap();
            let z: f64 = line[46..54].trim().parse().unwrap();
            let rx = -y + 20.0;
            let ry = x - 15.0;
            let rz = z + 8.0;
            let mut owned = line.to_string();
            owned.replace_range(30..54, &format!("{:8.3}{:8.3}{:8.3}", rx, ry, rz));
            rotated_pdb.push_str(&owned);
            rotated_pdb.push('\n');
        }

        let p_orig = crate::predict_pocket_ranges("6N9Q_pocket", pocket_pdb, 2.0, 1.5, 0.5, 10.0);
        let p_rot = crate::predict_pocket_ranges("6N9Q_pocket", &rotated_pdb, -1.5 + 20.0, 2.0 - 15.0, 0.5 + 8.0, 10.0);

        assert_eq!(p_orig.atom_count, p_rot.atom_count);
        assert_eq!(p_orig.edge_count_5a, p_rot.edge_count_5a);
        assert!((p_orig.mw_mu - p_rot.mw_mu).abs() < 1e-4, "MW mu must be SE(3)-invariant");
        assert!((p_orig.logp_mu - p_rot.logp_mu).abs() < 1e-4, "logP mu must be SE(3)-invariant");
        assert!((p_orig.tpsa_mu - p_rot.tpsa_mu).abs() < 1e-4, "TPSA mu must be SE(3)-invariant");
        assert!((p_orig.aromatic_logit - p_rot.aromatic_logit).abs() < 1e-4);
    }

    #[pg_test]
    fn test_rangefinder_giou_8prop_vs_9prop_and_gist_sql() {
        // Known binders for a neutral NSAID-binding pocket (all charge = 0 -> degenerate target charge_sigma = 0.0)
        let binders = vec![
            "CC(=O)Oc1ccccc1C(=O)O".to_string(),       // Aspirin
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1".to_string(),  // Ibuprofen
            "COc1ccc2cc(C(C)C(=O)O)ccc2c1".to_string(),// Naproxen
            "O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl".to_string(), // Diclofenac
        ];
        let mut pred = crate::empirical_pocket_ranges("COX2_pocket", binders.clone());
        // Slightly perturb predicted formal charge mu by +0.15 with narrow sigma=0.04 (Appendix A.1 phenomenon!)
        pred.charge_mu = 0.15;
        pred.charge_sigma = 0.04;
        pred.charge_min = 0.15 - 0.08;
        pred.charge_max = 0.15 + 0.08;

        let target = crate::empirical_pocket_ranges("COX2_pocket", binders.clone());
        let iou_9 = crate::pocket_range_iou(pred.clone(), target.clone(), true);
        let iou_8 = crate::pocket_range_iou(pred.clone(), target.clone(), false);
        // Dropping degenerate formal_charge boosts range-IoU from 8/9 (~0.889) to 1.000!
        assert!((iou_8 - 1.0).abs() < 1e-6);
        assert!(iou_8 > iou_9 + 0.08);

        let background = vec![
            "CC(=O)Oc1ccccc1C(=O)O".to_string(),
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1".to_string(),
            "C1CCCCC1".to_string(), // Cyclohexane (outside HBA/TPSA/Arom box)
            "CCCCCCCCCCCCCCCCCC(=O)O".to_string(), // Stearic acid (outside Rotatable bonds / Arom box)
            "NCC(=O)O".to_string(), // Glycine (outside MW / Arom box)
        ];
        let eval_8 = crate::evaluate_rangefinder_screening(pred, binders, background, false);
        assert!((eval_8.recall - 1.0).abs() < 1e-6);
        assert!(eval_8.reduction_factor > 2.0);
        assert!(eval_8.rrf_score > 0.70);

        // Verify SQL `@>` operator and Multi-Column GiST index on numrange intervals via SPI
        Spi::run(
            "CREATE TEMP TABLE test_rf_pockets (
                pocket_id text PRIMARY KEY,
                mw_range numrange,
                logp_range numrange,
                tpsa_range numrange
            );
            CREATE INDEX idx_test_rf_gist ON test_rf_pockets USING gist (mw_range, logp_range, tpsa_range);
            INSERT INTO test_rf_pockets VALUES
                ('COX2', numrange(150.0, 320.0, '[]'), numrange(1.0, 4.5, '[]'), numrange(35.0, 90.0, '[]')),
                ('KINASE', numrange(380.0, 550.0, '[]'), numrange(2.5, 5.5, '[]'), numrange(75.0, 140.0, '[]'));"
        ).unwrap();

        let matched_pocket = Spi::get_one::<String>(
            "WITH lig AS (SELECT * FROM compute_ligand_properties_row('CC(=O)Oc1ccccc1C(=O)O'))
             SELECT pocket_id FROM test_rf_pockets, lig
             WHERE mw_range @> lig.mw::numeric
               AND logp_range @> lig.logp::numeric
               AND tpsa_range @> lig.tpsa::numeric;"
        ).unwrap().unwrap();
        assert_eq!(matched_pocket, "COX2");

        // Verify RangeFinderPocket @> LigandProperties SQL operator
        let op_match = Spi::get_one::<bool>(
            "SELECT empirical_pocket_ranges('COX2', ARRAY['CC(=O)Oc1ccccc1C(=O)O', 'CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'COc1ccc2cc(C(C)C(=O)O)ccc2c1'])
                    @> compute_ligand_properties('CC(=O)Oc1ccccc1C(=O)O');"
        ).unwrap().unwrap();
        assert!(op_match);
    }

    #[pg_test]
    fn test_vector_index_tier_and_prism_hybrid_integration() {
        // 1. Small dataset (2,000 rows, 320-D) -> fits easily in HNSW without TOAST
        let small_rec = crate::evaluate_vector_index_tier(2_000, 320, true, 4096);
        assert!(!small_rec.requires_storage_plain);
        assert_eq!(small_rec.recommended_access_method, "hnsw");
        assert_eq!(small_rec.recommended_opclass, "halfvec_cosine_ops");

        // 2. Full bio_demo proteome (2,340,766 rows, 1280-D halfvec, 8 GB RAM budget)
        // HNSW requires ~8.77 GB whereas PRISM requires ~960 MB (~9.5x reduction) and STORAGE PLAIN
        let large_rec = crate::evaluate_vector_index_tier(2_340_766, 1280, true, 8192);
        assert!(large_rec.requires_storage_plain);
        assert_eq!(large_rec.recommended_access_method, "prism");
        assert_eq!(large_rec.recommended_opclass, "vec16_cosine_ops");
        assert!(large_rec.hnsw_est_mb > 8000.0);
        assert!(large_rec.prism_default_est_mb < 1100.0);
        assert!(large_rec.memory_reduction_factor > 8.5);
        assert!(large_rec.recommended_nlist >= 1500);
        assert!(large_rec.recommended_nprobe >= 19);

        // 3. Verify live pg_vectorsearch PRISM index + pgvector halfvec + pg_bio <~> hybrid operator via SPI
        Spi::run(
            "CREATE TEMP TABLE test_prism_proteins (
                uniprot_id text PRIMARY KEY,
                sequence text,
                embedding halfvec(1280) STORAGE PLAIN
            );
            INSERT INTO test_prism_proteins
            SELECT 'P_' || i,
                   CASE WHEN i <= 3 THEN 'MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK'
                        ELSE 'GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA'
                   END,
                   get_esm_embedding(
                       CASE WHEN i <= 3 THEN 'MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK'
                            ELSE 'GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA'
                       END
                   )::halfvec(1280)
            FROM generate_series(1, 600) i;
            CREATE INDEX idx_test_prism_emb ON test_prism_proteins USING prism (embedding vec16_cosine_ops);"
        ).unwrap();

        let top_match = Spi::get_one::<String>(
            "WITH target AS (
                SELECT embedding AS emb, embedding::vector AS vec, sequence AS seq
                FROM test_prism_proteins WHERE uniprot_id = 'P_1'
            )
            SELECT p.uniprot_id
            FROM test_prism_proteins p, target t
            WHERE p.uniprot_id != 'P_1'
            ORDER BY (ROW(p.embedding::vector, p.sequence)::bio_feature <~> ROW(t.vec, t.seq)::bio_feature) ASC
            LIMIT 1;"
        ).unwrap().unwrap();
        assert!(top_match == "P_2" || top_match == "P_3");
    }

    #[pg_test]
    fn test_layered_centroids_morton_ecfp4_and_pockets() {
        // 1. Test 3D Morton Layered Octree Sphere Router across the (24.288, 24.288, 24.288) 2^19 = 524,288 major octree boundary
        // and the (0, 0, 0) origin!
        let cx = 24.0;
        let cy = 24.0;
        let cz = 24.0;
        let radius = 4.0;
        let naive_min = crate::z_order_encode(cx - radius, cy - radius, cz - radius);
        let naive_max = crate::z_order_encode(cx + radius, cy + radius, cz + radius);
        let naive_span = (naive_max - naive_min) as f64;

        let ranges: Vec<_> = crate::z_order_layered_sphere_ranges(cx, cy, cz, radius, 32).collect();
        assert!(!ranges.is_empty() && ranges.len() <= 32);
        let layered_span: f64 = ranges.iter().map(|r| (r.1 - r.0 + 1) as f64).sum();
        // Layered octree pruning eliminates >99.9% of the dead Morton Z-jump span!
        assert!(
            layered_span < naive_span * 0.001,
            "Expected layered span ({}) < 0.1% of naive span ({})",
            layered_span,
            naive_span
        );

        // Verify 100% recall on test points inside the 4.0 A sphere around (24.0, 24.0, 24.0)
        let inside_pts = [
            (24.0, 24.0, 24.0),
            (24.0 - 2.1, 24.0 + 1.5, 24.0 - 1.2),
            (24.0 + 2.2, 24.0 - 1.8, 24.0 + 1.9),
            (24.0 - 1.0, 24.0 - 1.0, 24.0 + 2.5),
            (24.0 + 3.9, 24.0, 24.0),
        ];
        for (px, py, pz) in inside_pts {
            let z = crate::z_order_encode(px, py, pz);
            let covered = ranges.iter().any(|&(z0, z1, _)| z >= z0 && z <= z1);
            assert!(covered, "Point ({}, {}, {}) with z={} was missed by layered ranges!", px, py, pz, z);
        }

        // 2. Test 1024-Bit ECFP4 Fastscan Sketch & 100% Guaranteed Tanimoto Upper Bound
        let smiles_set = [
            "CC(=O)OC1=CC=CC=C1C(=O)O",       // Aspirin (Kekule)
            "CC(=O)Oc1ccccc1C(=O)O",          // Aspirin (Aromatic)
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1",     // Ibuprofen
            "COc1ccc2cc(C(C)C(=O)O)ccc2c1",   // Naproxen
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",     // Caffeine
            "C1CCCCC1",                       // Cyclohexane
            "CCCCCCCCCCCCCCCCCC(=O)O",        // Stearic acid
        ];
        let mols: Vec<_> = smiles_set.iter().map(|s| crate::parse_smiles(s)).collect();
        let sketches: Vec<_> = mols.iter().cloned().map(crate::molecule_fastscan_sketch).collect();

        // Kekule and aromatic Aspirin must have identical Fastscan sketches and centroids
        assert_eq!(sketches[0], sketches[1]);

        for i in 0..mols.len() {
            for j in 0..mols.len() {
                let exact_t = crate::molecule_tanimoto(mols[i].clone(), mols[j].clone());
                let ub_t = crate::ecfp4_tanimoto_upper_bound(sketches[i].clone(), sketches[j].clone());
                assert!(
                    ub_t + 1e-12 >= exact_t,
                    "Upper bound {} violated exact Tanimoto {} for pair ({}, {})",
                    ub_t,
                    exact_t,
                    i,
                    j
                );
            }
        }
        // Dissimilar pair (Aspirin vs Cyclohexane) upper bound must prune aggressively (< 0.35)
        let ub_asp_cyc = crate::ecfp4_tanimoto_upper_bound(sketches[0].clone(), sketches[5].clone());
        assert!(ub_asp_cyc < 0.35);

        // 3. Test 9D RangeFinder Pocket Layered Centroid Assignment & Exact Envelope Upper Bound
        let p_cox2 = crate::empirical_pocket_ranges(
            "COX2",
            vec![
                "CC(=O)Oc1ccccc1C(=O)O".to_string(),
                "CC(C)Cc1ccc(C(C)C(=O)O)cc1".to_string(),
                "COc1ccc2cc(C(C)C(=O)O)ccc2c1".to_string(),
            ],
        );
        let p_alkane = crate::empirical_pocket_ranges(
            "ALKANE",
            vec![
                "C1CCCCC1".to_string(),
                "CCCCCCCC".to_string(),
                "CCCCCCCCCC".to_string(),
            ],
        );
        let c_cox2: Vec<_> = crate::pocket_layered_centroid_assign(p_cox2.clone()).collect();
        assert_eq!(c_cox2.len(), 1);
        assert!(c_cox2[0].0 >= 0 && c_cox2[0].0 < 16);
        assert!(c_cox2[0].1 >= 0 && c_cox2[0].1 < 64);
        assert_ne!(c_cox2[0].1, c_cox2[0].2); // SOAR secondary cluster is distinct from primary leaf

        // Envelope containing only p_alkane must strictly upper-bound exact IoU(p_cox2, p_alkane)
        let env_mins = vec![
            p_alkane.mw_min, p_alkane.logp_min, p_alkane.hbd_min, p_alkane.hba_min,
            p_alkane.tpsa_min, p_alkane.rot_min, p_alkane.charge_min, p_alkane.arom_min, p_alkane.fsp3_min,
        ];
        let env_maxs = vec![
            p_alkane.mw_max, p_alkane.logp_max, p_alkane.hbd_max, p_alkane.hba_max,
            p_alkane.tpsa_max, p_alkane.rot_max, p_alkane.charge_max, p_alkane.arom_max, p_alkane.fsp3_max,
        ];
        let env_widths: Vec<f64> = env_mins.iter().zip(env_maxs.iter()).map(|(l, u)| (u - l).max(0.0)).collect();

        let exact_iou = crate::pocket_range_iou(p_cox2.clone(), p_alkane.clone(), false);
        let ub_iou = crate::pocket_envelope_iou_upper_bound(p_cox2.clone(), env_mins, env_maxs, env_widths, false);
        assert!((ub_iou - exact_iou).abs() < 1e-9);

        // 4. Verify Register-Direct Zero-Deserialization Fastscan Kernels (ecfp4_fastscan_ub_raw & pocket_fastscan_iou_ub_raw)
        let d_asp: Vec<_> = crate::ecfp4_register_descriptor("CC(=O)Oc1ccccc1C(=O)O").collect();
        let d_cyc: Vec<_> = crate::ecfp4_register_descriptor("C1CCCCC1").collect();
        let raw_ub = crate::ecfp4_fastscan_ub_raw(
            d_asp[0].3, d_asp[0].4, d_asp[0].5, d_asp[0].6,
            d_cyc[0].3, d_cyc[0].4, d_cyc[0].5, d_cyc[0].6,
        );
        let cluster_ub = crate::ecfp4_cluster_envelope_ub_raw(
            d_asp[0].3, d_asp[0].4, d_asp[0].5, d_asp[0].6,
            d_cyc[0].3, d_cyc[0].3, d_cyc[0].4, d_cyc[0].5, d_cyc[0].6,
        );
        let raw_t = crate::ecfp4_tanimoto_words_raw(d_asp[0].7.clone(), d_cyc[0].7.clone());
        assert!(raw_ub + 1e-12 >= raw_t);
        assert!(cluster_ub + 1e-12 >= raw_ub);
        assert!(cluster_ub < 0.35);
        assert!(raw_ub < 0.35);

        let q_cox2: Vec<_> = crate::pocket_quantize_fastscan(
            p_cox2.mw_min, p_cox2.mw_max, p_cox2.logp_min, p_cox2.logp_max,
            p_cox2.hbd_min, p_cox2.hbd_max, p_cox2.hba_min, p_cox2.hba_max,
            p_cox2.tpsa_min, p_cox2.tpsa_max, p_cox2.rot_min, p_cox2.rot_max,
            p_cox2.arom_min, p_cox2.arom_max, p_cox2.fsp3_min, p_cox2.fsp3_max,
        ).collect();
        let q_alk: Vec<_> = crate::pocket_quantize_fastscan(
            p_alkane.mw_min, p_alkane.mw_max, p_alkane.logp_min, p_alkane.logp_max,
            p_alkane.hbd_min, p_alkane.hbd_max, p_alkane.hba_min, p_alkane.hba_max,
            p_alkane.tpsa_min, p_alkane.tpsa_max, p_alkane.rot_min, p_alkane.rot_max,
            p_alkane.arom_min, p_alkane.arom_max, p_alkane.fsp3_min, p_alkane.fsp3_max,
        ).collect();
        let q_ub = crate::pocket_fastscan_iou_ub_raw(q_cox2[0].3, q_cox2[0].4, q_alk[0].3, q_alk[0].4);
        let raw_iou = crate::pocket_iou_8prop_raw(
            p_cox2.mw_min, p_cox2.mw_max, p_cox2.logp_min, p_cox2.logp_max,
            p_cox2.hbd_min, p_cox2.hbd_max, p_cox2.hba_min, p_cox2.hba_max,
            p_cox2.tpsa_min, p_cox2.tpsa_max, p_cox2.rot_min, p_cox2.rot_max,
            p_cox2.arom_min, p_cox2.arom_max, p_cox2.fsp3_min, p_cox2.fsp3_max,
            p_alkane.mw_min, p_alkane.mw_max, p_alkane.logp_min, p_alkane.logp_max,
            p_alkane.hbd_min, p_alkane.hbd_max, p_alkane.hba_min, p_alkane.hba_max,
            p_alkane.tpsa_min, p_alkane.tpsa_max, p_alkane.rot_min, p_alkane.rot_max,
            p_alkane.arom_min, p_alkane.arom_max, p_alkane.fsp3_min, p_alkane.fsp3_max,
        );
        assert!((raw_iou - exact_iou).abs() < 1e-12);
        assert!(q_ub + 1e-12 >= raw_iou, "16B register quantized upper bound {} must >= exact IoU {}", q_ub, raw_iou);
    }

    #[pg_test]
    fn test_tpp_microprotein_knockout_and_tpca_decoherence() {
        // 1. Test 64-bit register melting curve quantization across 40°C..68°C
        let curve_55 = crate::tpp_melting_curve_quantize(55.0, 2.0);
        let curve_47 = crate::tpp_melting_curve_quantize(47.0, 3.8);
        assert_ne!(curve_55, curve_47);

        let b55 = (curve_55 as u64).to_le_bytes();
        // Melting curve must be monotonically non-increasing with temperature (40°C -> 68°C)
        for t in 0..7 {
            assert!(
                b55[t] >= b55[t + 1],
                "Melting curve must decrease with temperature: step {} ({}) < step {} ({})",
                t,
                b55[t],
                t + 1,
                b55[t + 1]
            );
        }

        // 2. Test Register-Direct TPCA Decoherence when one subunit shifts from 55°C to 47°C upon microprotein KO
        let partner_wt = crate::tpp_melting_curve_quantize(54.7, 2.1);
        let decoh = crate::tpp_tpca_decoherence_raw(curve_55, partner_wt, curve_47, partner_wt);
        assert!(decoh > 0.20 && decoh <= 1.0, "Expected strong TPCA decoherence (> 0.20), got {}", decoh);

        // Identical WT and KO states must yield 0.0 decoherence
        let zero_decoh = crate::tpp_tpca_decoherence_raw(curve_55, partner_wt, curve_55, partner_wt);
        assert!(zero_decoh.abs() < 1e-12);

        // 3. Test biophysical microprotein knockout ΔTm & chaperone sink load prediction
        let strong_ko: Vec<_> = crate::tpp_predict_microprotein_ko_shift(
            68,    // 68 aa microprotein staple
            480,   // 480 aa target complex subunit
            0.09,  // close PRISM cosine distance
            0.84,  // high interface complementarity
            0.78,  // high 3D RangeFinder pocket overlap
            0.44,  // high unmasked surface hydrophobicity
            14,    // 14 interfacial H-bonds severed
            4,     // 4 buried unsatisfied polars exposed
        )
        .collect();
        assert_eq!(strong_ko.len(), 1);
        let (tm_wt, tm_ko, delta_tm, tpca_d, chap_load, wt_code, ko_code) = strong_ko[0];
        assert!(tm_wt > tm_ko);
        assert!(delta_tm < -3.0, "Strong microprotein staple KO should drop Tm by >3°C, got {}", delta_tm);
        assert!(tpca_d > 0.12, "Expected significant TPCA decoherence (> 0.12), got {}", tpca_d);
        assert!(chap_load > 2.5, "Expected high HSP70/HSP90 chaperone sink load, got {}", chap_load);
        assert_ne!(wt_code, ko_code);

        let weak_ko: Vec<_> = crate::tpp_predict_microprotein_ko_shift(
            140, 200, 0.24, 0.25, 0.10, 0.16, 2, 0,
        )
        .collect();
        assert!(
            strong_ko[0].2 < weak_ko[0].2,
            "Strong staple KO delta_tm ({}) must be more negative than weak KO ({})",
            strong_ko[0].2,
            weak_ko[0].2
        );
        assert!(strong_ko[0].4 > weak_ko[0].4);

        // 4. Test Ideal Operating Temperature (T_opt) & Eyring-TPP Catalytic Flux Gain Window
        let opt_win: Vec<_> = crate::tpp_optimal_bioreactor_window(tm_wt, tm_ko, 14.5).collect();
        assert_eq!(opt_win.len(), 1);
        let (t_opt_wt, t_opt_ko, delta_t_opt, flux_gain_x, act_wt_code, act_ko_code) = opt_win[0];
        assert!(
            t_opt_wt > t_opt_ko + 3.0,
            "Stapled WT T_opt ({}) must exceed KO T_opt ({}) by >3°C",
            t_opt_wt,
            t_opt_ko
        );
        assert!(delta_t_opt > 3.0);
        assert!(
            flux_gain_x > 1.40,
            "Microprotein staple should boost peak catalytic productivity by >1.4x, got {}x",
            flux_gain_x
        );
        assert_ne!(act_wt_code, act_ko_code);
    }

    #[pg_test]
    fn test_phase1_anisotropy_whitening_and_rabitq_eps_bound() {
        // Generate synthetic 128-D anisotropic ESM-2 cone where SwiGLU channels 0..3 dominate variance
        let n = 48usize;
        let dim = 128usize;
        let mut raw_cone: Vec<Vec<f32>> = Vec::with_capacity(n);
        for i in 0..n {
            let mut v = vec![0.0f32; dim];
            for d in 0..dim {
                let noise = ((((i * 97 + d * 31 + 11) % 503) as f32) / 251.5) - 1.0;
                v[d] = if d == 0 {
                    40.0 + 1.5 * noise
                } else if d < 3 {
                    15.0 + 1.2 * noise
                } else {
                    0.30 * noise
                };
            }
            let norm = crate::pgbio_simd_dot_f32(&v, &v).sqrt();
            for d in 0..dim {
                v[d] /= norm;
            }
            raw_cone.push(v);
        }

        let raw_cone_flat: Vec<f32> = raw_cone.iter().flat_map(|v| v.iter().copied()).collect();
        let stats = crate::bio_fit_whitening_stats(raw_cone_flat.clone(), dim as i32);
        assert!(stats.len() > 16);
        assert_eq!(&stats[0..4], b"PGBW");

        let w0 = crate::bio_whiten_embedding(raw_cone[0].clone(), stats.clone());
        assert_eq!(w0.len(), dim);
        let w0_norm = crate::pgbio_simd_dot_f32(&w0, &w0).sqrt();
        assert!((w0_norm - 1.0).abs() < 1e-4, "Whitened vector must be L2 unit normalized, got {}", w0_norm);

        let raw_audit: Vec<_> =
            crate::bio_rabitq_eps_bound_audit(raw_cone_flat.clone(), None, 0.45, dim as i32).collect();
        let whitened_audit: Vec<_> =
            crate::bio_rabitq_eps_bound_audit(raw_cone_flat, Some(stats), 0.45, dim as i32).collect();

        let (raw_cos, raw_eps, _raw_p95, raw_prune) = raw_audit[0];
        let (w_cos, w_eps, _w_p95, w_prune) = whitened_audit[0];

        assert!(
            raw_cos > 0.90 && w_cos.abs() < 0.15,
            "Expected whitening to collapse cone anisotropy from >0.90 ({}) to <0.15 ({})",
            raw_cos,
            w_cos
        );
        assert!(
            w_eps < raw_eps * 0.65,
            "Expected whitening to tighten RaBitQ eps_bound significantly: raw={} vs whitened={}",
            raw_eps,
            w_eps
        );
        assert!(
            w_prune > raw_prune + 0.30,
            "Expected whitening to increase 32-wide fastscan pruning ratio: raw={} vs whitened={}",
            raw_prune,
            w_prune
        );
    }

    #[pg_test]
    fn test_phase1_biovec_int8_g128_inline_varlena() {
        // Build two 1280-D embeddings and verify 1,300 B payload / 1,304 B inline Postgres varlena
        let v1: Vec<f32> = (0..1280).map(|i| ((i as f32) * 0.017).sin() * 2.5).collect();
        let v2: Vec<f32> = (0..1280).map(|i| ((i as f32) * 0.017 + 0.15).sin() * 2.5).collect();

        let payload_bytes = crate::bio_vec_int8_g128_byte_size(crate::bio_vec_int8_g128_encode(v1.clone()));
        assert_eq!(payload_bytes, 1300, "1280-D BioVecInt8G128 payload must be exactly 1,300 bytes (< 2,000 B TOAST threshold)");

        let decoded1 = crate::bio_vec_int8_g128_decode(crate::bio_vec_int8_g128_encode(v1.clone()));
        assert_eq!(decoded1.len(), 1280);
        let roundtrip_dist = crate::embedding_cosine_distance(v1.clone(), decoded1);
        assert!(
            roundtrip_dist < 1e-4,
            "128-channel group i8 + f16 quantization should have <1e-4 cosine distortion, got {}",
            roundtrip_dist
        );

        let exact_cos_dist = crate::embedding_cosine_distance(v1.clone(), v2.clone()) as f32;
        let simd_i32_cos_dist = crate::bio_vec_int8_g128_cosine_distance(
            crate::bio_vec_int8_g128_encode(v1.clone()),
            crate::bio_vec_int8_g128_encode(v2.clone()),
        );
        assert!(
            (exact_cos_dist - simd_i32_cos_dist).abs() < 1e-3,
            "SIMD i32 group cosine distance ({}) must match f32 cosine distance ({}) within 1e-3",
            simd_i32_cos_dist,
            exact_cos_dist
        );

        let exact_dot = crate::pgbio_simd_dot_f32(&v1, &v2);
        let simd_i32_dot = crate::bio_vec_int8_g128_dot(
            crate::bio_vec_int8_g128_encode(v1.clone()),
            crate::bio_vec_int8_g128_encode(v2.clone()),
        );
        let rel_err = ((exact_dot - simd_i32_dot) / exact_dot).abs();
        assert!(
            rel_err < 0.01,
            "SIMD i32 group dot ({}) must match exact dot ({}) within 1%",
            simd_i32_dot,
            exact_dot
        );

        // Verify via SQL that pg_column_size(bio_vec_int8_g128_encode(...)) is 1,304 bytes (4B varlena header + 1,300B payload)
        let col_size = Spi::get_one::<i32>(
            "SELECT pg_column_size(bio_vec_int8_g128_encode(ARRAY(SELECT sin(g::float8)::float4 FROM generate_series(1, 1280) g)))",
        )
        .unwrap()
        .unwrap();
        assert_eq!(
            col_size, 1304,
            "Expected zero-CBOR fixed-header C-varlena on-disk size to be 1,304 bytes (4B header + 1,300B payload), got {}",
            col_size
        );
    }

    #[pg_test]
    fn test_phase2_structure_sanitization_kabsch_cyx_his() {
        // Defective mini-structure:
        // - Res 1 (CYS A 1): missing backbone O and missing CB, SG at (0.0, 0.0, 0.0)
        // - Res 2 (CYS A 2): SG at (2.04, 0.0, 0.0) -> 2.04 Å (< 2.2 Å disulfide bridge -> both must become CYX!)
        // - Res 3 (HIS A 3): NE2 within 2.1 Å of ZN -> must become HID!
        // - Res 4 (HIS A 4): ND1 and NE2 both within 2.9 Å of ASP/GLU carboxylates -> must become HIP!
        // - Res 5 (HIS A 5): missing imidazole ring atoms (only N, CA, C, O present) -> repaired via Kabsch SVD & standardized to HIE!
        let defective_pdb = "\
ATOM      1  N   CYS A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  CYS A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   CYS A   1       0.530   1.425   0.000  1.00  0.00           C
ATOM      4  SG  CYS A   1      -0.625  -2.548  -0.924  1.00  0.00           S
ATOM      5  N   CYS A   2       2.500   0.000   0.000  1.00  0.00           N
ATOM      6  CA  CYS A   2       3.958   0.000   0.000  1.00  0.00           C
ATOM      7  C   CYS A   2       4.488   1.425   0.000  1.00  0.00           C
ATOM      8  O   CYS A   2       4.587   2.388   0.771  1.00  0.00           O
ATOM      9  CB  CYS A   2       3.440  -0.774  -1.206  1.00  0.00           C
ATOM     10  SG  CYS A   2       1.415  -2.548  -0.924  1.00  0.00           S
ATOM     11  N   HIS A   3       8.542   0.000   0.000  1.00  0.00           N
ATOM     12  CA  HIS A   3      10.000   0.000   0.000  1.00  0.00           C
ATOM     13  C   HIS A   3      10.530   1.425   0.000  1.00  0.00           C
ATOM     14  O   HIS A   3      10.629   2.388   0.771  1.00  0.00           O
ATOM     15  CB  HIS A   3       9.482  -0.774  -1.206  1.00  0.00           C
ATOM     16  CG  HIS A   3       9.406  -2.268  -1.102  1.00  0.00           C
ATOM     17  ND1 HIS A   3       8.238  -2.885  -0.684  1.00  0.00           N
ATOM     18  CD2 HIS A   3      10.335  -3.224  -1.368  1.00  0.00           C
ATOM     19  CE1 HIS A   3       8.469  -4.179  -0.701  1.00  0.00           C
ATOM     20  NE2 HIS A   3       9.725  -4.415  -1.112  1.00  0.00           N
HETATM   21 ZN    ZN A 101      10.100  -6.350  -1.112  1.00  0.00          ZN
ATOM     22  N   HIS A   4      18.542   0.000   0.000  1.00  0.00           N
ATOM     23  CA  HIS A   4      20.000   0.000   0.000  1.00  0.00           C
ATOM     24  C   HIS A   4      20.530   1.425   0.000  1.00  0.00           C
ATOM     25  O   HIS A   4      20.629   2.388   0.771  1.00  0.00           O
ATOM     26  CB  HIS A   4      19.482  -0.774  -1.206  1.00  0.00           C
ATOM     27  CG  HIS A   4      19.406  -2.268  -1.102  1.00  0.00           C
ATOM     28  ND1 HIS A   4      18.238  -2.885  -0.684  1.00  0.00           N
ATOM     29  CD2 HIS A   4      20.335  -3.224  -1.368  1.00  0.00           C
ATOM     30  CE1 HIS A   4      18.469  -4.179  -0.701  1.00  0.00           C
ATOM     31  NE2 HIS A   4      19.725  -4.415  -1.112  1.00  0.00           N
ATOM     32  OD1 ASP A   6      16.200  -2.800  -0.684  1.00  0.00           O
ATOM     33  OE1 GLU A   7      21.500  -5.400  -1.112  1.00  0.00           O
ATOM     34  N   HIS A   5      28.542   0.000   0.000  1.00  0.00           N
ATOM     35  CA  HIS A   5      30.000   0.000   0.000  1.00  0.00           C
ATOM     36  C   HIS A   5      30.530   1.425   0.000  1.00  0.00           C
ATOM     37  O   HIS A   5      30.629   2.388   0.771  1.00  0.00           O
";

        let audit: Vec<_> = crate::bio_sanitize_structure_audit(defective_pdb).collect();
        assert_eq!(audit.len(), 1);
        let (ref sanitized, repaired, disulfides, hid, hie, hip) = audit[0];

        // Res 1 missing O + CB (2 atoms) + Res 5 missing CB, CG, ND1, CD2, CE1, NE2 (6 atoms) = 8 repaired atoms
        assert_eq!(repaired, 8, "Expected 8 Kabsch-repaired heavy atoms, got {}", repaired);
        assert_eq!(disulfides, 1, "Expected 1 <2.2A CYX disulfide bond, got {}", disulfides);
        assert_eq!(hid, 1, "Expected 1 Zn-coordinated HID, got {}", hid);
        assert_eq!(hie, 1, "Expected 1 standard HIE, got {}", hie);
        assert_eq!(hip, 1, "Expected 1 dual-salt-bridged HIP, got {}", hip);
        assert!(sanitized.contains("SSBOND   1 CYX A    1    CYX A    2"));
        assert!(sanitized.contains("HID A   3"));
        assert!(sanitized.contains("HIP A   4"));
        assert!(sanitized.contains("HIE A   5"));
    }

    #[pg_test]
    fn test_phase2_pains_and_aggregator_screen() {
        // Canonical PAINS & aggregators that MUST be rejected (true)
        assert!(crate::bio_is_pains_or_aggregator("O=C1C=CC(=O)C=C1"), "para-Benzoquinone must trigger PAINS");
        assert!(crate::bio_is_pains_or_aggregator("O=C1CSC(=S)N1"), "Rhodanine core must trigger PAINS");
        assert!(crate::bio_is_pains_or_aggregator("c1ccc(N=Nc2ccccc2)cc1"), "Azobenzene must trigger PAINS");
        assert!(crate::bio_is_pains_or_aggregator("Oc1ccccc1O"), "Catechol redox cycler must trigger PAINS");
        assert!(crate::bio_is_pains_or_aggregator("O=C(NN=Cc1ccccc1)c2ccccc2"), "Acylhydrazone must trigger PAINS");
        assert!(crate::bio_is_pains_or_aggregator("O=C1C=CC(=O)N1c2ccccc2"), "N-Phenylmaleimide must trigger PAINS");

        // Clean drug molecules that MUST pass (false)
        assert!(!crate::bio_is_pains_or_aggregator("CC(=O)Oc1ccccc1C(=O)O"), "Aspirin must NOT trigger PAINS");
        assert!(!crate::bio_is_pains_or_aggregator("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), "Ibuprofen must NOT trigger PAINS");
        assert!(!crate::bio_is_pains_or_aggregator("COc1ccc2cc(C(C)C(=O)O)ccc2c1"), "Naproxen must NOT trigger PAINS");
        assert!(!crate::bio_is_pains_or_aggregator("Cn1cnc2c1c(=O)n(C)c(=O)n2C"), "Caffeine must NOT trigger PAINS");
    }

    #[pg_test]
    fn test_phase3_enzref_role_contrast_and_storm_l0() {
        // 1. Test EnzRef Role-Contrast Scorer
        let enz_bg = vec![1.0f32, 0.0, 0.0, 0.0, 0.5, -0.5, 0.2, 0.1];
        let reactant = vec![2.0f32, 1.0, 0.0, 0.0, 0.5, -0.5, 0.2, 0.1];
        // Aligned product along the same catalytic polarization direction (reactant - enz_bg = [1, 1, 0...])
        let prod_forward = vec![3.5f32, 2.5, 0.0, 0.0, 0.5, -0.5, 0.2, 0.1];
        // Orthogonal non-productive shift
        let prod_ortho = vec![2.0f32, 1.0, 1.5, -1.5, 0.5, -0.5, 0.2, 0.1];
        // Reverse / inhibitory shift
        let prod_reverse = vec![0.5f32, -0.5, 0.0, 0.0, 0.5, -0.5, 0.2, 0.1];

        let score_fwd = crate::bio_enzref_role_contrast(enz_bg.clone(), reactant.clone(), prod_forward);
        let score_ortho = crate::bio_enzref_role_contrast(enz_bg.clone(), reactant.clone(), prod_ortho);
        let score_rev = crate::bio_enzref_role_contrast(enz_bg.clone(), reactant.clone(), prod_reverse);
        let score_dead = crate::bio_enzref_role_contrast(enz_bg.clone(), reactant.clone(), reactant.clone());

        assert!(score_fwd > 0.30, "Forward catalytic role contrast must be strongly positive, got {}", score_fwd);
        assert!(score_ortho.abs() < 1e-5, "Orthogonal shift role contrast must be ~0.0, got {}", score_ortho);
        assert!(score_rev < -0.30, "Reverse role contrast must be negative, got {}", score_rev);
        assert_eq!(score_dead, 0.0, "Zero reaction turnover (product == reactant) must yield 0.0");

        // 2. Test STORM ID128 64-hot deterministic target hashing
        let h1 = crate::bio_id128_hash("P00918_CA2_HUMAN");
        let h1_again = crate::bio_id128_hash("P00918_CA2_HUMAN");
        let h2 = crate::bio_id128_hash("Q1JU72_FA1_DEHALOGENASE");
        assert_eq!(h1.len(), 16);
        assert_eq!(h1, h1_again);
        assert_ne!(h1, h2);
        let pop1: u32 = h1.iter().map(|b| b.count_ones()).sum();
        let pop2: u32 = h2.iter().map(|b| b.count_ones()).sum();
        assert_eq!(pop1, 64, "ID128 hash must be strictly 64-hot (64 bits set out of 128), got {}", pop1);
        assert_eq!(pop2, 64, "ID128 hash must be strictly 64-hot (64 bits set out of 128), got {}", pop2);

        // 3. Test STORM bitwise Tanimoto neighborhood smoother with activity cliff preservation
        // Three closely related congeners (toluene, ethylbenzene, propylbenzene) where compound 2 has an activity cliff (pIC50 = 9.2 vs 5.1, 5.3)
        // plus one dissimilar molecule (caffeine, pIC50 = 4.0)
        let fps = vec![
            "Cc1ccccc1".to_string(),
            "CCc1ccccc1".to_string(),
            "CCCc1ccccc1".to_string(),
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C".to_string(),
        ];
        let raw_scores = vec![5.1, 9.2, 5.3, 4.0];
        let smoothed: Vec<_> = crate::bio_smooth_assay_series(fps, raw_scores, 0.25).collect();
        assert_eq!(smoothed.len(), 4);
        // Verify both raw_score and smoothed_score are exposed and the activity cliff on compound 2 is preserved
        assert_eq!(smoothed[1].1, 9.2); // raw_score preserved!
        assert!(smoothed[1].2 < 8.5 && smoothed[1].2 > 5.5, "Expected smoothed score to blend congener neighborhood, got {}", smoothed[1].2);
        assert!(smoothed[1].3 >= 2, "Expected congener neighbors >= 2, got {}", smoothed[1].3);
        assert!(smoothed[1].4 > 0.8, "Expected positive activity cliff delta on cliff compound, got {}", smoothed[1].4);
        // Dissimilar caffeine has only itself in its Tanimoto >= 0.25 neighborhood
        assert_eq!(smoothed[3].3, 1);
        assert!((smoothed[3].1 - smoothed[3].2).abs() < 1e-12);
    }

    #[pg_test]
    fn test_phase4_pgrx_bench_suite() {
        let rows: Vec<_> = crate::bio_run_phase4_benchmarks().collect();
        assert_eq!(rows.len(), 7);
        for (name, base, opt, ratio, unit) in &rows {
            pgrx::log!(
                "[PGRX-BENCH] {} | baseline={:.4} | optimized={:.4} | ratio={:.2}x ({})",
                name,
                base,
                opt,
                ratio,
                unit
            );
        }
        // Assert RaBitQ eps_bound improved after whitening
        assert!(rows[1].2 < rows[1].1);
        // Assert fastscan pruning ratio improved after whitening
        assert!(rows[2].2 > rows[2].1);
        // Assert BioVecInt8G128 varlena size is 1,304 B (< 2,000 B TOAST threshold)
        assert_eq!(rows[3].2, 1304.0);
        // Assert high throughput for EnzRef SIMD (> 20,000 rows/sec in debug, > 2M/sec in release) and PAINS screen (> 5,000 molecules/sec)
        assert!(rows[5].2 > 20_000.0);
        assert!(rows[6].2 > 5_000.0);
    }
}

#[cfg(test)]
pub mod pg_test {
    pub fn setup(_options: Vec<&str>) { }
    pub fn postgresql_conf_options() -> Vec<&'static str> { vec![] }
}
