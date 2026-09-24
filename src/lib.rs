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
// 2. HIGH-DIMENSIONAL VECTOR EMBEDDINGS (ESM AI MODELS)
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

/// Represents a single non-zero interaction in an AI Attention Map.
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
// TESTS
// =====================================================================

#[cfg(any(test, feature = "pg_test"))]
#[pg_schema]
mod tests {
    use pgrx::prelude::*;
    use crate::{ResidueCoord, z_order_encode, residue_z_index, embedding_cosine_distance};
    use crate::{AttentionEntry, SparseAttentionMap, get_top_interacting_residues};

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
        assert!(crate::prosite_match("ASAAAR", pattern) == true); 
        assert!(crate::prosite_match("ATAAAR", pattern) == true); 
        assert!(crate::prosite_match("AGGGGD", pattern) == false); 
        
        let pattern2 = "{ST}-x-[RK]";
        assert!(crate::prosite_match("AAAR", pattern2) == true); 
        assert!(crate::prosite_match("ASAR", pattern2) == false); 
    }
}

#[cfg(test)]
pub mod pg_test {
    pub fn setup(_options: Vec<&str>) { }
    pub fn postgresql_conf_options() -> Vec<&'static str> { vec![] }
}
