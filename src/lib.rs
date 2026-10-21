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
// 5. SMALL MOLECULE CHEMISTRY (Cheminformatics)
// =====================================================================

/// Represents a small molecule parsed from a SMILES string.
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct Molecule {
    pub smiles: String,
    pub is_valid: bool,
}

/// Parses and validates a SMILES chemical string.
#[pg_extern(immutable, parallel_safe)]
pub fn parse_smiles(smiles_str: &str) -> Molecule {
    let mut builder = purr::graph::Builder::new();
    let is_valid = purr::read::read(smiles_str, &mut builder, None).is_ok();
    Molecule {
        smiles: smiles_str.to_string(),
        is_valid,
    }
}

/// Computes a very basic boolean "substructure match".
/// Note: In a production extension, this would do a full molecular graph isomorphism check.
#[pg_extern(immutable, parallel_safe)]
pub fn smiles_contains(target: Molecule, substructure: &str) -> bool {
    if !target.is_valid {
        return false;
    }
    target.smiles.contains(substructure)
}

/// Native Morgan (ECFP) Fingerprint generation in Rust
#[pg_extern(immutable, parallel_safe)]
pub fn smiles_to_fingerprint(smiles: &str) -> Vec<f32> {
    let mut builder = purr::graph::Builder::new();
    if purr::read::read(smiles, &mut builder, None).is_err() {
        return vec![0.0; 1024]; // Return empty on parse fail
    }
    let atoms = match builder.build() {
        Ok(a) => a,
        Err(_) => return vec![0.0; 1024]
    };
    let mut fp = vec![0.0; 1024];

    let mut hashes: Vec<u64> = atoms.iter().enumerate().map(|(i, a)| {
        let mut h = (i as u64).wrapping_add(1) * 13;
        let kind_str = format!("{:?}", a.kind);
        for byte in kind_str.bytes() {
            h = h.wrapping_mul(31).wrapping_add(byte as u64);
        }
        h = h.wrapping_mul(31).wrapping_add(a.bonds.len() as u64);
        // Set bit for radius 0
        fp[(h as usize) % 1024] = 1.0;
        h
    }).collect();

    for _radius in 1..=2 {
        let mut new_hashes = hashes.clone();
        for (i, a) in atoms.iter().enumerate() {
            let mut neighbor_hashes: Vec<u64> = a.bonds.iter().map(|b| {
                let neighbor_id = b.tid;
                let bond_hash = match b.kind {
                    purr::feature::BondKind::Single => 1,
                    purr::feature::BondKind::Double => 2,
                    purr::feature::BondKind::Triple => 3,
                    purr::feature::BondKind::Quadruple => 4,
                    purr::feature::BondKind::Aromatic => 5,
                    purr::feature::BondKind::Up => 6,
                    purr::feature::BondKind::Down => 7,
                    purr::feature::BondKind::Elided => 8,
                };
                hashes[neighbor_id].wrapping_mul(17).wrapping_add(bond_hash)
            }).collect();
            
            neighbor_hashes.sort_unstable();
            
            let mut h = hashes[i];
            for nh in neighbor_hashes {
                h = h.wrapping_mul(31).wrapping_add(nh);
            }
            new_hashes[i] = h;
            // Set bit
            fp[(h as usize) % 1024] = 1.0;
        }
        hashes = new_hashes;
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
-- 5. Create Tanimoto similarity operator for float arrays
CREATE OPERATOR % (
    LEFTARG = real[],
    RIGHTARG = real[],
    PROCEDURE = tanimoto_similarity,
    COMMUTATOR = %
);
"#,
    name = "pg_bio_tanimoto_operator",
    requires = [tanimoto_similarity]
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
    use noodles::bgzf;
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
    let (pdb, _errors) = pdbtbx::open(filepath).unwrap_or_else(|e| panic!("Failed to open PDB file {}: {:?}", filepath, e));
    
    let mut records = Vec::new();
    for residue in pdb.residues() {
        let res_name = residue.name().map(|n| n.to_string()).unwrap_or_default();
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
    
    pgrx::iter::TableIterator::new(records.into_iter())
}


// =====================================================================
// NEW: STRUCTURAL BIOLOGY (PDB Parsing & RMSD)
// =====================================================================

#[pg_extern]
pub fn bio_parse_pdb(pdb_text: &str) -> TableIterator<'static, (name!(atom_name, String), name!(res_name, String), name!(res_id, i32), name!(x, f64), name!(y, f64), name!(z, f64))> {
    use std::io::{BufReader, Cursor};
    let reader = BufReader::new(Cursor::new(pdb_text.as_bytes()));
    let (pdb, _errors) = pdbtbx::ReadOptions::new().set_format(pdbtbx::Format::Pdb).read_raw(reader).unwrap();
    
    let mut results = Vec::new();
    for model in pdb.models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                for atom in residue.atoms() {
                    results.push((
                        atom.name().to_string(),
                        residue.name().unwrap_or("").to_string(),
                        residue.serial_number() as i32,
                        atom.pos().0,
                        atom.pos().1,
                        atom.pos().2
                    ));
                }
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
}

#[cfg(test)]
pub mod pg_test {
    pub fn setup(_options: Vec<&str>) { }
    pub fn postgresql_conf_options() -> Vec<&'static str> { vec![] }
}
