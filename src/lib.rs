use pgrx::prelude::*;
use serde::{Deserialize, Serialize};

pgrx::pg_module_magic!();

/// A custom PostgreSQL type representing a 3D physical coordinate
/// of an amino acid residue in space (e.g., from a PDB file).
#[derive(PostgresType, Serialize, Deserialize, Debug, Clone)]
pub struct ResidueCoord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub name: String, // e.g., "TRP", "GLY"
}

/// Helper to interleave bits for 3D Z-Order (Morton Coding).
/// Takes a 21-bit integer and spaces its bits out by 3.
#[inline]
fn split_by_3(x: u64) -> u64 {
    let mut res = 0;
    for i in 0..21 {
        res |= ((x >> i) & 1) << (3 * i);
    }
    res
}

/// Computes the 1D Z-Order curve index (Morton Code) for a 3D coordinate.
/// Maps floats within a [-500.0, 500.0] Ångstrom bounding box into a 63-bit integer.
#[pg_extern(immutable, parallel_safe)]
pub fn z_order_encode(x: f64, y: f64, z: f64) -> i64 {
    // Normalize coordinates to [0, 2^21 - 1].
    // Shift by +500 to make them positive, scale by 1000 for sub-Ångstrom precision (0.001 Å).
    let map = |v: f64| -> u64 {
        let v = v.clamp(-500.0, 500.0) + 500.0;
        let scaled = (v * 1000.0) as u64; // Max value ~ 1,000,000, which fits in 21 bits (2,097,151)
        scaled & 0x1FFFFF 
    };

    let xx = split_by_3(map(x));
    let yy = split_by_3(map(y));
    let zz = split_by_3(map(z));

    // Interleave X, Y, and Z bits: (Z_n, Y_n, X_n, ..., Z_0, Y_0, X_0)
    (xx | (yy << 1) | (zz << 2)) as i64
}

/// A Postgres function to extract the Z-order index directly from our custom type.
/// By making it IMMUTABLE, Postgres can use this function to build B-Tree indexes!
#[pg_extern(immutable, parallel_safe)]
pub fn residue_z_index(coord: ResidueCoord) -> i64 {
    z_order_encode(coord.x, coord.y, coord.z)
}

/// Computes the exact 3D Euclidean distance (in Ångstroms).
#[pg_extern(immutable, parallel_safe)]
pub fn distance_angstroms(a: ResidueCoord, b: ResidueCoord) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// A mocked function that demonstrates fetching an ESM embedding vector.
#[pg_extern]
pub fn get_esm_embedding(_sequence: &str) -> Vec<f32> {
    vec![0.01; 1280] 
}

#[cfg(any(test, feature = "pg_test"))]
#[pg_schema]
mod tests {
    use pgrx::prelude::*;
    use crate::{ResidueCoord, z_order_encode, residue_z_index};

    #[pg_test]
    fn test_distance() {
        let r1 = ResidueCoord { x: 0.0, y: 0.0, z: 0.0, name: "ALA".to_string() };
        let r2 = ResidueCoord { x: 3.0, y: 4.0, z: 0.0, name: "GLY".to_string() };
        let d = crate::distance_angstroms(r1, r2);
        assert!((d - 5.0).abs() < 1e-6);
    }

    #[pg_test]
    fn test_z_order_encoding() {
        // Two atoms very close to each other should have similar Z-order prefixes
        let r1 = ResidueCoord { x: 10.0, y: 10.0, z: 10.0, name: "TRP".to_string() };
        let r2 = ResidueCoord { x: 10.1, y: 10.1, z: 10.1, name: "PHE".to_string() };
        
        let z1 = residue_z_index(r1);
        let z2 = residue_z_index(r2);
        
        // Z-values should be close since the coordinates are spatially close
        assert!(z1 > 0);
        assert!(z2 > 0);
        
        // Ensure a distant atom has a drastically different Z-order
        let r3 = ResidueCoord { x: 300.0, y: 10.0, z: 10.0, name: "LEU".to_string() };
        let z3 = residue_z_index(r3);
        assert!((z1 - z2).abs() < (z1 - z3).abs());
    }
}

#[cfg(test)]
pub mod pg_test {
    pub fn setup(_options: Vec<&str>) { }
    pub fn postgresql_conf_options() -> Vec<&'static str> { vec![] }
}
