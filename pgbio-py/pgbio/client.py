import psycopg
import json
from typing import List, Optional
from .models import Protein, Atom, SpatialBounds, PreorganizationProfile, EnergyDecompositionAudit

class PgBioClient:
    def __init__(self, connection_string: str):
        self.conn_string = connection_string

    def _get_conn(self):
        return psycopg.connect(self.conn_string)

    def find_homologues(self, query_sequence: str, limit: int = 5) -> List[Protein]:
        """
        Find proteins with similar ESM embeddings natively in the database.
        """
        query = """
            SELECT uniprot_id, name, sequence, 
                   embedding_cosine_distance(embedding, get_esm_embedding(%s)::halfvec) as distance
            FROM proteins
            ORDER BY distance ASC
            LIMIT %s;
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (query_sequence, limit))
                for row in cur.fetchall():
                    results.append(Protein(
                        uniprot_id=row[0], 
                        name=row[1], 
                        sequence=row[2], 
                        embedding_distance=row[3]
                    ))
        return results

    def align_sequence(self, target_sequence: str, threshold: int = 50, limit: int = 5) -> List[dict]:
        """
        Uses pg_bio's native rust-bio Smith-Waterman alignment to search proteins by local alignment.
        """
        query = """
            SELECT uniprot_id, name, sequence,
                   smith_waterman_score(sequence, %s, 3, -1, -2) as alignment_score
            FROM proteins
            WHERE smith_waterman_score(sequence, %s, 3, -1, -2) >= %s
            ORDER BY alignment_score DESC
            LIMIT %s;
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (target_sequence, target_sequence, threshold, limit))
                for row in cur.fetchall():
                    results.append({
                        "uniprot_id": row[0],
                        "name": row[1],
                        "sequence": row[2],
                        "alignment_score": row[3]
                    })
        return results

    def find_atoms_in_radius(self, target_x: float, target_y: float, target_z: float, radius: float) -> List[Atom]:
        """
        Uses pg_bio's massive Z-Order B-Tree spatial index to instantly find atoms within a 3D bounding box.
        """
        # We define a bounding box first for the Z-order index
        bounds = SpatialBounds(
            min_x=target_x - radius, max_x=target_x + radius,
            min_y=target_y - radius, max_y=target_y + radius,
            min_z=target_z - radius, max_z=target_z + radius
        )
        
        query = """
            WITH bounds AS (
                SELECT 
                    residue_z_index(%(min_coord)s::ResidueCoord) as min_z,
                    residue_z_index(%(max_coord)s::ResidueCoord) as max_z
            )
            SELECT atom_id, uniprot_id, coord
            FROM protein_atoms, bounds 
            WHERE z_index BETWEEN min_z AND max_z
            AND distance_angstroms(coord, %(tgt_coord)s::ResidueCoord) <= %(radius)s;
        """
        
        min_coord = json.dumps({"x": bounds.min_x, "y": bounds.min_y, "z": bounds.min_z, "name": ""})
        max_coord = json.dumps({"x": bounds.max_x, "y": bounds.max_y, "z": bounds.max_z, "name": ""})
        tgt_coord = json.dumps({"x": target_x, "y": target_y, "z": target_z, "name": ""})
        
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, {
                    "min_coord": min_coord,
                    "max_coord": max_coord,
                    "tgt_coord": tgt_coord,
                    "radius": radius
                })
                for row in cur.fetchall():
                    atom_id, uniprot_id, coord_str = row
                    coord = json.loads(coord_str) if isinstance(coord_str, str) else coord_str
                    results.append(Atom(
                        atom_id=atom_id,
                        uniprot_id=uniprot_id,
                        x=coord.get('x', 0),
                        y=coord.get('y', 0),
                        z=coord.get('z', 0),
                        name=coord.get('name', '')
                    ))
        return results

    def find_interacting_residues(self, protein_id: str, target_residue_index: int, limit: int = 5) -> List[int]:
        """
        Uses pg_bio's native Sparse Attention Maps to traverse Attention weights instantly.
        """
        query = """
            SELECT get_top_interacting_residues(attention_data, %s, %s)
            FROM protein_attention_maps
            WHERE uniprot_id = %s;
        """
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (target_residue_index, limit, protein_id))
                row = cur.fetchone()
                if row and row[0]:
                    return row[0]
        return []

    def parse_pdb_file(self, filepath: str) -> List[dict]:
        """
        Parse a PDB or mmCIF file natively using the pg_bio Rust parser.
        """
        query = """
        SELECT atom_id, residue, x, y, z 
        FROM parse_pdb_file(%s);
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (filepath,))
                for row in cur.fetchall():
                    results.append({
                        "atom_id": row[0],
                        "residue": row[1],
                        "x": row[2],
                        "y": row[3],
                        "z": row[4]
                    })
        return results

    def parse_vcf(self, filepath: str) -> List[dict]:
        """
        Parse a VCF file directly using the pg_bio native rust noodles parser.
        """
        query = """
        SELECT chrom, pos, id, ref_allele, alt_allele, qual, filter, info 
        FROM parse_vcf(%s);
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (filepath,))
                for row in cur.fetchall():
                    results.append({
                        "chrom": row[0],
                        "pos": row[1],
                        "id": row[2],
                        "ref": row[3],
                        "alt": row[4],
                        "qual": row[5],
                        "filter": row[6],
                        "info": row[7]
                    })
        return results

    def calculate_tanimoto_similarity(self, smiles_a: str, smiles_b: str) -> float:
        """
        Uses pg_bio's native Tanimoto similarity operator `%` to compare two SMILES strings.
        """
        query = "SELECT smiles_to_fingerprint(%s) %% smiles_to_fingerprint(%s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (smiles_a, smiles_b))
                row = cur.fetchone()
                if row:
                    return float(row[0])
        return 0.0

    def find_similar_compounds(self, target_smiles: str, threshold: float = 0.85, limit: int = 10) -> List[dict]:
        """
        Queries a 'library_compounds' table for compounds similar to the target SMILES.
        """
        query = """
            SELECT compound_name, smiles, 
                   (smiles_to_fingerprint(smiles) %% smiles_to_fingerprint(%s)) AS similarity
            FROM library_compounds
            WHERE (smiles_to_fingerprint(smiles) %% smiles_to_fingerprint(%s)) >= %s
            ORDER BY similarity DESC
            LIMIT %s;
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (target_smiles, target_smiles, threshold, limit))
                for row in cur.fetchall():
                    results.append({
                        "compound_name": row[0],
                        "smiles": row[1],
                        "similarity": float(row[2])
                    })
        return results

    def has_coordinated_water_shell(
        self, active_site_coords: List[dict], metal_ion: str = "ZN", max_distance: float = 2.4
    ) -> bool:
        """
        Checks whether an active-site coordinate array has an explicit metal cofactor
        and a coordinated catalytic water shell within `max_distance` Angstroms.
        """
        payload = [json.dumps(c) for c in active_site_coords]
        query = "SELECT has_coordinated_water_shell(%s::ResidueCoord[], %s, %s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (payload, metal_ion, max_distance))
                row = cur.fetchone()
                return bool(row[0]) if row else False

    def calculate_theozyme_rmsd(
        self, site_atoms: List[dict], motif_name: str = "4qhp_zn_pentad"
    ) -> float:
        """
        Computes the Kabsch-superimposed active-site RMSD (in Angstroms) against a named
        theozyme template (e.g. '4qhp_zn_pentad') or custom JSON coordinate array.
        """
        payload = [json.dumps(c) for c in site_atoms]
        query = "SELECT theozyme_rmsd(%s::ResidueCoord[], %s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (payload, motif_name))
                row = cur.fetchone()
                return float(row[0]) if row else 999.0

    def calculate_attack_angle(self, atom_a: dict, vertex_b: dict, atom_c: dict) -> float:
        """
        Computes the 3-atom bond/nucleophilic attack angle (in degrees) with `vertex_b` as vertex.
        E.g. Burgi-Dunitz angle between nucleophilic water O, carbonyl C, and carbonyl O.
        """
        query = "SELECT bio_angle(%s::ResidueCoord, %s::ResidueCoord, %s::ResidueCoord);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (json.dumps(atom_a), json.dumps(vertex_b), json.dumps(atom_c)),
                )
                row = cur.fetchone()
                return float(row[0]) if row else 0.0

    def score_interface_complementarity(
        self, pocket_embedding: List[float], peptide_embedding: List[float]
    ) -> float:
        """
        Scores two-sided enzyme-substrate interface complementarity in [0.0, 1.0].
        """
        query = "SELECT interface_complementarity(%s::real[], %s::real[]);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (pocket_embedding, peptide_embedding))
                row = cur.fetchone()
                return float(row[0]) if row else 0.0

    def check_strand_pairing(
        self, guide_strand_coords: List[dict], backbone_coords: List[dict]
    ) -> bool:
        """
        Checks antiparallel/parallel beta-strand pairing register between protease
        guide-strand coordinates and peptide substrate backbone coordinates.
        """
        guide_payload = [json.dumps(c) for c in guide_strand_coords]
        sub_payload = [json.dumps(c) for c in backbone_coords]
        query = "SELECT strand_pairing_check(%s::ResidueCoord[], %s::ResidueCoord[]);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (guide_payload, sub_payload))
                row = cur.fetchone()
                return bool(row[0]) if row else False

    def screen_metalloprotease_designs(
        self,
        motif_name: str = "4qhp_zn_pentad",
        max_rmsd: float = 0.95,
        min_preorg_fraction: float = 0.80,
    ) -> List[dict]:
        """
        Screens candidate metalloprotease designs in-database using explicit cofactor/water
        coordination, Kabsch theozyme RMSD, Burgi-Dunitz attack angles, and MD ensemble preorganization.
        """
        query = """
            SELECT design_id, series, scaffold_family,
                   theozyme_rmsd(active_site_coords, %s) AS rmsd,
                   bio_angle(water_o, substrate_c, substrate_o) AS attack_angle,
                   bio_distance(glu_oe, water_h) AS proton_transfer_dist,
                   ensemble_fraction_below(proton_transfer_dist_hist, 3.0) AS preorg_frac,
                   digested_fraction_12h
            FROM metalloprotease_designs
            WHERE has_coordinated_water_shell(active_site_coords, 'ZN', 2.4) = TRUE
              AND theozyme_rmsd(active_site_coords, %s) <= %s
              AND bio_angle(water_o, substrate_c, substrate_o) BETWEEN 100.0 AND 115.0
              AND bio_distance(glu_oe, water_h) < 3.0
              AND ensemble_fraction_below(proton_transfer_dist_hist, 3.0) >= %s
            ORDER BY rmsd ASC;
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (motif_name, motif_name, max_rmsd, min_preorg_fraction))
                for row in cur.fetchall():
                    results.append({
                        "design_id": row[0],
                        "series": row[1],
                        "scaffold_family": row[2],
                        "rmsd": float(row[3]),
                        "attack_angle": float(row[4]),
                        "proton_transfer_dist": float(row[5]),
                        "preorg_frac": float(row[6]),
                        "digested_fraction_12h": float(row[7]),
                    })
        return results

    def audit_mmpbsa_decomposition(
        self,
        model: str,
        internal: List[float],
        vdw: List[float],
        electrostatic: List[float],
        polar_solv: List[float],
        nonpolar_solv: List[float],
        reported_mm_plus_polar: List[float],
        global_ggas: float,
        global_polar: float,
        global_nonpolar_net: float,
        tolerance: float = 0.005,
    ) -> EnergyDecompositionAudit:
        """
        Runs a native O(N) Rust thermodynamic double-entry audit comparing per-residue
        MM/GBSA or MM/PBSA decomposition arrays against global endpoint energy terms.
        """
        query = """
            SELECT audit_mmpbsa_decomposition(
                %s, %s::float8[], %s::float8[], %s::float8[],
                %s::float8[], %s::float8[], %s::float8[],
                %s, %s, %s, %s
            )::text;
        """
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (
                        model,
                        internal,
                        vdw,
                        electrostatic,
                        polar_solv,
                        nonpolar_solv,
                        reported_mm_plus_polar,
                        global_ggas,
                        global_polar,
                        global_nonpolar_net,
                        tolerance,
                    ),
                )
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return EnergyDecompositionAudit(**data)

    def discover_cryptic_hotspots(
        self,
        sequence_length: int,
        residue_ids: List[int],
        mm_plus_polar: List[float],
        known_pocket_residues: List[int],
        ligand_node_id: int = 999,
        energy_cutoff: float = -0.05,
        min_stabilization: float = 0.50,
        limit: int = 5,
    ) -> List[int]:
        """
        Builds a CSR SparseAttentionMap from per-residue MM+Polar energies and traverses
        it to return stabilizing peripheral/cryptic hotspots outside `known_pocket_residues`.
        """
        query = """
            SELECT discover_cryptic_hotspots(
                build_energy_hotspot_map(%s, %s::int4[], %s::float8[], %s, %s),
                %s, %s::int4[], %s, %s
            );
        """
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (
                        sequence_length,
                        residue_ids,
                        mm_plus_polar,
                        ligand_node_id,
                        energy_cutoff,
                        ligand_node_id,
                        known_pocket_residues,
                        min_stabilization,
                        limit,
                    ),
                )
                row = cur.fetchone()
                return list(row[0]) if row and row[0] else []

    def calculate_allosteric_mechanism_index(
        self,
        sequence_length: int,
        residue_ids: List[int],
        mm_plus_polar: List[float],
        hydrophobic_residues: List[int],
        polar_ionic_residues: List[int],
        ligand_node_id: int = 999,
        energy_cutoff: float = -0.05,
    ) -> float:
        """
        Computes the normalized Polar/Ionic-to-Hydrophobic Allosteric Mechanism Index in [-1.0, +1.0].
        """
        query = """
            SELECT allosteric_mechanism_index(
                build_energy_hotspot_map(%s, %s::int4[], %s::float8[], %s, %s),
                %s, %s::int4[], %s::int4[]
            );
        """
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (
                        sequence_length,
                        residue_ids,
                        mm_plus_polar,
                        ligand_node_id,
                        energy_cutoff,
                        ligand_node_id,
                        hydrophobic_residues,
                        polar_ionic_residues,
                    ),
                )
                row = cur.fetchone()
                return float(row[0]) if row else 0.0


