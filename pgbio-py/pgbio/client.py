import psycopg
import json
from typing import List, Optional
from .models import (
    Protein,
    Atom,
    SpatialBounds,
    PreorganizationProfile,
    EnergyDecompositionAudit,
    ExpressionCassette,
    GlobularFoldingProfile,
    RasaProfile,
    HBondAudit,
    ChaperoneEvaluation,
    ChaperoneCoreExampleAudit,
    LigandProperties,
    RangeFinderPocket,
    RangeFinderEvaluation,
)

class PgBioClient:
    def __init__(self, connection_string: str):
        self.conn_string = connection_string

    def _get_conn(self):
        return psycopg.connect(self.conn_string)

    def find_homologues(self, query_sequence: str, limit: int = 5) -> List[Protein]:
        """
        Find proteins with similar ESM embeddings natively in the database using the HNSW `<=>` operator.
        """
        query = """
            SELECT uniprot_id, name, sequence, 
                   (embedding <=> get_esm_embedding(%s)::halfvec) as distance
            FROM proteins
            WHERE embedding IS NOT NULL
            ORDER BY embedding <=> get_esm_embedding(%s)::halfvec ASC
            LIMIT %s;
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (query_sequence, query_sequence, limit))
                for row in cur.fetchall():
                    results.append(Protein(
                        uniprot_id=row[0], 
                        name=row[1], 
                        sequence=row[2], 
                        embedding_distance=float(row[3]) if row[3] is not None else None
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
        Automatically applies Hückel 4n+2 aromaticity perception, salt stripping, and tautomer normalization.
        """
        query = "SELECT parse_smiles(%s) %% parse_smiles(%s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (smiles_a, smiles_b))
                row = cur.fetchone()
                if row:
                    return float(row[0])
        return 0.0

    def get_molecule_wl_hash(self, smiles: str) -> int:
        """
        Computes the canonical 64-bit Weisfeiler-Lehman (WL) graph hash after Hückel 4n+2
        aromaticity perception, salt stripping, and prototropic tautomer normalization.
        """
        query = "SELECT molecule_wl_hash(%s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (smiles,))
                row = cur.fetchone()
                return int(row[0]) if row and row[0] is not None else 0

    def molecules_are_identical(self, smiles_a: str, smiles_b: str) -> bool:
        """
        Checks whether two SMILES strings (e.g. Kekulé vs. aromatic `C1=CC=CC=C1` vs `c1ccccc1`,
        or salted `CC(=O)[O-].[Na+]` vs `CC(=O)O`) compile to the exact same normalized molecular graph.
        """
        query = "SELECT molecules_are_identical(%s, %s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (smiles_a, smiles_b))
                row = cur.fetchone()
                return bool(row[0]) if row else False

    def contains_substructure(self, target_smiles: str, query_smiles: str) -> bool:
        """
        Uses pg_bio's `@>` operator (1024-bit ECFP4 screening + exact VF2 subgraph isomorphism recheck)
        with Hückel 4n+2 aromaticity normalization to test if `target_smiles` contains `query_smiles`.
        """
        query = "SELECT parse_smiles(%s) @> parse_smiles(%s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (target_smiles, query_smiles))
                row = cur.fetchone()
                return bool(row[0]) if row else False

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

    def reverse_translate(self, protein_seq: str, host_organism: str = "e_coli") -> str:
        """
        Reverse-translates an amino acid sequence into a host-optimized DNA coding sequence (CDS) natively in SQL.
        """
        query = "SELECT bio_reverse_translate(%s, %s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (protein_seq, host_organism))
                row = cur.fetchone()
                return str(row[0]) if row else ""

    def calculate_codon_adaptation_index(self, dna_seq: str, host_organism: str = "e_coli") -> float:
        """
        Computes the Codon Adaptation Index (CAI in [0.0, 1.0]) of a DNA sequence for a target expression host.
        """
        query = "SELECT bio_codon_adaptation_index(%s, %s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (dna_seq, host_organism))
                row = cur.fetchone()
                return float(row[0]) if row else 0.0

    def design_expression_cassette(
        self, protein_seq: str, host_organism: str = "e_coli", promoter_type: str = "T7"
    ) -> ExpressionCassette:
        """
        Designs a complete wet-lab expression cassette (Promoter + RBS/Kozak + Host-Optimized CDS + Stop + Terminator)
        and audits GC content, CAI score, and BsaI Golden Gate restriction site safety natively in PostgreSQL.
        """
        query = "SELECT bio_design_expression_cassette(%s, %s, %s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (protein_seq, host_organism, promoter_type))
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return ExpressionCassette(**data)

    def evaluate_globular_compactness(self, pdb_text: str) -> GlobularFoldingProfile:
        """
        Evaluates the 3D Globular Compactness Gate (mean pLDDT, Radius of Gyration Rg, Flory compactness ratio,
        and gyration tensor principal-axis elongation) directly from PDB text inside PostgreSQL.
        """
        query = "SELECT evaluate_globular_compactness(%s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (pdb_text,))
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return GlobularFoldingProfile(**data)

    def simulate_defluorination_kinetics(
        self,
        kcat_per_hr: float,
        km_u_m: float,
        substrate_u_m: float,
        hydrophobin_anchor_boost: float,
        calcium_m_m: float,
        time_hours: float,
    ) -> dict:
        """
        Runs the in-database Michaelis-Menten + Hydrophobin surface adsorption + Fluorapatite Ca5(PO4)3F
        biomineralization Digital Twin for a given time point.
        """
        query = """
            SELECT time_hr, fluoride_ppm, fluorapatite_mg_per_l, fractional_conversion
            FROM simulate_defluorination_kinetics(%s, %s, %s, %s, %s, %s);
        """
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (
                        kcat_per_hr,
                        km_u_m,
                        substrate_u_m,
                        hydrophobin_anchor_boost,
                        calcium_m_m,
                        time_hours,
                    ),
                )
                row = cur.fetchone()
                if not row:
                    return {}
                return {
                    "time_hr": float(row[0]),
                    "fluoride_ppm": float(row[1]),
                    "fluorapatite_mg_per_l": float(row[2]),
                    "fractional_conversion": float(row[3]),
                }

    def calculate_rasa(self, pdb_text: str, probe_radius: float = 1.4) -> RasaProfile:
        """
        Computes Shrake-Rupley SASA and Tien (2013) Relative Accessible Surface Area (RASA)
        natively in Rust inside PostgreSQL (Chaperone `rasa` block).
        """
        query = "SELECT bio_rasa(%s, %s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (pdb_text, probe_radius))
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return RasaProfile(**data)

    def audit_hydrogen_bonds(
        self, pdb_text: str, binder_chain: str = "A", target_chain: str = "B"
    ) -> HBondAudit:
        """
        Audits the 3D Hydrogen-Bond network (backbone, sidechain, cross-interface) and counts
        Buried Unsatisfied Polar Atoms (BUNS) natively in PostgreSQL (Chaperone `hbond` block).
        """
        query = "SELECT audit_hydrogen_bonds(%s, %s, %s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (pdb_text, binder_chain, target_chain))
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return HBondAudit(**data)

    def evaluate_chaperone_design(
        self,
        candidate_pdb: str,
        reference_backbone_pdb: str,
        binder_chain: str = "A",
        target_chain: str = "B",
    ) -> ChaperoneEvaluation:
        """
        Executes Chaperone's post-oracle evaluation suite (`rmsd` + `rasa` + `hbond` + `globular`)
        in a single fused Rust pass inside PostgreSQL.
        """
        query = "SELECT evaluate_chaperone_design(%s, %s, %s, %s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (candidate_pdb, reference_backbone_pdb, binder_chain, target_chain),
                )
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return ChaperoneEvaluation(**data)

    def compute_chaperone_lineage_hash(
        self, stage_name: str, parent_id: str, sequence: str, cycle_num: int = 0
    ) -> str:
        """
        Generates a deterministic content-addressable lineage node ID for Chaperone
        PipelineDataHandler rows across multi-stage workflows and cycles.
        """
        query = "SELECT chaperone_lineage_hash(%s, %s, %s, %s);"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (stage_name, parent_id, sequence, cycle_num))
                row = cur.fetchone()
                return str(row[0]) if row else ""

    def list_chaperone_core_examples(self) -> List[dict]:
        """
        Returns the 8 official RosettaCommons Chaperone Core Examples
        (`nanoparticle-design`, `peptide-binder`, `cortisol-binder`, `dna-binder`,
        `enzyme-glycosidase`, `enzyme-protease`, `ribozyme-hammerhead`, `multi-mpnn`)
        with their ready-to-evaluate 3D reference and candidate PDB structures.
        """
        query = """
            SELECT example_name, category, pipeline_blocks, target_modality,
                   binder_chain, target_chain, ligand_smiles, sequence,
                   reference_pdb, candidate_pdb
            FROM chaperone_core_examples();
        """
        results = []
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query)
                for row in cur.fetchall():
                    results.append(
                        {
                            "example_name": row[0],
                            "category": row[1],
                            "pipeline_blocks": row[2],
                            "target_modality": row[3],
                            "binder_chain": row[4],
                            "target_chain": row[5],
                            "ligand_smiles": row[6],
                            "sequence": row[7],
                            "reference_pdb": row[8],
                            "candidate_pdb": row[9],
                        }
                    )
        return results

    def evaluate_chaperone_core_example(
        self,
        example_name: str,
        candidate_pdb: str,
        reference_pdb: str,
        ligand_smiles: str = "",
    ) -> ChaperoneCoreExampleAudit:
        """
        Evaluates any of the 8 official RosettaCommons Chaperone Core Examples natively in PostgreSQL.
        """
        query = "SELECT evaluate_chaperone_core_example(%s, %s, %s, %s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (example_name, candidate_pdb, reference_pdb, ligand_smiles),
                )
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return ChaperoneCoreExampleAudit(**data)

    def compute_ligand_properties(self, smiles: str) -> LigandProperties:
        """
        Computes the 9 RangeFinder physicochemical properties (MW, logP, HBD, HBA, TPSA,
        rotatable bonds, formal charge, aromatic rings, Fsp3) directly from a SMILES string in Rust.
        """
        query = "SELECT compute_ligand_properties(%s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (smiles,))
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return LigandProperties(**data)

    def predict_pocket_ranges(
        self,
        pocket_id: str,
        pdb_text: str,
        center_x: float,
        center_y: float,
        center_z: float,
        radius_angstroms: float = 10.0,
    ) -> RangeFinderPocket:
        """
        Predicts the 9 RangeFinder ligand property intervals [mu - 2|sigma|, mu + 2|sigma|]
        and aromatic presence probability directly from a 10 Å protein pocket in PDB format.
        """
        query = "SELECT predict_pocket_ranges(%s, %s, %s, %s, %s, %s)::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (pocket_id, pdb_text, center_x, center_y, center_z, radius_angstroms),
                )
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return RangeFinderPocket(**data)

    def empirical_pocket_ranges(
        self, pocket_id: str, smiles_list: List[str]
    ) -> RangeFinderPocket:
        """
        Computes empirical RangeFinder target intervals [mu - 2*sigma, mu + 2*sigma]
        over a list of known binding ligands for a pocket.
        """
        query = "SELECT empirical_pocket_ranges(%s, %s::text[])::text;"
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (pocket_id, smiles_list))
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return RangeFinderPocket(**data)

    def evaluate_rangefinder_screening(
        self,
        pocket_id: str,
        true_binder_smiles: List[str],
        background_library_smiles: List[str],
        include_formal_charge: bool = False,
    ) -> RangeFinderEvaluation:
        """
        Evaluates a RangeFinder pocket's Range-IoU, GIoU loss, Library Reduction (1/r_p),
        Recall, and Reduction-Recall F-Score (RRF) inside PostgreSQL.
        """
        query = """
            SELECT evaluate_rangefinder_screening(
                empirical_pocket_ranges(%s, %s::text[]),
                %s::text[],
                %s::text[],
                %s
            )::text;
        """
        with self._get_conn() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    query,
                    (
                        pocket_id,
                        true_binder_smiles,
                        true_binder_smiles,
                        background_library_smiles,
                        include_formal_charge,
                    ),
                )
                row = cur.fetchone()
                data = json.loads(row[0]) if isinstance(row[0], str) else row[0]
                return RangeFinderEvaluation(**data)





