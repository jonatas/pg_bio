import psycopg
import json
from typing import List, Optional
from .models import Protein, Atom, SpatialBounds

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
