from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Optional
from pgbio import PgBioClient
import os

app = FastAPI(
    title="pg_bio API",
    description="Backend engine for pg_bio computational biology database.",
    version="1.0.0"
)

# Connect to pg_bio (fallback to default local connection if not set)
DB_URI = os.getenv("DATABASE_URL", "postgresql://localhost:28818/bio_demo")
client = PgBioClient(DB_URI)

class SequenceQuery(BaseModel):
    sequence: str
    limit: int = 5

class AtomResponse(BaseModel):
    atom_id: int
    uniprot_id: str
    x: float
    y: float
    z: float
    name: str

class ResearchWorkflowRequest(BaseModel):
    sequence: str

@app.get("/")
def health_check():
    return {"status": "ok", "message": "pg_bio backend engine is running."}

@app.post("/api/homologues")
def find_homologues(query: SequenceQuery):
    """Find functionally/structurally similar proteins via K-mer vector embeddings."""
    try:
        proteins = client.find_homologues(query.sequence, limit=query.limit)
        return {"results": [{"uniprot_id": p.uniprot_id, "name": p.name, "distance": p.embedding_distance} for p in proteins]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/proteins/{protein_id}/attention/{residue_idx}")
def get_attention_network(protein_id: str, residue_idx: int, limit: int = 5):
    """Traverse the Sparse Attention Map to find physically/evolutionarily interacting amino acids."""
    try:
        interactions = client.find_interacting_residues(protein_id, residue_idx, limit)
        if not interactions:
            return {"protein_id": protein_id, "target_residue": residue_idx, "interactions": [], "message": "No attention data found."}
        return {"protein_id": protein_id, "target_residue": residue_idx, "interacting_residues": interactions}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/spatial/radius")
def find_atoms_in_radius(x: float, y: float, z: float, radius: float = 5.0):
    """Use Z-Order spatial indexing to instantly find atoms within a 3D bounding box across the entire database."""
    try:
        atoms = client.find_atoms_in_radius(x, y, z, radius)
        return {"center": {"x": x, "y": y, "z": z}, "radius": radius, "atoms_found": len(atoms), "atoms": atoms}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/research/workflow")
def full_research_workflow(query: ResearchWorkflowRequest):
    """
    Demonstrates a full research pipeline:
    1. Find closest structural homologue to the input sequence.
    2. Check the attention map for its most critical node (e.g. residue 50).
    3. Find the 3D binding pocket (atoms within 5A) around that critical node.
    """
    # Step 1: Homology
    proteins = client.find_homologues(query.sequence, limit=1)
    if not proteins:
        raise HTTPException(status_code=404, detail="No homologues found in database.")
    
    top_protein = proteins[0]
    
    # Step 2: Attention Map (Mocking target residue 50 for the sake of the demo, 
    # since we seeded attention maps with random residues)
    target_residue = 50
    interacting_residues = client.find_interacting_residues(top_protein.uniprot_id, target_residue, limit=5)
    
    # Step 3: Spatial binding pocket
    # Let's find the coordinate of the target residue first to use as our center
    center_coord = None
    with client._get_conn() as conn:
        with conn.cursor() as cur:
            # We loosely map residue 50 to an atom (assuming 1 atom per residue in our naive mock, or just taking the 50th atom)
            cur.execute("""
                SELECT coord FROM protein_atoms 
                WHERE uniprot_id = %s 
                ORDER BY atom_id ASC 
                OFFSET %s LIMIT 1
            """, (top_protein.uniprot_id, target_residue - 1))
            row = cur.fetchone()
            if row:
                import json
                coord_str = row[0]
                center_coord = json.loads(coord_str) if isinstance(coord_str, str) else coord_str

    pocket_atoms = []
    if center_coord:
        pocket_atoms = client.find_atoms_in_radius(
            center_coord.get('x', 0), 
            center_coord.get('y', 0), 
            center_coord.get('z', 0), 
            radius=6.0
        )
        
    return {
        "step_1_homology": {
            "query_sequence": query.sequence,
            "best_match": top_protein.name,
            "match_id": top_protein.uniprot_id,
            "distance": top_protein.embedding_distance
        },
        "step_2_attention": {
            "target_residue": target_residue,
            "strongly_interacting_residues": interacting_residues,
            "insight": "These residues likely form the functional active site based on Evolutionary weights."
        },
        "step_3_spatial_pocket": {
            "center_coordinate": center_coord,
            "atoms_in_6A_radius": len(pocket_atoms),
            "insight": "Instantly retrieved the 3D binding pocket using Z-Order indexing."
        }
    }
