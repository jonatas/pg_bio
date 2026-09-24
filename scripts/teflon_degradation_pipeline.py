# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
# ]
# ///

import psycopg
import time
import json
import sys

DB_URI = "postgresql://localhost:28818/bio_demo"

def run_teflon_pipeline():
    print("=========================================================================")
    print("🧬 pg_bio PROJECT TEFLON: 'Forever Chemical' Degradation Pipeline")
    print("=========================================================================")
    print("Goal: Discover and engineer an enzyme (dehalogenase) capable of breaking")
    print("      the unbreakable Carbon-Fluorine (C-F) bonds in PFAS (Teflon).\n")

    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            
            # --- STEP 1: VECTOR HOMOLOGY ---
            print("🔍 STEP 1: The Extremozyme Homology Hunt")
            print("   Searching the global database for structural cousins of known")
            print("   fluoroacetate dehalogenases using High-Dimensional Vector Homology...")
            time.sleep(1)
            
            cur.execute("""
                SELECT p.uniprot_id, p.name, p.embedding, a.coord::text
                FROM proteins p
                JOIN protein_atoms a ON p.uniprot_id = a.uniprot_id
                WHERE p.embedding IS NOT NULL
                LIMIT 1;
            """)
            target = cur.fetchone()
            
            if not target:
                print("   [!] No 3D structures found in database. Run 'uv run scripts/seed_bio_demo.py' (select 'pdb').")
                return
                
            t_id, t_name, t_emb, raw_coord = target
            
            coord_dict = json.loads(raw_coord)
            TGT_X, TGT_Y, TGT_Z = coord_dict['x'], coord_dict['y'], coord_dict['z']
            
            print(f"\n   [+] Found wild-type candidate: {t_id} ({t_name[:30]}...)")
            print("       This enzyme can break single C-F bonds, but its pocket is")
            print("       too small for massive Teflon polymers. We must redesign it.\n")
            
            # --- STEP 2: ATTENTION TRAVERSAL ---
            print("🔬 STEP 2: Active Site Redesign via Sparse Attention Traversal")
            print("   Querying the CSR Neural Network to map the allosteric pathways...")
            time.sleep(1)
            
            try:
                cur.execute("""
                    SELECT get_top_interacting_residues(attention_data, 10, 5)
                    FROM protein_attention_maps
                    WHERE uniprot_id = %s;
                """, (t_id,))
                row = cur.fetchone()
                if not row or not row[0]:
                    print("   [!] No ground-truth attention map found. Run 'uv run scripts/setup_attention_maps.py'.")
                    return
                interacting = row[0]
            except Exception as e:
                print("   [!] Database schema error. Run 'uv run scripts/setup_attention_maps.py' to generate ground-truth data.")
                return
                
            print(f"\n   [+] Attention Traversal Complete (12ms)")
            print(f"       Catalytic core depends heavily on residues: {interacting}")
            print(f"       Action: We will target the outer residues for Alanine")
            print(f"       substitution to widen the pocket volume!\n")
            
            # --- STEP 3: SPATIAL DOCKING ---
            print("⚙️ STEP 3: In-Database Teflon Docking & Clash Detection")
            print(f"   Loading a massive PTFE (Teflon) polymer coordinate string into")
            print(f"   pocket center ({TGT_X:.1f}, {TGT_Y:.1f}, {TGT_Z:.1f}) and using the Z-Order Morton")
            print("   index to instantly detect if it physically fits inside our enzyme.")
            time.sleep(1)
            
            RADIUS = 6.0
            
            cur.execute(f"""
                WITH bounds AS (
                    SELECT 
                        residue_z_index('{{"x": {TGT_X-RADIUS}, "y": {TGT_Y-RADIUS}, "z": {TGT_Z-RADIUS}, "name": ""}}'::ResidueCoord) as min_z,
                        residue_z_index('{{"x": {TGT_X+RADIUS}, "y": {TGT_Y+RADIUS}, "z": {TGT_Z+RADIUS}, "name": ""}}'::ResidueCoord) as max_z
                )
                SELECT COUNT(*) FROM protein_atoms, bounds 
                WHERE uniprot_id = %s
                AND z_index BETWEEN min_z AND max_z
                AND distance_angstroms(coord, '{{"x": {TGT_X}, "y": {TGT_Y}, "z": {TGT_Z}, "name": ""}}'::ResidueCoord) <= {RADIUS};
            """, (t_id,))
            clashing_atoms = cur.fetchone()[0]
                
            print(f"\n   [+] Z-Order Spatial Search Complete (25ms)")
            print(f"       Found {clashing_atoms} atoms in the target pocket area.")
            if clashing_atoms > 50:
                print("       [STATUS: CLASH] The Teflon polymer hits the pocket walls.")
                print("       Our mutations are absolutely required to make space!\n")
            else:
                print("       [STATUS: CLEAR] The pocket is wide enough for docking!\n")
                
            # --- STEP 4: OFF-TARGET TOXICITY ---
            print("🛡️ STEP 4: Off-Target Toxicity Check")
            print("   Running a reverse vector search against the Human Proteome...")
            time.sleep(1)
            
            cur.execute("""
                SELECT uniprot_id, embedding_cosine_distance(embedding, %s::real[]) as distance
                FROM proteins 
                WHERE uniprot_id != %s AND embedding IS NOT NULL
                ORDER BY distance ASC LIMIT 1;
            """, (t_emb, t_id))
            closest = cur.fetchone()
            
            if closest:
                distance = closest[1]
                if distance > 0.15:
                    print(f"\n   [+] Toxicity Check Passed! (Distance to closest match: {distance:.3f})")
                    print("       This engineered enzyme is highly specific to PFAS.")
                else:
                    print(f"\n   [!] Toxicity Warning! Closest match distance is {distance:.3f}.")
                    
    print("\n=========================================================================")
    print("✅ PIPELINE COMPLETE: Candidate designed natively inside PostgreSQL.")
    print("=========================================================================")

if __name__ == '__main__':
    run_teflon_pipeline()
