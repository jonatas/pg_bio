# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
# ]
# ///

import psycopg

DB_URI = "postgresql://localhost:28818/bio_demo"

# Human Geranylgeranyl pyrophosphate synthase (GGPS1) sequence
GG_SEQ = "MEKTQETVQRILLEPYKYLLQLPGKQVRTKLSQAFNHWLKVPEDKLQIIIEVTEMLHNASLLIDDIEDNSKLRRGFPVAHSIYGIPSVINSANYVYFLGLEKVLTLDHPDAVKLFTRQLLELHQGQGLDIYWRDNYTCPTEEEYKAMVLQKTGGLFGLAVGLMQLFSDYKEDLKPLLNTLGLFFQIRDDYANLHSKEYSENKSFCEDLTEGKFSFPTIHAIWSRPESTQVQNILRQRTENIDIKKYCVHYLEDVGSFEYTRNTLKELEAKAYKQIDARGGNPELVALVKHLSKMFKEENE"

def run_discovery():
    print("=========================================================================")
    print("🌿 GERANYLGERANYL (GGPP) DISCOVERY PIPELINE")
    print("=========================================================================")
    print("Target: Geranylgeranyl Pyrophosphate Synthase (GGPS1)")
    
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # 1. Vector Homology Search
            print("\n🔍 STEP 1: Vector Homology Search")
            print("   Computing K-mer embedding for GGPS1 and finding the closest")
            print("   structural matches in our local PDB database of 495 proteins...")
            
            cur.execute("""
                WITH target AS (
                    SELECT get_esm_embedding(%s) as emb
                )
                SELECT p.uniprot_id, p.name, embedding_cosine_distance(p.embedding, t.emb) as distance
                FROM proteins p, target t
                WHERE EXISTS (SELECT 1 FROM protein_attention_maps m WHERE m.uniprot_id = p.uniprot_id)
                ORDER BY distance ASC
                LIMIT 3;
            """, (GG_SEQ,))
            
            matches = cur.fetchall()
            print("\n   Top Homologues Found:")
            for m in matches:
                print(f"   -> {m[0]} ({m[1]}) - Distance: {m[2]:.4f}")
                
            best_match = matches[0][0]
            print(f"\n   We will proceed with {best_match} for active site analysis.")
            
            # 2. Attention Traversal
            print("\n🔬 STEP 2: Allosteric Pathway Traversal")
            print("   Using Sparse Attention Maps to find which residues control the active site...")
            
            cur.execute("""
                SELECT get_top_interacting_residues(attention_data, 50, 5)
                FROM protein_attention_maps
                WHERE uniprot_id = %s;
            """, (best_match,))
            row = cur.fetchone()
            
            if row and row[0]:
                interacting = row[0]
                print(f"   [+] Residue 50 interacts heavily with: {interacting}")
            else:
                print("   [!] No attention map found for this protein. (Run setup_attention_maps.py)")
                return
                
            # 3. Spatial Binding Pocket
            print("\n⚙️ STEP 3: Z-Order Spatial Pocket Extraction")
            print("   Extracting the 3D binding pocket coordinates...")
            
            cur.execute("""
                SELECT coord::text FROM protein_atoms WHERE uniprot_id = %s LIMIT 1;
            """, (best_match,))
            atom_row = cur.fetchone()
            
            import json
            if atom_row:
                coord = json.loads(atom_row[0])
                X, Y, Z = coord['x'], coord['y'], coord['z']
                RADIUS = 8.0
                
                cur.execute(f"""
                    WITH bounds AS (
                        SELECT 
                            residue_z_index('{{"x": {X-RADIUS}, "y": {Y-RADIUS}, "z": {Z-RADIUS}, "name": ""}}'::ResidueCoord) as min_z,
                            residue_z_index('{{"x": {X+RADIUS}, "y": {Y+RADIUS}, "z": {Z+RADIUS}, "name": ""}}'::ResidueCoord) as max_z
                    )
                    SELECT atom_id, coord::text FROM protein_atoms, bounds 
                    WHERE uniprot_id = %s
                    AND z_index BETWEEN min_z AND max_z
                    AND distance_angstroms(coord, '{{"x": {X}, "y": {Y}, "z": {Z}, "name": ""}}'::ResidueCoord) <= {RADIUS};
                """, (best_match,))
                
                atoms = cur.fetchall()
                print(f"   [+] Found {len(atoms)} atoms within an {RADIUS}Å radius!")
                print(f"       Example coordinates in pocket:")
                for a in atoms[:3]:
                    print(f"       - Atom {a[0]}: {a[1]}")
                    
    print("\n=========================================================================")
    print("✅ ANALYSIS COMPLETE")
    print("=========================================================================")

if __name__ == '__main__':
    run_discovery()
