import psycopg

DB_URI = "postgresql://localhost:28818/bio_demo"

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            print("--- VALIDATION REPORT ---")
            
            # 1. Total counts
            cur.execute("SELECT COUNT(*) FROM proteins;")
            total_proteins = cur.fetchone()[0]
            print(f"Total proteins: {total_proteins}")
            
            cur.execute("SELECT COUNT(*) FROM protein_atoms;")
            total_atoms = cur.fetchone()[0]
            print(f"Total atoms: {total_atoms}")
            
            # 2. Proteins without atoms
            cur.execute("""
                SELECT COUNT(*) FROM proteins p 
                LEFT JOIN protein_atoms a ON p.uniprot_id = a.uniprot_id 
                WHERE a.uniprot_id IS NULL;
            """)
            empty_proteins = cur.fetchone()[0]
            print(f"Proteins with 0 atoms: {empty_proteins}")
            
            # 3. Proteins with missing sequences or embeddings
            cur.execute("SELECT COUNT(*) FROM proteins WHERE sequence IS NULL OR sequence = '';")
            missing_seqs = cur.fetchone()[0]
            print(f"Proteins missing sequence: {missing_seqs}")
            
            cur.execute("SELECT COUNT(*) FROM proteins WHERE embedding IS NULL;")
            missing_emb = cur.fetchone()[0]
            print(f"Proteins missing embedding: {missing_emb}")
            
            # 4. Atoms with missing coordinates or Z-index
            cur.execute("SELECT COUNT(*) FROM protein_atoms WHERE coord IS NULL;")
            missing_coord = cur.fetchone()[0]
            print(f"Atoms missing coordinates: {missing_coord}")
            
            cur.execute("SELECT COUNT(*) FROM protein_atoms WHERE z_index IS NULL;")
            missing_z = cur.fetchone()[0]
            print(f"Atoms missing Z-index: {missing_z}")

            # 5. Embedding lengths
            cur.execute("SELECT MIN(array_length(embedding, 1)), MAX(array_length(embedding, 1)) FROM proteins WHERE embedding IS NOT NULL;")
            emb_min, emb_max = cur.fetchone()
            print(f"Embedding length min/max: {emb_min} / {emb_max}")

except Exception as e:
    print(f"Error connecting to database: {e}")
