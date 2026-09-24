import psycopg

DB_URI = "postgresql://localhost:28818/bio_demo"

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            print("--- 1. Molecular Weight ---")
            cur.execute("""
                SELECT uniprot_id, molecular_weight(sequence) as mass
                FROM proteins
                WHERE molecular_weight(sequence) BETWEEN 15000 AND 50000
                LIMIT 3;
            """)
            for row in cur.fetchall():
                print(row)

            print("\n--- 2. PROSITE Match ---")
            cur.execute("""
                SELECT uniprot_id, substring(sequence from 1 for 20) as seq_start
                FROM proteins
                WHERE prosite_match(sequence, '[ST]-x(2)-[RK]')
                LIMIT 3;
            """)
            matches = cur.fetchall()
            for row in matches:
                print(row)
            
            print("\n--- 3. Center of Mass ---")
            # Get an ID that actually has atoms
            cur.execute("""
                SELECT p.uniprot_id 
                FROM proteins p
                JOIN protein_atoms a ON p.uniprot_id = a.uniprot_id
                WHERE prosite_match(p.sequence, '[ST]-x(2)-[RK]')
                LIMIT 1;
            """)
            res = cur.fetchone()
            if res:
                target_id = res[0]
                print(f"Target: {target_id}")
                cur.execute(f"""
                    SELECT center_of_mass(array_agg(coord))
                    FROM protein_atoms
                    WHERE uniprot_id = '{target_id}';
                """)
                print(cur.fetchone()[0])
            
except Exception as e:
    print(f"Error: {e}")
