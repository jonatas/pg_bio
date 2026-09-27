import psycopg
DB_URI = "postgresql://localhost:28818/bio_demo"

def run():
    try:
        with psycopg.connect(DB_URI) as conn:
            with conn.cursor() as cur:
                print("\n--- SCHEMA FOR PROTEINS ---")
                cur.execute("SELECT column_name, data_type FROM information_schema.columns WHERE table_name = 'proteins';")
                for row in cur.fetchall(): print(row)
                
                print("\n--- VECTOR QUERY WITH CAST ---")
                cur.execute("""
                    WITH target_protein AS (
                        SELECT uniprot_id, embedding FROM proteins 
                        WHERE name ILIKE '%receptor%' AND embedding IS NOT NULL LIMIT 1
                    )
                    SELECT 
                        p.uniprot_id, 
                        substring(p.name from 1 for 45) as name, 
                        embedding_cosine_distance(p.embedding::real[], t.embedding::real[]) as distance
                    FROM proteins p, target_protein t
                    WHERE p.uniprot_id != t.uniprot_id AND p.embedding IS NOT NULL
                    ORDER BY distance ASC
                    LIMIT 5;
                """)
                for row in cur.fetchall(): print(row)

                print("\n--- TEST ATOMS ---")
                cur.execute("SELECT * FROM protein_atoms LIMIT 1;")
                print(cur.fetchall())
                
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    run()
