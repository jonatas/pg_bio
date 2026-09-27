import psycopg
DB_URI = "postgresql://localhost:28818/bio_demo"

def run():
    try:
        with psycopg.connect(DB_URI) as conn:
            with conn.cursor() as cur:
                print("--- TABLES ---")
                cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='public';")
                for row in cur.fetchall(): print(row[0])
                
                print("\n--- VECTOR QUERY ---")
                cur.execute("""
                    WITH target_protein AS (
                        SELECT uniprot_id, embedding FROM proteins 
                        WHERE name ILIKE '%dopamine receptor%' AND embedding IS NOT NULL LIMIT 1
                    )
                    SELECT 
                        p.uniprot_id, 
                        substring(p.name from 1 for 45) as name, 
                        embedding_cosine_distance(p.embedding, t.embedding) as distance
                    FROM proteins p, target_protein t
                    WHERE p.uniprot_id != t.uniprot_id AND p.embedding IS NOT NULL
                    ORDER BY distance ASC
                    LIMIT 5;
                """)
                for row in cur.fetchall(): print(row)
                
    except Exception as e:
        print(f"Error: {e}")

run()
