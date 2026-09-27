import psycopg
DB_URI = "postgresql://localhost:28818/bio_demo"

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT uniprot_id, name FROM proteins WHERE name ILIKE '%Argonaute%' LIMIT 1;")
            target = cur.fetchone()
            if target:
                t_id, t_name = target
                print(f"Target: {t_id} ({t_name})")
                query = """
                    SELECT uniprot_id, name, distance, hybrid_score 
                    FROM pg_bio_search_homologs(%s, p_max_distance := 0.5, p_limit := 3);
                """
                cur.execute(query, (t_id,))
                for row in cur.fetchall():
                    print(row)
            else:
                print("No target found.")
except Exception as e:
    print(f"Error: {e}")
