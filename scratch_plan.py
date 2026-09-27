import psycopg
DB_URI = "postgresql://localhost:28818/bio_demo"

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # 1. Fetch a target ID
            cur.execute("SELECT uniprot_id FROM proteins WHERE name ILIKE '%Argonaute%' LIMIT 1;")
            t_id = cur.fetchone()[0]
            
            # 2. Explain Analyze with the WHERE inside the index block
            print("--- WITH WHERE INSIDE ---")
            cur.execute("EXPLAIN ANALYZE SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM proteins WHERE uniprot_id = %s)) FROM proteins p WHERE p.uniprot_id != %s ORDER BY p.embedding <=> (SELECT embedding FROM proteins WHERE uniprot_id = %s) ASC LIMIT 50;", (t_id, t_id, t_id))
            for row in cur.fetchall():
                print(row[0])
                
            # 3. Explain Analyze with the WHERE moved to outer block
            print("\n--- WITH WHERE OUTSIDE ---")
            cur.execute("EXPLAIN ANALYZE SELECT * FROM (SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM proteins WHERE uniprot_id = %s)) as dist FROM proteins p ORDER BY p.embedding <=> (SELECT embedding FROM proteins WHERE uniprot_id = %s) ASC LIMIT 50) c WHERE c.uniprot_id != %s;", (t_id, t_id, t_id))
            for row in cur.fetchall():
                print(row[0])
except Exception as e:
    print(f"Error: {e}")
