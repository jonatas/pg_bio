import psycopg
DB_URI = "postgresql://localhost:28818/bio_demo"
with psycopg.connect(DB_URI, autocommit=True) as conn:
    with conn.cursor() as cur:
        try:
            cur.execute("DROP EXTENSION pg_bio CASCADE;")
            print("Dropped extension.")
        except Exception as e:
            print("Drop failed:", e)
        try:
            cur.execute("CREATE EXTENSION pg_bio;")
            print("Created extension.")
        except Exception as e:
            print("Create failed:", e)
