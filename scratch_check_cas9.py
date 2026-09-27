import psycopg
DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT uniprot_id, name, embedding, sequence FROM proteins WHERE name ILIKE '%Cas9%' LIMIT 1;")
            row = cur.fetchone()
            if row:
                uid, name, emb, seq = row
                print(f"Bait: {uid} - {name}")
                cur.execute("""
                WITH closest AS (
                    SELECT uniprot_id, name, sequence, embedding,
                           (embedding <=> %s) as dist
                    FROM proteins
                    ORDER BY embedding <=> %s ASC
                    LIMIT 200
                )
                SELECT c.uniprot_id, c.name, c.dist,
                       (ROW(c.embedding, c.sequence)::bio_feature <~> ROW(%s, %s)::bio_feature) as hybrid_score
                FROM closest c
                WHERE c.name ILIKE %s
                  AND c.dist <= 0.35
                ORDER BY hybrid_score ASC
                LIMIT 3;
                """, (emb, emb, emb, seq, '%uncharacterized%'))
                
                for r in cur.fetchall():
                    print(r)
            else:
                print("No Cas9 found.")
except Exception as e:
    print(e)
