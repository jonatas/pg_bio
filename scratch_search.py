import psycopg
DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # Get the embedding for bait
            cur.execute("SELECT embedding, sequence FROM proteins WHERE uniprot_id = 'A0A1Y3GFM2'")
            row = cur.fetchone()
            if row:
                emb, seq = row
                print("Running Oversample and Refine query...")
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
                LIMIT 5;
                """, (emb, emb, emb, seq, '%uncharacterized%'))
                
                results = cur.fetchall()
                for r in results:
                    print(r)
            else:
                print("Bait not found")
except Exception as e:
    print(e)
