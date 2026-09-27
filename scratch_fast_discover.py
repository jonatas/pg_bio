import psycopg

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

def find_discoveries(family, limit=3):
    print(f"\n--- Scanning for {family} ---")
    try:
        with psycopg.connect(DB_URI) as conn:
            with conn.cursor() as cur:
                # Find a few known targets
                cur.execute("SELECT uniprot_id, name, embedding, sequence FROM proteins WHERE name ILIKE %s LIMIT 5", (f"%{family}%",))
                targets = cur.fetchall()
                for target in targets:
                    t_id, t_name, emb, seq = target
                    print(f"Bait: {t_id} - {t_name.split(' OS=')[0]}")
                    
                    cur.execute("""
                        WITH closest AS (
                            SELECT uniprot_id, name, sequence, embedding,
                                   (embedding <=> %s) as dist
                            FROM proteins
                            ORDER BY embedding <=> %s ASC
                            LIMIT 100
                        )
                        SELECT c.uniprot_id, c.name, c.dist,
                               (ROW(c.embedding, c.sequence)::bio_feature <~> ROW(%s, %s)::bio_feature) as hybrid_score
                        FROM closest c
                        WHERE c.name ILIKE %s
                          AND c.dist <= 0.35
                        ORDER BY hybrid_score ASC
                        LIMIT 1;
                    """, (emb, emb, emb, seq, '%uncharacterized%'))
                    
                    hit = cur.fetchone()
                    if hit:
                        u_id, u_name, dist, h_score = hit
                        print(f"  -> HIT: {u_id} ({u_name.split(' OS=')[0]}) | Dist: {dist:.4f} | Score: {h_score:.4f}")
                    else:
                        print("  -> No close uncharacterized matches found.")
    except Exception as e:
        print("Error:", e)

find_discoveries("Kinase")
find_discoveries("Polymerase")
