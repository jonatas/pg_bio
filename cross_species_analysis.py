from pgbio import PgBioClient
import psycopg

# Connect to the local pg_bio instance
client = PgBioClient("postgresql://localhost:28818/bio_demo")

query = """
    WITH human_proteins AS (
        SELECT uniprot_id, name, sequence, embedding
        FROM proteins
        WHERE name LIKE '%_HUMAN%'
        LIMIT 100
    ),
    mouse_proteins AS (
        SELECT uniprot_id, name, sequence, embedding
        FROM proteins
        WHERE name LIKE '%_MOUSE%'
    )
    SELECT 
        h.name as human_name,
        m.name as mouse_name,
        embedding_cosine_distance(h.embedding, m.embedding) as distance
    FROM human_proteins h
    CROSS JOIN mouse_proteins m
    ORDER BY distance ASC
    LIMIT 10;
"""

print("Searching for the most structurally conserved proteins between Human and Mouse...")
print("Executing vector embedding cross-join...")

with client._get_conn() as conn:
    with conn.cursor() as cur:
        cur.execute("SELECT count(*) FROM proteins WHERE name LIKE '%_HUMAN%'")
        h_count = cur.fetchone()[0]
        cur.execute("SELECT count(*) FROM proteins WHERE name LIKE '%_MOUSE%'")
        m_count = cur.fetchone()[0]
        print(f"Current DB Stats: {h_count} Human proteins, {m_count} Mouse proteins.")
        
        if m_count > 0:
            cur.execute(query)
            results = cur.fetchall()
            print("\nTop 10 Conserved Human-Mouse Protein Pairs:")
            for row in results:
                print(f"Distance {row[2]:.2e} | {row[0]} <---> {row[1]}")
        else:
            print("\nWaiting for Mouse data to finish loading...")
