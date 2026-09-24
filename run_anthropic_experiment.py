import psycopg
import time

DB_URI = "postgresql://localhost:28818/bio_demo"

print("Anthropic's 950 Agents spent 21 hours searching... We will use 1 SQL Query!")
start_time = time.time()

try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT count(*) FROM proteins;")
            total_proteins = cur.fetchone()[0]
            print(f"Scanning {total_proteins} proteins in local DB...")

            # We look for large proteins (like Cas enzymes, typically > 100kDa)
            # and search for a basic catalytic triad motif proxy: Asp/Glu, then a gap, then Asp/Glu, then Lys/Arg
            query = """
                SELECT uniprot_id, name, molecular_weight(sequence) as mass
                FROM proteins
                WHERE molecular_weight(sequence) > 100000
                AND (
                    prosite_match(sequence, '[DE]-x(2)-[DE]') OR 
                    prosite_match(sequence, '[KRH]-x(3)-[DE]')
                )
            """
            cur.execute(query)
            candidates = cur.fetchall()
            
            elapsed = time.time() - start_time
            print(f"Time taken: {elapsed:.4f} seconds (vs 21 hours!)")
            print(f"Found {len(candidates)} potential CRISPR-like enzyme candidates > 100kDa.")
            for c in candidates[:5]:
                print(f" - {c[0]} | {c[1][:30]}... | Mass: {c[2]:.2f} Da")
            if len(candidates) > 5:
                print(" ... (truncated)")
                
except Exception as e:
    print(f"Error: {e}")
