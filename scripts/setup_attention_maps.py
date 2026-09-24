# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "pgbio @ file://./pgbio-py",
# ]
# ///
import psycopg
import random
from pgbio import PgBioClient

DB_URI = "postgresql://localhost:28818/bio_demo"

def setup_attention_table():
    print("Connecting to 'bio_demo' to setup Sparse Attention Maps...")
    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS protein_attention_maps (
                    id SERIAL PRIMARY KEY,
                    uniprot_id VARCHAR(20) REFERENCES proteins(uniprot_id),
                    attention_data SparseAttentionMap
                );
            """)
            print("Table 'protein_attention_maps' created.")

            # Let's seed an attention map for a known protein, e.g., '1F5X' if it exists.
            cur.execute("SELECT uniprot_id, length(sequence) FROM proteins WHERE name LIKE 'PDB Structure%' LIMIT 5;")
            targets = cur.fetchall()
            
            for uniprot_id, seq_len in targets:
                print(f"Generating Sparse Attention Map for {uniprot_id} (Length: {seq_len})...")
                
                sources = []
                targets_list = []
                weights = []
                
                # Simulate sparse interactions (only ~5 interactions per residue)
                for i in range(1, seq_len + 1):
                    # Interacts with itself (diagonal)
                    sources.append(i)
                    targets_list.append(i)
                    weights.append(1.0)
                    
                    # Interacts with 4 random distant residues
                    for _ in range(4):
                        sources.append(i)
                        targets_list.append(random.randint(1, seq_len))
                        weights.append(random.uniform(0.5, 0.9))
                
                print(f"  -> Generated {len(sources)} non-zero interactions.")
                
                cur.execute("""
                    INSERT INTO protein_attention_maps (uniprot_id, attention_data)
                    VALUES (%s, create_sparse_map(%s::integer, %s::int[], %s::int[], %s::real[]))
                """, (uniprot_id, seq_len, sources, targets_list, weights))
            
            print("Successfully seeded Sparse Attention Maps!")

if __name__ == "__main__":
    setup_attention_table()
