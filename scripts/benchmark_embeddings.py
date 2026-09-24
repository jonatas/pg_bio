# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "numpy",
#     "scipy",
# ]
# ///

import time
import psycopg
import numpy as np
from scipy.spatial.distance import cosine

DB_URI = "postgresql://localhost:28818/bio_demo"

def run_benchmark():
    print("--- BENCHMARK: Vector Math (pg_bio vs Python) ---")
    
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # Fetch a real embedding from the database to query with
            cur.execute("SELECT embedding FROM proteins WHERE embedding IS NOT NULL LIMIT 1;")
            row = cur.fetchone()
            if not row:
                print("No embeddings found. Seed database first.")
                return
            target_emb_list = row[0]
            target_emb = np.array(target_emb_list, dtype=np.float32)
            
            # Fetch all for Python bench
            cur.execute("SELECT uniprot_id, embedding FROM proteins WHERE embedding IS NOT NULL;")
            all_records = cur.fetchall()
            if not all_records:
                return
                
            db_vectors = {r[0]: np.array(r[1], dtype=np.float32) for r in all_records}
            print(f"Loaded {len(db_vectors)} real biological vectors.\n")
            
            # 1. PYTHON BENCHMARK
            print("[1] Python Pipeline (NumPy/SciPy)")
            start = time.perf_counter()
            best_id, best_dist = None, float('inf')
            for uid, emb in db_vectors.items():
                dist = cosine(target_emb, emb)
                if dist < best_dist and dist > 0.0001: # Avoid self-match
                    best_dist = dist
                    best_id = uid
            py_time = time.perf_counter() - start
            print(f"  Best Match: {best_id} (Distance: {best_dist:.4f})")
            print(f"  Execution Time: {py_time*1000:.2f} ms\n")
            
            # 2. PG_BIO BENCHMARK
            print("[2] pg_bio Pipeline (In-Database Rust)")
            start = time.perf_counter()
            emb_str = "{" + ",".join(map(str, target_emb_list)) + "}"
            cur.execute(f"""
                SELECT uniprot_id, embedding_cosine_distance(embedding, '{emb_str}'::real[]) as distance
                FROM proteins 
                WHERE embedding IS NOT NULL AND embedding_cosine_distance(embedding, '{emb_str}'::real[]) > 0.0001
                ORDER BY distance ASC LIMIT 1;
            """)
            pg_match = cur.fetchone()
            pg_time = time.perf_counter() - start
            
            print(f"  Best Match: {pg_match[0]} (Distance: {pg_match[1]:.4f})")
            print(f"  Execution Time: {pg_time*1000:.2f} ms\n")
            
            # RESULT
            if pg_time > 0:
                print(f"--- pg_bio is {py_time/pg_time:.1f}x faster! ---")

if __name__ == '__main__':
    run_benchmark()
