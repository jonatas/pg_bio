# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "psycopg>=3.1",
#     "numpy",
#     "scipy",
# ]
# ///

import psycopg
import numpy as np
import time
from scipy.spatial.distance import cosine

DB_URI = "postgresql://localhost:28818/bio_demo"

# ==============================================================================
# Python Adapters
# ==============================================================================

class PgBioAdapter:
    """A fast Python adapter utilizing the pg_bio PostgreSQL extension natively."""
    def __init__(self, db_uri):
        self.db_uri = db_uri
        
    def find_similar_proteins(self, target_embedding, limit=5):
        """Uses pg_bio's native `embedding_cosine_distance` to search inside the DB."""
        with psycopg.connect(self.db_uri) as conn:
            with conn.cursor() as cur:
                # The heavy vector math happens inside Rust/Postgres!
                query = """
                    SELECT name, embedding_cosine_distance(embedding, %s::real[]) as distance 
                    FROM proteins 
                    ORDER BY distance ASC 
                    LIMIT %s;
                """
                # psycopg automatically converts numpy arrays to Postgres REAL[] arrays
                cur.execute(query, (target_embedding.tolist(), limit))
                return cur.fetchall()


class PythonBruteForceAdapter:
    """The standard Data Science approach (pulling data into Python to do math)."""
    def __init__(self, db_uri):
        self.db_uri = db_uri
        
    def find_similar_proteins(self, target_embedding, limit=5):
        """Fetches all embeddings over the network and uses Numpy for math."""
        with psycopg.connect(self.db_uri) as conn:
            with conn.cursor() as cur:
                # 1. I/O Bottleneck: Fetching thousands of vectors over the network
                cur.execute("SELECT name, embedding FROM proteins WHERE embedding IS NOT NULL;")
                rows = cur.fetchall()
        
        # 2. Compute Bottleneck: Calculating cosine distance using SciPy/Numpy
        results = []
        for name, emb in rows:
            dist = cosine(target_embedding, emb)
            results.append((name, dist))
            
        # 3. Sort in python
        results.sort(key=lambda x: x[1])
        return results[:limit]


# ==============================================================================
# Benchmark
# ==============================================================================

def run_benchmark():
    print("Initializing Database Adapters...")
    adapter_pg = PgBioAdapter(DB_URI)
    adapter_py = PythonBruteForceAdapter(DB_URI)

    # Generate a mock 320-dimensional embedding (Simulating an ESM-2 vector query)
    target_emb = np.random.rand(320).astype(np.float32)

    print("\n--- Standard Python / Numpy Data Science Approach ---")
    print("Pulling all vectors from DB and calculating in Pandas/Numpy...")
    start = time.perf_counter()
    py_results = adapter_py.find_similar_proteins(target_emb)
    py_time = (time.perf_counter() - start) * 1000
    print(f"⏱️ Time: {py_time:.2f} ms")
    
    print("\n--- pg_bio (Postgres + Rust Native) Approach ---")
    print("Pushing the query down to the Database engine...")
    start = time.perf_counter()
    pg_results = adapter_pg.find_similar_proteins(target_emb)
    pg_time = (time.perf_counter() - start) * 1000
    print(f"⏱️ Time: {pg_time:.2f} ms")
    
    speedup = py_time / pg_time if pg_time > 0 else 0
    print(f"\n🚀 RESULT: pg_bio is {speedup:.1f}x faster than standard Python!")
    print("Why? Standard Python approaches suffer from massive I/O serialization penalties.")
    print("You have to transmit megabytes of float arrays over the connection just to do math on them.")
    print("pg_bio processes the math securely inside the C/Rust database engine memory, returning only the 5 closest matches!")

if __name__ == "__main__":
    run_benchmark()
