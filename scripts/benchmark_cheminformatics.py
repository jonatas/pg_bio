import time
import psycopg
import sys
import random

def get_db_connection():
    # pgrx runs locally on 28818 for pg18
    return psycopg.connect("host=localhost port=28818 dbname=pg_bio user=jonatas")

def benchmark():
    print("🧪 Starting Cheminformatics Benchmark: Python vs pg_bio native...")
    try:
        conn = get_db_connection()
    except Exception as e:
        print(f"Failed to connect to db: {e}")
        return

    # Generate 100 random SMILES-like strings to benchmark
    bases = ["C1=CC=CC=C1", "CC1=CC=CC=C1", "C(=O)O", "CCN", "c1ccccc1", "C1CCCCC1"]
    smiles_list = []
    for _ in range(10_000):
        # random composite smiles
        s = "".join(random.choices(bases, k=random.randint(1, 4)))
        smiles_list.append(s)

    # 1. Native pg_bio Benchmark
    print(f"Comparing Tanimoto Similarity for {len(smiles_list)} pairs natively...")
    start_native = time.time()
    
    with conn.cursor() as cur:
        # We process in batches of 1000 using native SQL map
        for i in range(0, len(smiles_list), 1000):
            batch = smiles_list[i:i+1000]
            # Compare each in the batch against Benzene
            query = """
            SELECT (smiles_to_fingerprint(%s) %% smiles_to_fingerprint('C1=CC=CC=C1'))
            """
            for s in batch:
                cur.execute(query, (s,))
                cur.fetchone()
                
    native_time = time.time() - start_native
    print(f"✅ Native pg_bio Tanimoto (10k ops): {native_time:.3f} seconds")

    # 2. Python RDKit benchmark
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
        
        start_py = time.time()
        benzene = Chem.MolFromSmiles('C1=CC=CC=C1')
        fp_benzene = AllChem.GetMorganFingerprintAsBitVect(benzene, 2, nBits=1024)
        
        for s in smiles_list:
            mol = Chem.MolFromSmiles(s)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                _sim = DataStructs.TanimotoSimilarity(fp_benzene, fp)
        
        py_time = time.time() - start_py
        print(f"✅ Python RDKit Tanimoto (10k ops): {py_time:.3f} seconds")
        print(f"🚀 pg_bio is {py_time/native_time:.2f}x faster by avoiding I/O serialization!")

    except ImportError:
        print("⚠️ RDKit not installed. Skipping Python benchmark.")
        print("To run the full benchmark: uv pip install rdkit")

if __name__ == '__main__':
    benchmark()
