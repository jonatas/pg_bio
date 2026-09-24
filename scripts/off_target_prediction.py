# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "pgbio @ file://./pgbio-py",
# ]
# ///
import psycopg
import math

DB_URI = "postgresql://localhost:28818/bio_demo"

def run_exploration():
    print("🔬 REAL-WORLD SCENARIO: Off-Target Drug Effect Prediction\n")
    
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # Find a clinically significant protein (e.g., a Receptor or Kinase target)
            cur.execute("""
                SELECT uniprot_id, name, embedding 
                FROM proteins 
                WHERE name ILIKE '%receptor%' OR name ILIKE '%kinase%'
                LIMIT 1;
            """)
            target = cur.fetchone()

            t_id, t_name, t_emb = target
            print(f"🎯 DRUG TARGET (The protein our drug is designed to bind to):")
            print(f"   ID: {t_id}")
            print(f"   Name: {t_name}")
            print(f"   Math: Our target is a point in a 320-dimensional space.\n")
            
            print("🚀 RUNNING pg_bio VECTOR SEARCH ACROSS THE HUMAN PROTEOME...")
            print("   (Looking for structurally identical proteins that our drug might accidentally hit...)\n")
            
            # Use pg_bio to find the nearest vectors mathematically
            cur.execute("""
                SELECT 
                    uniprot_id, 
                    substring(name from 1 for 65) as name, 
                    embedding_cosine_distance(embedding, %s::real[]) as distance
                FROM proteins 
                WHERE uniprot_id != %s AND embedding IS NOT NULL
                ORDER BY distance ASC
                LIMIT 5;
            """, (t_emb, t_id))
            
            print("⚠️ POTENTIAL OFF-TARGET SIDE EFFECTS DISCOVERED:")
            for row in cur.fetchall():
                # Cosine distance to similarity percentage
                distance = row[2]
                similarity_pct = (1 - (distance / 2.0)) * 100
                print(f"   [{row[0]}] {row[1]}")
                print(f"       ↳ Cosine Distance: {distance:.4f} ({similarity_pct:.1f}% Structural Match)")
            
            print("\n📈 THE MATH BEHIND THIS:")
            print("   Instead of matching letters (A, C, T, G), pg_bio calculates the Cosine Distance")
            print("   between two 320-dimensional arrays (A and B).")
            print("   Equation: 1 - ( (A • B) / (||A|| * ||B||) )")
            print("   A distance close to 0 means the proteins fold and function identically,")
            print("   making them massive risks for off-target drug side effects!")

if __name__ == '__main__':
    run_exploration()
