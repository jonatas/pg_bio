# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
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
            
            # Use pg_bio to find the nearest structural and sequence homologs
            cur.execute("""
                SELECT uniprot_id, substring(name from 1 for 65) as short_name, distance, hybrid_score
                FROM pg_bio_search_homologs(%s, p_max_distance := 0.2, p_limit := 5);
            """, (t_id,))
            
            print("⚠️ POTENTIAL OFF-TARGET SIDE EFFECTS DISCOVERED:")
            for row in cur.fetchall():
                u_id, u_name, distance, h_score = row
                similarity_pct = (1 - (distance / 2.0)) * 100
                
                warning = "  🚨 SEVERE HYBRID MATCH!" if h_score < 0.1 else ""
                
                print(f"   [{u_id}] {u_name}{warning}")
                print(f"       ↳ Structural Distance: {distance:.4f} ({similarity_pct:.1f}% Match) | Hybrid Score: {h_score:.4f}")
            
            print("\n📈 THE MATH BEHIND THIS:")
            print("   Instead of matching letters (A, C, T, G), pg_bio calculates the Cosine Distance")
            print("   between two 320-dimensional arrays (A and B).")
            print("   Equation: 1 - ( (A • B) / (||A|| * ||B||) )")
            print("   A distance close to 0 means the proteins fold and function identically,")
            print("   making them massive risks for off-target drug side effects!")

if __name__ == '__main__':
    run_exploration()
