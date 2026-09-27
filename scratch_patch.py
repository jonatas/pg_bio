import sys

with open("scripts/off_target_prediction.py", "r") as f:
    content = f.read()

old_query = """            # Use pg_bio to find the nearest vectors mathematically
            cur.execute(\"\"\"
                SELECT 
                    uniprot_id, 
                    substring(name from 1 for 65) as name, 
                    embedding_cosine_distance(embedding, %s::real[]) as distance
                FROM proteins 
                WHERE uniprot_id != %s AND embedding IS NOT NULL
                ORDER BY distance ASC
                LIMIT 5;
            \"\"\", (t_emb, t_id))
            
            print("⚠️ POTENTIAL OFF-TARGET SIDE EFFECTS DISCOVERED:")
            for row in cur.fetchall():
                # Cosine distance to similarity percentage
                distance = row[2]
                similarity_pct = (1 - (distance / 2.0)) * 100
                print(f"   [{row[0]}] {row[1]}")
                print(f"       ↳ Cosine Distance: {distance:.4f} ({similarity_pct:.1f}% Structural Match)")"""

new_query = """            # Use pg_bio to find the nearest structural and sequence homologs
            cur.execute(\"\"\"
                SELECT uniprot_id, substring(name from 1 for 65) as short_name, distance, hybrid_score
                FROM pg_bio_search_homologs(%s, p_max_distance := 0.2, p_limit := 5);
            \"\"\", (t_id,))
            
            print("⚠️ POTENTIAL OFF-TARGET SIDE EFFECTS DISCOVERED:")
            for row in cur.fetchall():
                u_id, u_name, distance, h_score = row
                similarity_pct = (1 - (distance / 2.0)) * 100
                
                warning = "  🚨 SEVERE HYBRID MATCH!" if h_score < 0.1 else ""
                
                print(f"   [{u_id}] {u_name}{warning}")
                print(f"       ↳ Structural Distance: {distance:.4f} ({similarity_pct:.1f}% Match) | Hybrid Score: {h_score:.4f}")"""

content = content.replace(old_query, new_query)

with open("scripts/off_target_prediction.py", "w") as f:
    f.write(content)
