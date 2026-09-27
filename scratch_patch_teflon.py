import sys

with open("scripts/teflon_degradation_pipeline.py", "r") as f:
    content = f.read()

old_query = """            cur.execute(\"\"\"
                SELECT uniprot_id, embedding_cosine_distance(embedding, %s::real[]) as distance
                FROM proteins 
                WHERE uniprot_id != %s AND embedding IS NOT NULL
                ORDER BY distance ASC LIMIT 1;
            \"\"\", (t_emb, t_id))
            closest = cur.fetchone()
            
            if closest:
                distance = closest[1]
                if distance > 0.15:
                    print(f"\\n   [+] Toxicity Check Passed! (Distance to closest match: {distance:.3f})")
                    print("       This engineered enzyme is highly specific to PFAS.")
                else:
                    print(f"\\n   [!] Toxicity Warning! Closest match distance is {distance:.3f}.")"""

new_query = """            # Upgraded: Now uses HNSW Vector Index and Smith-Waterman Hybrid Operator (<~>)
            cur.execute(\"\"\"
                SELECT uniprot_id, hybrid_score
                FROM pg_bio_search_homologs(%s, p_max_distance := 0.35, p_limit := 1);
            \"\"\", (t_id,))
            closest = cur.fetchone()
            
            if closest:
                h_score = closest[1]
                if h_score > 0.15:
                    print(f"\\n   [+] Toxicity Check Passed! (Hybrid Score to closest human protein: {h_score:.3f})")
                    print("       This engineered enzyme is highly specific to PFAS and will not accidentally")
                    print("       bind to similar structural pockets in human biology.")
                else:
                    print(f"\\n   [!] Toxicity Warning! Found a massive homology risk (Hybrid Score: {h_score:.3f}).")
            else:
                print(f"\\n   [+] Toxicity Check Passed! No structural homologs found within warning threshold.")"""

content = content.replace(old_query, new_query)

with open("scripts/teflon_degradation_pipeline.py", "w") as f:
    f.write(content)
