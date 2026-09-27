import psycopg

DB_URL = "postgresql://localhost:28818/bio_demo"

def validate_teaflon_pipeline():
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            print("Running Validation Suite for Teaflon Multiomics Pipeline...")
            
            # --- Test 1: Genomics (Insertion Site) ---
            cur.execute("""
                SELECT safe.feature_name 
                FROM human_genome safe
                JOIN human_genome promoter ON safe.chromosome = promoter.chromosome
                WHERE safe.feature_type = 'safe_harbor'
                  AND promoter.feature_name = 'KRT14_promoter'
                  AND GREATEST(
                      lower(safe.genomic_range) - upper(promoter.genomic_range),
                      lower(promoter.genomic_range) - upper(safe.genomic_range), 0
                  ) < 5000
            """)
            result_genomics = cur.fetchone()
            assert result_genomics is not None, "Genomics Test Failed: No safe harbor found."
            assert result_genomics[0] == "SafeLocus_B", f"Genomics Test Failed: Expected SafeLocus_B, got {result_genomics[0]}"
            print("✅ Test 1 Passed: Genomics GiST Range Query successfully targets SafeLocus_B.")
            
            # --- Test 2: Proteomics (Z-Order Folding) ---
            cur.execute("""
                SELECT COUNT(*) FROM teaflon_atoms
                WHERE z_curve_index BETWEEN z_order_encode(5.0, -5.0, 20.0) AND z_order_encode(15.0, 5.0, 30.0)
                  AND x BETWEEN 5.0 AND 15.0
                  AND y BETWEEN -5.0 AND 5.0
                  AND z BETWEEN 20.0 AND 30.0
                  AND residue_name = 'ARG_BIND'
            """)
            result_proteomics = cur.fetchone()
            assert result_proteomics[0] > 0, "Proteomics Test Failed: Z-Order curve found no binding pockets."
            print(f"✅ Test 2 Passed: Z-Order Curve successfully matched {result_proteomics[0]} binding atoms.")
            
            # --- Test 3: Transcriptomics (CSR Dosage) ---
            cur.execute("""
                WITH unnested_expr AS (
                    SELECT cell_type, expr_counts.expr
                    FROM single_cell_expression e
                    CROSS JOIN LATERAL unnest(e.expressed_gene_ids) WITH ORDINALITY AS gene_ids(g_id, i1)
                    JOIN LATERAL unnest(e.expression_counts) WITH ORDINALITY AS expr_counts(expr, i2) ON gene_ids.i1 = expr_counts.i2
                    WHERE gene_ids.g_id = 9999
                )
                SELECT cell_type, AVG(expr) FROM unnested_expr GROUP BY cell_type
            """)
            results_transcript = cur.fetchall()
            
            keratin_dosage = next((row[1] for row in results_transcript if "Keratinocyte" in row[0]), 0)
            bone_dosage = next((row[1] for row in results_transcript if "Osteoblast" in row[0]), 0)
            
            assert 15.0 <= keratin_dosage <= 50.0, f"Transcriptomics Test Failed: Keratin dosage {keratin_dosage} out of bounds."
            assert bone_dosage < 2.0, f"Transcriptomics Test Failed: Bone dosage {bone_dosage} is dangerously high (leakage)."
            print(f"✅ Test 3 Passed: Transcriptomics CSR arrays verified Dosage. (Nails: {keratin_dosage:.2f}, Bones: {bone_dosage:.2f})")
            
            print("\n🎉 All Teaflon Multiomics pipelines successfully validated!")

if __name__ == "__main__":
    validate_teaflon_pipeline()
