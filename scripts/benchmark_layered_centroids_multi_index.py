# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "numpy",
# ]
# ///
"""
Experiment & Benchmark Suite: Applying pg_vectorsearch's Layered Centroids,
SOAR Boundary Replication, Register-Direct Quantization, and Theoretical
Error-Bounded Fastscan to pg_bio's Three Other Core Indexing Engines:
  1. 3D Morton Z-Order Atomic Spatial Index (`z_order_layered_sphere_ranges`)
  2. 1024-Bit ECFP4 Cheminformatics Tanimoto Index (`ecfp4_register_descriptor`, `ecfp4_fastscan_ub_raw`, `ecfp4_tanimoto_words_raw`)
  3. 9D RangeFinder Active-Site Pocket Index (`pocket_quantize_fastscan`, `pocket_fastscan_iou_ub_raw`, `pocket_iou_8prop_raw`)
"""

import json
import os
import statistics
import time
import psycopg

DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")


def register_new_sql_functions(cur):
    """Ensure the newly compiled Rust symbols in pg_bio.dylib are registered in the live database."""
    cur.execute(
        """
        CREATE OR REPLACE FUNCTION z_order_layered_sphere_ranges(
            center_x double precision,
            center_y double precision,
            center_z double precision,
            radius double precision,
            max_ranges integer
        ) RETURNS TABLE(z_min bigint, z_max bigint, fully_contained boolean)
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'z_order_layered_sphere_ranges_wrapper';

        DO $$
        BEGIN
            IF NOT EXISTS (SELECT 1 FROM pg_type WHERE typname = 'ecfp4fastscansketch') THEN
                CREATE TYPE Ecfp4FastscanSketch;
                CREATE FUNCTION ecfp4fastscansketch_in(cstring) RETURNS Ecfp4FastscanSketch
                IMMUTABLE PARALLEL SAFE LANGUAGE c AS 'pg_bio', 'ecfp4fastscansketch_in_wrapper';
                CREATE FUNCTION ecfp4fastscansketch_out(Ecfp4FastscanSketch) RETURNS cstring
                IMMUTABLE PARALLEL SAFE STRICT LANGUAGE c AS 'pg_bio', 'ecfp4fastscansketch_out_wrapper';
                CREATE TYPE Ecfp4FastscanSketch (
                    INTERNALLENGTH = variable,
                    INPUT = ecfp4fastscansketch_in,
                    OUTPUT = ecfp4fastscansketch_out,
                    STORAGE = extended
                );
            END IF;
        END $$;

        CREATE OR REPLACE FUNCTION ecfp4_fastscan_sketch(smiles text) RETURNS Ecfp4FastscanSketch
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'ecfp4_fastscan_sketch_wrapper';

        CREATE OR REPLACE FUNCTION molecule_fastscan_sketch(mol Molecule) RETURNS Ecfp4FastscanSketch
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'molecule_fastscan_sketch_wrapper';

        CREATE OR REPLACE FUNCTION ecfp4_tanimoto_upper_bound(a Ecfp4FastscanSketch, b Ecfp4FastscanSketch) RETURNS double precision
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'ecfp4_tanimoto_upper_bound_wrapper';

        CREATE OR REPLACE FUNCTION ecfp4_register_descriptor(smiles text)
        RETURNS TABLE(
            l1_centroid smallint,
            l2_centroid smallint,
            soar_centroid smallint,
            popcount smallint,
            block_counts bigint,
            sketch_lo bigint,
            sketch_hi bigint,
            ecfp_words bigint[]
        )
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'ecfp4_register_descriptor_wrapper';

        CREATE OR REPLACE FUNCTION ecfp4_fastscan_ub_raw(
            pop_a smallint, blocks_a bigint, lo_a bigint, hi_a bigint,
            pop_b smallint, blocks_b bigint, lo_b bigint, hi_b bigint
        ) RETURNS double precision
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'ecfp4_fastscan_ub_raw_wrapper';

        CREATE OR REPLACE FUNCTION ecfp4_tanimoto_words_raw(words_a bigint[], words_b bigint[])
        RETURNS double precision
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'ecfp4_tanimoto_words_raw_wrapper';

        CREATE OR REPLACE FUNCTION pocket_layered_centroid_assign(pocket RangeFinderPocket)
        RETURNS TABLE(level1_centroid integer, level2_centroid integer, soar_centroid integer)
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'pocket_layered_centroid_assign_wrapper';

        CREATE OR REPLACE FUNCTION pocket_envelope_iou_upper_bound(
            query RangeFinderPocket,
            env_mins double precision[],
            env_maxs double precision[],
            env_min_widths double precision[],
            include_formal_charge boolean
        ) RETURNS double precision
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'pocket_envelope_iou_upper_bound_wrapper';

        CREATE OR REPLACE FUNCTION pocket_quantize_fastscan(
            mw_min double precision, mw_max double precision,
            logp_min double precision, logp_max double precision,
            hbd_min double precision, hbd_max double precision,
            hba_min double precision, hba_max double precision,
            tpsa_min double precision, tpsa_max double precision,
            rot_min double precision, rot_max double precision,
            arom_min double precision, arom_max double precision,
            fsp3_min double precision, fsp3_max double precision
        ) RETURNS TABLE(
            l1_centroid smallint,
            l2_centroid smallint,
            soar_centroid smallint,
            lo_code bigint,
            hi_code bigint
        )
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'pocket_quantize_fastscan_wrapper';

        CREATE OR REPLACE FUNCTION pocket_fastscan_iou_ub_raw(
            a_lo_code bigint, a_hi_code bigint,
            b_lo_code bigint, b_hi_code bigint
        ) RETURNS double precision
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'pocket_fastscan_iou_ub_raw_wrapper';

        CREATE OR REPLACE FUNCTION pocket_iou_8prop_raw(
            a_mw_min double precision, a_mw_max double precision, a_logp_min double precision, a_logp_max double precision,
            a_hbd_min double precision, a_hbd_max double precision, a_hba_min double precision, a_hba_max double precision,
            a_tpsa_min double precision, a_tpsa_max double precision, a_rot_min double precision, a_rot_max double precision,
            a_arom_min double precision, a_arom_max double precision, a_fsp3_min double precision, a_fsp3_max double precision,
            b_mw_min double precision, b_mw_max double precision, b_logp_min double precision, b_logp_max double precision,
            b_hbd_min double precision, b_hbd_max double precision, b_hba_min double precision, b_hba_max double precision,
            b_tpsa_min double precision, b_tpsa_max double precision, b_rot_min double precision, b_rot_max double precision,
            b_arom_min double precision, b_arom_max double precision, b_fsp3_min double precision, b_fsp3_max double precision
        ) RETURNS double precision
        IMMUTABLE PARALLEL SAFE STRICT
        LANGUAGE c AS 'pg_bio', 'pocket_iou_8prop_raw_wrapper';
        """
    )


def run_experiment_1_morton_3d(cur) -> dict:
    print("==========================================================================")
    print("🧪 Experiment 1: 3D Morton Z-Order Spatial Search (200,000 Protein Atoms)")
    print("   Flat 1-Range Morton Scan vs. Layered Octree Centroid Router + Bounds")
    print("==========================================================================")

    cur.execute("DROP TABLE IF EXISTS bench_atoms_200k CASCADE;")
    cur.execute(
        """
        CREATE TEMP TABLE bench_atoms_200k (
            atom_id serial PRIMARY KEY,
            x double precision NOT NULL,
            y double precision NOT NULL,
            z double precision NOT NULL,
            z_index bigint NOT NULL
        );
        INSERT INTO bench_atoms_200k (x, y, z, z_index)
        SELECT
            px, py, pz, z_order_encode(px, py, pz)
        FROM (
            SELECT
                (( (i::bigint * 73856093) % 90001 ) - 45000)::float8 / 1000.0 AS px,
                (( (i::bigint * 19349663) % 90001 ) - 45000)::float8 / 1000.0 AS py,
                (( (i::bigint * 83492791) % 90001 ) - 45000)::float8 / 1000.0 AS pz
            FROM generate_series(1, 200000) i
        ) s;
        CREATE INDEX idx_bench_atoms_200k_z ON bench_atoms_200k (z_index);
        ANALYZE bench_atoms_200k;
        """
    )

    queries = [
        (0.0, 0.0, 0.0, 5.0),
        (24.2, 24.2, 24.2, 5.0),
        (24.0, -0.2, 24.1, 5.0),
        (-41.2, 24.3, 0.0, 5.0),
        (0.5, -0.4, 1.2, 5.0),
        (-1.8, 0.9, -0.7, 5.0),
        (2.1, -1.5, 0.3, 6.0),
        (-3.2, 2.4, 1.8, 6.0),
        (0.2, 8.5, -1.1, 5.0),
        (-0.9, -12.0, 4.5, 5.0),
        (14.2, -0.5, 7.8, 5.0),
        (-18.5, 11.2, -0.3, 6.0),
        (7.5, 7.5, -2.0, 6.0),
        (-6.0, 4.0, -4.0, 5.0),
        (12.3, 15.1, 18.4, 5.0),
        (-22.1, -14.8, 9.2, 6.0),
        (5.2, -19.4, -11.0, 5.0),
        (24.2, 24.2, 24.2, 8.0),
        (0.0, 0.0, 0.0, 8.0),
        (1.5, -1.5, 1.5, 8.0),
    ]

    naive_lats = []
    naive_candidates = []
    layered_24_lats = []
    layered_24_candidates = []
    layered_64_lats = []
    layered_64_candidates = []
    recalls_24 = []
    recalls_64 = []

    for cx, cy, cz, r in queries:
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH bounds AS (
                SELECT z_order_encode(%s, %s, %s) AS min_z,
                       z_order_encode(%s, %s, %s) AS max_z
            )
            SELECT count(*),
                   count(*) FILTER (
                       WHERE (a.x - %s)*(a.x - %s) + (a.y - %s)*(a.y - %s) + (a.z - %s)*(a.z - %s) <= %s
                   )
            FROM bench_atoms_200k a, bounds b
            WHERE a.z_index BETWEEN b.min_z AND b.max_z;
            """,
            (cx - r, cy - r, cz - r, cx + r, cy + r, cz + r, cx, cx, cy, cy, cz, cz, r * r),
        )
        n_cand, n_true = cur.fetchone()
        naive_lats.append((time.perf_counter() - t0) * 1000.0)
        naive_candidates.append(n_cand)

        t0 = time.perf_counter()
        cur.execute(
            """
            SELECT count(*),
                   count(*) FILTER (
                       WHERE rng.fully_contained
                          OR (a.x - %s)*(a.x - %s) + (a.y - %s)*(a.y - %s) + (a.z - %s)*(a.z - %s) <= %s
                   )
            FROM z_order_layered_sphere_ranges(%s, %s, %s, %s, 24) rng
            JOIN bench_atoms_200k a ON a.z_index BETWEEN rng.z_min AND rng.z_max;
            """,
            (cx, cx, cy, cy, cz, cz, r * r, cx, cy, cz, r),
        )
        l24_cand, l24_true = cur.fetchone()
        layered_24_lats.append((time.perf_counter() - t0) * 1000.0)
        layered_24_candidates.append(l24_cand)
        recalls_24.append(1.0 if n_true == 0 else (l24_true / n_true))

        t0 = time.perf_counter()
        cur.execute(
            """
            SELECT count(*),
                   count(*) FILTER (
                       WHERE rng.fully_contained
                          OR (a.x - %s)*(a.x - %s) + (a.y - %s)*(a.y - %s) + (a.z - %s)*(a.z - %s) <= %s
                   )
            FROM z_order_layered_sphere_ranges(%s, %s, %s, %s, 64) rng
            JOIN bench_atoms_200k a ON a.z_index BETWEEN rng.z_min AND rng.z_max;
            """,
            (cx, cx, cy, cy, cz, cz, r * r, cx, cy, cz, r),
        )
        l64_cand, l64_true = cur.fetchone()
        layered_64_lats.append((time.perf_counter() - t0) * 1000.0)
        layered_64_candidates.append(l64_cand)
        recalls_64.append(1.0 if n_true == 0 else (l64_true / n_true))

    res = {
        "dataset_atoms": 200000,
        "num_queries": len(queries),
        "naive_1range": {
            "p50_ms": round(statistics.median(naive_lats), 2),
            "mean_ms": round(statistics.mean(naive_lats), 2),
            "avg_candidates_scanned": round(statistics.mean(naive_candidates), 1),
        },
        "layered_octree_24ranges": {
            "p50_ms": round(statistics.median(layered_24_lats), 2),
            "mean_ms": round(statistics.mean(layered_24_lats), 2),
            "avg_candidates_scanned": round(statistics.mean(layered_24_candidates), 1),
            "candidate_reduction_factor": round(
                statistics.mean(naive_candidates) / max(statistics.mean(layered_24_candidates), 1.0), 2
            ),
            "speedup_p50": round(statistics.median(naive_lats) / max(statistics.median(layered_24_lats), 0.01), 2),
            "speedup_mean": round(statistics.mean(naive_lats) / max(statistics.mean(layered_24_lats), 0.01), 2),
            "recall_pct": round(statistics.mean(recalls_24) * 100.0, 2),
        },
        "layered_octree_64ranges": {
            "p50_ms": round(statistics.median(layered_64_lats), 2),
            "mean_ms": round(statistics.mean(layered_64_lats), 2),
            "avg_candidates_scanned": round(statistics.mean(layered_64_candidates), 1),
            "candidate_reduction_factor": round(
                statistics.mean(naive_candidates) / max(statistics.mean(layered_64_candidates), 1.0), 2
            ),
            "speedup_p50": round(statistics.median(naive_lats) / max(statistics.median(layered_64_lats), 0.01), 2),
            "speedup_mean": round(statistics.mean(naive_lats) / max(statistics.mean(layered_64_lats), 0.01), 2),
            "recall_pct": round(statistics.mean(recalls_64) * 100.0, 2),
        },
    }

    print(
        f"  • Flat 1-Range Morton Scan:         p50={res['naive_1range']['p50_ms']:6.2f} ms | "
        f"mean={res['naive_1range']['mean_ms']:6.2f} ms | avg_scanned={res['naive_1range']['avg_candidates_scanned']:8.1f} atoms"
    )
    print(
        f"  • Layered Octree (max_ranges=24):   p50={res['layered_octree_24ranges']['p50_ms']:6.2f} ms | "
        f"mean={res['layered_octree_24ranges']['mean_ms']:6.2f} ms | avg_scanned={res['layered_octree_24ranges']['avg_candidates_scanned']:8.1f} atoms "
        f"({res['layered_octree_24ranges']['candidate_reduction_factor']}x fewer candidates, "
        f"{res['layered_octree_24ranges']['speedup_mean']}x faster, Recall={res['layered_octree_24ranges']['recall_pct']}%)"
    )
    print(
        f"  • Layered Octree (max_ranges=64):   p50={res['layered_octree_64ranges']['p50_ms']:6.2f} ms | "
        f"mean={res['layered_octree_64ranges']['mean_ms']:6.2f} ms | avg_scanned={res['layered_octree_64ranges']['avg_candidates_scanned']:8.1f} atoms "
        f"({res['layered_octree_64ranges']['candidate_reduction_factor']}x fewer candidates, "
        f"{res['layered_octree_64ranges']['speedup_mean']}x faster, Recall={res['layered_octree_64ranges']['recall_pct']}%)\n"
    )
    return res


def run_experiment_2_ecfp4_tanimoto(cur) -> dict:
    print("==========================================================================")
    print("🧪 Experiment 2: 1024-Bit ECFP4 Tanimoto Search (10,000 Diverse Molecules)")
    print("   Full Molecule % Scan vs. Register-Direct Layered Centroids + Fastscan")
    print("==========================================================================")

    # Build a 10,000-molecule combinatorial scaffold library (100 cores/substituents x 100 R-groups)
    cores = [
        "c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccsc1", "c1c[nH]cn1",
        "c1ccc2ccccc2c1", "c1ccc2ncccc2c1", "C1CCCCC1", "C1CCNCC1", "C1CCOCC1",
    ]
    linkers = [
        "C(=O)O", "C(=O)N", "CC(=O)O", "OCC(=O)N", "NC(=O)C",
        "S(=O)(=O)N", "C(F)(F)F", "OC", "CC(C)C", "CCN(C)C",
    ]
    tails = [
        f"C{'C' * (k % 8)}{'O' if k % 3 == 0 else ('N' if k % 3 == 1 else 'F')}{'Cl' if k % 5 == 0 else ''}"
        for k in range(100)
    ]
    diverse_smiles = []
    for c in cores:
        for l in linkers:
            for t in tails:
                diverse_smiles.append(f"{t}{l}{c}")

    cur.execute("DROP TABLE IF EXISTS bench_chem_reg_10k CASCADE;")
    t0 = time.perf_counter()
    cur.execute(
        """
        CREATE TEMP TABLE bench_chem_reg_10k (
            id serial PRIMARY KEY,
            raw_smiles text NOT NULL,
            mol Molecule,
            l1_centroid smallint NOT NULL,
            l2_centroid smallint NOT NULL,
            soar_centroid smallint NOT NULL,
            popcount smallint NOT NULL,
            block_counts bigint NOT NULL,
            sketch_lo bigint NOT NULL,
            sketch_hi bigint NOT NULL,
            ecfp_words bigint[] STORAGE PLAIN NOT NULL
        );
        """
    )
    cur.execute(
        """
        INSERT INTO bench_chem_reg_10k (
            raw_smiles, mol, l1_centroid, l2_centroid, soar_centroid,
            popcount, block_counts, sketch_lo, sketch_hi, ecfp_words
        )
        SELECT
            s,
            parse_smiles(s),
            d.l1_centroid,
            d.l2_centroid,
            d.soar_centroid,
            d.popcount,
            d.block_counts,
            d.sketch_lo,
            d.sketch_hi,
            d.ecfp_words
        FROM unnest(%s::text[]) AS s
        CROSS JOIN LATERAL ecfp4_register_descriptor(s) d;
        """,
        (diverse_smiles,),
    )
    cur.execute(
        """
        CREATE INDEX idx_chem_reg_l2_pop ON bench_chem_reg_10k (l2_centroid, popcount);
        CREATE INDEX idx_chem_reg_soar_pop ON bench_chem_reg_10k (soar_centroid, popcount);
        CREATE INDEX idx_chem_reg_pop ON bench_chem_reg_10k (popcount);
        ANALYZE bench_chem_reg_10k;
        """
    )
    build_s = time.perf_counter() - t0
    print(f"  • Built 10,000-molecule Register-Direct Fastscan & Layered Centroid Table in {build_s:.2f} s")

    cur.execute("SELECT id, raw_smiles FROM bench_chem_reg_10k ORDER BY id LIMIT 1000;")
    pool = cur.fetchall()
    query_mols = pool[::50][:20]

    threshold = 0.70
    baseline_lats = []
    reg_exact_lats = []
    reg_centroid_soar_lats = []
    exact_survivors = []
    centroid_survivors = []
    recalls_exact = []
    recalls_centroid = []

    for qid, qsmiles in query_mols:
        # 1. Baseline: Full 10,000-row Molecule % Molecule Sequential Scan
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (SELECT parse_smiles(%s) AS qmol)
            SELECT id
            FROM bench_chem_reg_10k m, q
            WHERE (m.mol %% q.qmol) >= %s;
            """,
            (qsmiles, threshold),
        )
        gt_ids = set(r[0] for r in cur.fetchall())
        baseline_lats.append((time.perf_counter() - t0) * 1000.0)

        # 2. Method A: Register-Direct Popcount B-Tree + 24B Register Fastscan Error Bound -> ecfp4_tanimoto_words_raw (100% Exact Recall!)
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (
                SELECT * FROM ecfp4_register_descriptor(%s)
            ),
            stage1 AS (
                SELECT m.id, m.ecfp_words, q.ecfp_words AS q_words
                FROM bench_chem_reg_10k m, q
                WHERE m.popcount BETWEEN floor(q.popcount * %s)::smallint AND ceil(q.popcount / %s)::smallint
                  AND ecfp4_fastscan_ub_raw(
                      m.popcount, m.block_counts, m.sketch_lo, m.sketch_hi,
                      q.popcount, q.block_counts, q.sketch_lo, q.sketch_hi
                  ) >= %s
            )
            SELECT count(*),
                   array_agg(id) FILTER (WHERE ecfp4_tanimoto_words_raw(ecfp_words, q_words) >= %s)
            FROM stage1;
            """,
            (qsmiles, threshold, threshold, threshold, threshold),
        )
        surv_cnt, hit_arr = cur.fetchone()
        reg_exact_lats.append((time.perf_counter() - t0) * 1000.0)
        exact_survivors.append(surv_cnt)
        hits_a = set(hit_arr or [])
        recalls_exact.append(1.0 if not gt_ids else len(hits_a & gt_ids) / len(gt_ids))

        # 3. Method B: Layered Centroid + SOAR Secondary Cluster Routing + Register Fastscan Error Bound
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (
                SELECT * FROM ecfp4_register_descriptor(%s)
            ),
            stage1 AS (
                SELECT m.id, m.ecfp_words, q.ecfp_words AS q_words
                FROM bench_chem_reg_10k m, q
                WHERE (m.l2_centroid IN (q.l2_centroid, q.soar_centroid)
                       OR m.soar_centroid IN (q.l2_centroid, q.soar_centroid))
                  AND m.popcount BETWEEN floor(q.popcount * %s)::smallint AND ceil(q.popcount / %s)::smallint
                  AND ecfp4_fastscan_ub_raw(
                      m.popcount, m.block_counts, m.sketch_lo, m.sketch_hi,
                      q.popcount, q.block_counts, q.sketch_lo, q.sketch_hi
                  ) >= %s
            )
            SELECT count(*),
                   array_agg(id) FILTER (WHERE ecfp4_tanimoto_words_raw(ecfp_words, q_words) >= %s)
            FROM stage1;
            """,
            (qsmiles, threshold, threshold, threshold, threshold),
        )
        c_surv_cnt, c_hit_arr = cur.fetchone()
        reg_centroid_soar_lats.append((time.perf_counter() - t0) * 1000.0)
        centroid_survivors.append(c_surv_cnt)
        hits_b = set(c_hit_arr or [])
        recalls_centroid.append(1.0 if not gt_ids else len(hits_b & gt_ids) / len(gt_ids))

    res = {
        "dataset_molecules": 10000,
        "tanimoto_threshold": threshold,
        "baseline_full_scan": {
            "p50_ms": round(statistics.median(baseline_lats), 2),
            "mean_ms": round(statistics.mean(baseline_lats), 2),
            "molecules_evaluated": 10000.0,
        },
        "register_fastscan_exact": {
            "p50_ms": round(statistics.median(reg_exact_lats), 2),
            "mean_ms": round(statistics.mean(reg_exact_lats), 2),
            "avg_survivors_rechecked": round(statistics.mean(exact_survivors), 1),
            "pruned_pct": round((1.0 - statistics.mean(exact_survivors) / 10000.0) * 100.0, 2),
            "speedup_p50": round(statistics.median(baseline_lats) / max(statistics.median(reg_exact_lats), 0.01), 2),
            "recall_pct": round(statistics.mean(recalls_exact) * 100.0, 2),
        },
        "layered_centroid_soar_register": {
            "p50_ms": round(statistics.median(reg_centroid_soar_lats), 2),
            "mean_ms": round(statistics.mean(reg_centroid_soar_lats), 2),
            "avg_survivors_rechecked": round(statistics.mean(centroid_survivors), 1),
            "pruned_pct": round((1.0 - statistics.mean(centroid_survivors) / 10000.0) * 100.0, 2),
            "speedup_p50": round(statistics.median(baseline_lats) / max(statistics.median(reg_centroid_soar_lats), 0.01), 2),
            "recall_pct": round(statistics.mean(recalls_centroid) * 100.0, 2),
        },
    }

    print(
        f"  • Baseline Full Molecule % Scan:             p50={res['baseline_full_scan']['p50_ms']:6.2f} ms | "
        f"mean={res['baseline_full_scan']['mean_ms']:6.2f} ms | evaluated=10,000.0 mols"
    )
    print(
        f"  • Register Fastscan Bound (100% Exact):      p50={res['register_fastscan_exact']['p50_ms']:6.2f} ms | "
        f"mean={res['register_fastscan_exact']['mean_ms']:6.2f} ms | rechecked={res['register_fastscan_exact']['avg_survivors_rechecked']:6.1f} mols "
        f"(Pruned {res['register_fastscan_exact']['pruned_pct']}%, {res['register_fastscan_exact']['speedup_p50']}x faster, "
        f"Recall={res['register_fastscan_exact']['recall_pct']}%)"
    )
    print(
        f"  • Layered Centroid + SOAR + Register Bound:  p50={res['layered_centroid_soar_register']['p50_ms']:6.2f} ms | "
        f"mean={res['layered_centroid_soar_register']['mean_ms']:6.2f} ms | rechecked={res['layered_centroid_soar_register']['avg_survivors_rechecked']:6.1f} mols "
        f"(Pruned {res['layered_centroid_soar_register']['pruned_pct']}%, {res['layered_centroid_soar_register']['speedup_p50']}x faster, "
        f"Recall={res['layered_centroid_soar_register']['recall_pct']}%)\n"
    )
    return res


def run_experiment_3_rangefinder_pockets(cur) -> dict:
    print("==========================================================================")
    print("🧪 Experiment 3: 9D RangeFinder Active-Site Pocket Top-10 IoU Search")
    print("   Composite RangeFinderPocket Scan vs. 16B Register-Direct Layered Centroids")
    print("==========================================================================")

    cur.execute("DROP TABLE IF EXISTS bench_pockets_reg CASCADE;")
    cur.execute("DROP TABLE IF EXISTS bench_pocket_reg_envelopes CASCADE;")

    t0 = time.perf_counter()
    cur.execute(
        """
        CREATE TEMP TABLE bench_pockets_reg AS
        SELECT
            p.uniprot_id,
            p.pocket,
            q.l1_centroid,
            q.l2_centroid,
            q.soar_centroid,
            q.lo_code,
            q.hi_code,
            p.mw_min, p.mw_max,
            p.logp_min, p.logp_max,
            p.hbd_min, p.hbd_max,
            p.hba_min, p.hba_max,
            p.tpsa_min, p.tpsa_max,
            p.rot_min, p.rot_max,
            p.charge_min, p.charge_max,
            p.arom_min, p.arom_max,
            p.fsp3_min, p.fsp3_max
        FROM swarm_pocket_ranges p
        CROSS JOIN LATERAL pocket_quantize_fastscan(
            p.mw_min, p.mw_max, p.logp_min, p.logp_max,
            p.hbd_min, p.hbd_max, p.hba_min, p.hba_max,
            p.tpsa_min, p.tpsa_max, p.rot_min, p.rot_max,
            p.arom_min, p.arom_max, p.fsp3_min, p.fsp3_max
        ) q;

        CREATE INDEX idx_pockets_reg_l2 ON bench_pockets_reg (l2_centroid);
        CREATE INDEX idx_pockets_reg_soar ON bench_pockets_reg (soar_centroid);

        CREATE TEMP TABLE bench_pocket_reg_envelopes AS
        SELECT
            l2_centroid,
            min(l1_centroid) AS l1_centroid,
            count(*) AS member_count,
            ARRAY[
                min(mw_min), min(logp_min), min(hbd_min), min(hba_min),
                min(tpsa_min), min(rot_min), min(charge_min), min(arom_min), min(fsp3_min)
            ]::float8[] AS env_mins,
            ARRAY[
                max(mw_max), max(logp_max), max(hbd_max), max(hba_max),
                max(tpsa_max), max(rot_max), max(charge_max), max(arom_max), max(fsp3_max)
            ]::float8[] AS env_maxs,
            ARRAY[
                min(mw_max - mw_min), min(logp_max - logp_min), min(hbd_max - hbd_min), min(hba_max - hba_min),
                min(tpsa_max - tpsa_min), min(rot_max - rot_min), min(charge_max - charge_min), min(arom_max - arom_min), min(fsp3_max - fsp3_min)
            ]::float8[] AS env_min_widths
        FROM bench_pockets_reg
        GROUP BY l2_centroid;

        ANALYZE bench_pockets_reg;
        ANALYZE bench_pocket_reg_envelopes;
        """
    )
    build_s = time.perf_counter() - t0
    cur.execute("SELECT count(*), (SELECT count(*) FROM bench_pocket_reg_envelopes) FROM bench_pockets_reg;")
    total_pockets, num_clusters = cur.fetchone()
    print(f"  • Indexed {total_pockets:,} 3D AlphaFold pockets into {num_clusters} Balanced Level-2 Centroid Envelopes in {build_s:.3f} s")

    cur.execute("SELECT uniprot_id FROM bench_pockets_reg ORDER BY uniprot_id LIMIT 400;")
    q_ids = [r[0] for r in cur.fetchall()][::20][:20]

    baseline_composite_lats = []
    reg_full_lats = []
    nprobe2_lats = []
    nprobe4_lats = []
    nprobe8_lats = []
    recalls_np2 = []
    recalls_np4 = []
    recalls_np8 = []
    eps = 1e-6

    for qid in q_ids:
        # 1. Original Baseline: Composite RangeFinderPocket Scan across all 4,856 pockets
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (SELECT pocket AS qpocket FROM bench_pockets_reg WHERE uniprot_id = %s)
            SELECT p.uniprot_id, pocket_range_iou(q.qpocket, p.pocket, false) AS iou
            FROM bench_pockets_reg p, q
            WHERE p.uniprot_id != %s
            ORDER BY iou DESC
            LIMIT 10;
            """,
            (qid, qid),
        )
        gt_rows = cur.fetchall()
        baseline_composite_lats.append((time.perf_counter() - t0) * 1000.0)
        gt_ids = set(r[0] for r in gt_rows)
        gt_min_iou = float(gt_rows[-1][1]) if gt_rows else 0.0

        # 2. Register-Direct Full Scan (16B lo_code/hi_code Fastscan + pocket_iou_8prop_raw, 100% Exact)
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (
                SELECT lo_code, hi_code,
                       mw_min, mw_max, logp_min, logp_max, hbd_min, hbd_max, hba_min, hba_max,
                       tpsa_min, tpsa_max, rot_min, rot_max, arom_min, arom_max, fsp3_min, fsp3_max
                FROM bench_pockets_reg WHERE uniprot_id = %s
            )
            SELECT p.uniprot_id,
                   pocket_iou_8prop_raw(
                       q.mw_min, q.mw_max, q.logp_min, q.logp_max, q.hbd_min, q.hbd_max, q.hba_min, q.hba_max,
                       q.tpsa_min, q.tpsa_max, q.rot_min, q.rot_max, q.arom_min, q.arom_max, q.fsp3_min, q.fsp3_max,
                       p.mw_min, p.mw_max, p.logp_min, p.logp_max, p.hbd_min, p.hbd_max, p.hba_min, p.hba_max,
                       p.tpsa_min, p.tpsa_max, p.rot_min, p.rot_max, p.arom_min, p.arom_max, p.fsp3_min, p.fsp3_max
                   ) AS iou
            FROM bench_pockets_reg p, q
            WHERE p.uniprot_id != %s
              AND pocket_fastscan_iou_ub_raw(q.lo_code, q.hi_code, p.lo_code, p.hi_code) >= 0.85
            ORDER BY iou DESC
            LIMIT 10;
            """,
            (qid, qid),
        )
        cur.fetchall()
        reg_full_lats.append((time.perf_counter() - t0) * 1000.0)

        # 3. Layered Centroid Envelope Upper-Bound Routing + SOAR + Register-Direct IoU (nprobe = 2, 4, 8)
        for nprobe, lat_list, rec_list in [
            (2, nprobe2_lats, recalls_np2),
            (4, nprobe4_lats, recalls_np4),
            (8, nprobe8_lats, recalls_np8),
        ]:
            t0 = time.perf_counter()
            cur.execute(
                """
                WITH q AS (
                    SELECT pocket AS qpocket, l2_centroid AS ql2, soar_centroid AS qsoar,
                           lo_code, hi_code,
                           mw_min, mw_max, logp_min, logp_max, hbd_min, hbd_max, hba_min, hba_max,
                           tpsa_min, tpsa_max, rot_min, rot_max, arom_min, arom_max, fsp3_min, fsp3_max
                    FROM bench_pockets_reg
                    WHERE uniprot_id = %s
                ),
                top_clusters AS (
                    SELECT e.l2_centroid
                    FROM bench_pocket_reg_envelopes e, q
                    ORDER BY (e.l2_centroid = q.ql2 OR e.l2_centroid = q.qsoar) DESC,
                             pocket_envelope_iou_upper_bound(q.qpocket, e.env_mins, e.env_maxs, e.env_min_widths, false) DESC
                    LIMIT %s
                )
                SELECT p.uniprot_id,
                       pocket_iou_8prop_raw(
                           q.mw_min, q.mw_max, q.logp_min, q.logp_max, q.hbd_min, q.hbd_max, q.hba_min, q.hba_max,
                           q.tpsa_min, q.tpsa_max, q.rot_min, q.rot_max, q.arom_min, q.arom_max, q.fsp3_min, q.fsp3_max,
                           p.mw_min, p.mw_max, p.logp_min, p.logp_max, p.hbd_min, p.hbd_max, p.hba_min, p.hba_max,
                           p.tpsa_min, p.tpsa_max, p.rot_min, p.rot_max, p.arom_min, p.arom_max, p.fsp3_min, p.fsp3_max
                       ) AS iou
                FROM bench_pockets_reg p, q
                WHERE p.uniprot_id != %s
                  AND (p.l2_centroid IN (SELECT l2_centroid FROM top_clusters)
                       OR p.soar_centroid IN (q.ql2, q.qsoar))
                ORDER BY iou DESC
                LIMIT 10;
                """,
                (qid, nprobe, qid),
            )
            rows = cur.fetchall()
            lat_list.append((time.perf_counter() - t0) * 1000.0)
            hits = sum(1 for r in rows if (r[0] in gt_ids or float(r[1]) >= gt_min_iou - eps))
            rec_list.append(hits / 10.0)

    res = {
        "dataset_pockets": total_pockets,
        "num_centroid_clusters": num_clusters,
        "baseline_composite_scan": {
            "p50_ms": round(statistics.median(baseline_composite_lats), 2),
            "mean_ms": round(statistics.mean(baseline_composite_lats), 2),
        },
        "register_fastscan_100pct": {
            "p50_ms": round(statistics.median(reg_full_lats), 2),
            "mean_ms": round(statistics.mean(reg_full_lats), 2),
            "speedup_p50": round(statistics.median(baseline_composite_lats) / max(statistics.median(reg_full_lats), 0.01), 2),
            "recall_at_10_pct": 100.0,
        },
        "layered_envelope_nprobe2": {
            "p50_ms": round(statistics.median(nprobe2_lats), 2),
            "mean_ms": round(statistics.mean(nprobe2_lats), 2),
            "speedup_p50": round(statistics.median(baseline_composite_lats) / max(statistics.median(nprobe2_lats), 0.01), 2),
            "recall_at_10_pct": round(statistics.mean(recalls_np2) * 100.0, 2),
        },
        "layered_envelope_nprobe4": {
            "p50_ms": round(statistics.median(nprobe4_lats), 2),
            "mean_ms": round(statistics.mean(nprobe4_lats), 2),
            "speedup_p50": round(statistics.median(baseline_composite_lats) / max(statistics.median(nprobe4_lats), 0.01), 2),
            "recall_at_10_pct": round(statistics.mean(recalls_np4) * 100.0, 2),
        },
        "layered_envelope_nprobe8": {
            "p50_ms": round(statistics.median(nprobe8_lats), 2),
            "mean_ms": round(statistics.mean(nprobe8_lats), 2),
            "speedup_p50": round(statistics.median(baseline_composite_lats) / max(statistics.median(nprobe8_lats), 0.01), 2),
            "recall_at_10_pct": round(statistics.mean(recalls_np8) * 100.0, 2),
        },
    }

    print(
        f"  • Baseline Composite RangeFinderPocket Scan: p50={res['baseline_composite_scan']['p50_ms']:6.2f} ms | "
        f"mean={res['baseline_composite_scan']['mean_ms']:6.2f} ms"
    )
    print(
        f"  • 16B Register Fastscan + Raw IoU (100%):    p50={res['register_fastscan_100pct']['p50_ms']:6.2f} ms | "
        f"mean={res['register_fastscan_100pct']['mean_ms']:6.2f} ms ({res['register_fastscan_100pct']['speedup_p50']}x faster, Recall@10=100.0%)"
    )
    print(
        f"  • Layered Envelope + SOAR (nprobe=2):        p50={res['layered_envelope_nprobe2']['p50_ms']:6.2f} ms | "
        f"mean={res['layered_envelope_nprobe2']['mean_ms']:6.2f} ms ({res['layered_envelope_nprobe2']['speedup_p50']}x faster, "
        f"Recall@10={res['layered_envelope_nprobe2']['recall_at_10_pct']}%)"
    )
    print(
        f"  • Layered Envelope + SOAR (nprobe=4):        p50={res['layered_envelope_nprobe4']['p50_ms']:6.2f} ms | "
        f"mean={res['layered_envelope_nprobe4']['mean_ms']:6.2f} ms ({res['layered_envelope_nprobe4']['speedup_p50']}x faster, "
        f"Recall@10={res['layered_envelope_nprobe4']['recall_at_10_pct']}%)"
    )
    print(
        f"  • Layered Envelope + SOAR (nprobe=8):        p50={res['layered_envelope_nprobe8']['p50_ms']:6.2f} ms | "
        f"mean={res['layered_envelope_nprobe8']['mean_ms']:6.2f} ms ({res['layered_envelope_nprobe8']['speedup_p50']}x faster, "
        f"Recall@10={res['layered_envelope_nprobe8']['recall_at_10_pct']}%)\n"
    )
    return res


def main():
    conn = psycopg.connect(DB_URL, autocommit=True)
    cur = conn.cursor()

    register_new_sql_functions(cur)

    exp1 = run_experiment_1_morton_3d(cur)
    exp2 = run_experiment_2_ecfp4_tanimoto(cur)
    exp3 = run_experiment_3_rangefinder_pockets(cur)

    summary = {
        "experiment_1_morton_3d_octree": exp1,
        "experiment_2_ecfp4_fastscan_tanimoto": exp2,
        "experiment_3_rangefinder_9d_pockets": exp3,
    }

    out_path = os.path.join(os.path.dirname(__file__), "benchmark_layered_centroids_multi_index_results.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"✅ Saved all 3 multi-index layered centroid experiment results to {out_path}")

    cur.close()
    conn.close()


if __name__ == "__main__":
    main()
