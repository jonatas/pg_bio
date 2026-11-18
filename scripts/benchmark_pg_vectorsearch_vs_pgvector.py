# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "numpy",
# ]
# ///
"""
Comprehensive Benchmark: pgvector (HNSW, IVFFlat) vs. pg_vectorsearch (PRISM)
and Combined pg_bio Multi-Stage Homology Re-Ranking (<~> Hybrid Operator).
"""

import json
import os
import statistics
import time
import psycopg

DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
SAMPLE_SIZE = int(os.getenv("BENCH_SAMPLE_SIZE", "100000"))
NUM_QUERIES = int(os.getenv("BENCH_NUM_QUERIES", "30"))


def get_size_bytes(cur, relname: str) -> int:
    cur.execute("SELECT pg_relation_size(%s::regclass);", (relname,))
    return cur.fetchone()[0]


def get_table_storage_breakdown(cur, table_name: str) -> dict:
    cur.execute(
        """
        SELECT
            pg_relation_size(c.oid, 'main') AS main_bytes,
            COALESCE(pg_relation_size(c.reltoastrelid), 0) AS toast_bytes,
            pg_total_relation_size(c.oid) - pg_indexes_size(c.oid) AS table_total_bytes
        FROM pg_class c
        WHERE c.relname = %s;
        """,
        (table_name,),
    )
    row = cur.fetchone()
    return {
        "main_mb": round(row[0] / (1024 * 1024), 2),
        "toast_mb": round(row[1] / (1024 * 1024), 2),
        "table_total_mb": round(row[2] / (1024 * 1024), 2),
    }


def run_benchmark():
    print(f"======================================================================")
    print(f"🧬 pg_bio Vector Engine Benchmark: pgvector vs. pg_vectorsearch (PRISM)")
    print(f"   Dataset: {SAMPLE_SIZE:,} real UniProt proteins (1280-D ESM-2 halfvec)")
    print(f"   Queries: {NUM_QUERIES} proteins evaluated for Recall@10, Recall@50 & Hybrid <~>")
    print(f"======================================================================\n")

    conn = psycopg.connect(DB_URL, autocommit=True)
    cur = conn.cursor()

    cur.execute("CREATE EXTENSION IF NOT EXISTS vector;")
    cur.execute("CREATE EXTENSION IF NOT EXISTS pg_vectorsearch;")
    cur.execute("CREATE EXTENSION IF NOT EXISTS pg_bio;")
    cur.execute("SET maintenance_work_mem = '2GB';")
    cur.execute("SET max_parallel_maintenance_workers = 4;")

    # 1. Compare STORAGE EXTERNAL (default, TOASTed) vs STORAGE PLAIN (inline)
    print("[Phase 1] Preparing benchmark tables & testing TOAST vs. STORAGE PLAIN...")
    cur.execute("DROP TABLE IF EXISTS bench_proteins_external CASCADE;")
    cur.execute("DROP TABLE IF EXISTS bench_proteins_plain CASCADE;")

    t0 = time.perf_counter()
    cur.execute(
        f"""
        CREATE TABLE bench_proteins_external (
            uniprot_id varchar(20) PRIMARY KEY,
            name text,
            sequence text,
            embedding halfvec(1280)
        );
        INSERT INTO bench_proteins_external
        SELECT uniprot_id, name, sequence, embedding
        FROM proteins
        WHERE embedding IS NOT NULL AND sequence IS NOT NULL
        LIMIT {SAMPLE_SIZE};
        ANALYZE bench_proteins_external;
        """
    )
    t_ext = time.perf_counter() - t0

    t0 = time.perf_counter()
    cur.execute(
        f"""
        CREATE TABLE bench_proteins_plain (
            uniprot_id varchar(20) PRIMARY KEY,
            name text,
            sequence text,
            embedding halfvec(1280) STORAGE PLAIN
        );
        INSERT INTO bench_proteins_plain
        SELECT uniprot_id, name, sequence, embedding
        FROM bench_proteins_external;
        ANALYZE bench_proteins_plain;
        """
    )
    t_plain = time.perf_counter() - t0

    ext_storage = get_table_storage_breakdown(cur, "bench_proteins_external")
    plain_storage = get_table_storage_breakdown(cur, "bench_proteins_plain")
    print(f"  • STORAGE EXTERNAL (default): Main={ext_storage['main_mb']} MB, TOAST={ext_storage['toast_mb']} MB, Total={ext_storage['table_total_mb']} MB (load: {t_ext:.2f}s)")
    print(f"  • STORAGE PLAIN (inline):     Main={plain_storage['main_mb']} MB, TOAST={plain_storage['toast_mb']} MB, Total={plain_storage['table_total_mb']} MB (load: {t_plain:.2f}s)")

    # Pick deterministic query proteins with diverse embeddings
    cur.execute(
        f"""
        SELECT uniprot_id
        FROM bench_proteins_plain
        WHERE length(sequence) BETWEEN 120 AND 450
        ORDER BY md5(uniprot_id)
        LIMIT {NUM_QUERIES};
        """
    )
    query_ids = [r[0] for r in cur.fetchall()]

    # Measure PRISM build & query latency on STORAGE EXTERNAL vs STORAGE PLAIN
    print("\n[Phase 1B] Measuring TOAST (EXTERNAL) vs. Inline (PLAIN) impact on PRISM reranking...")
    t0 = time.perf_counter()
    cur.execute("CREATE INDEX idx_ext_prism ON bench_proteins_external USING prism (embedding vec16_cosine_ops);")
    prism_ext_build_s = time.perf_counter() - t0
    ext_lats = []
    cur.execute("SET prism.nprobe = 20;")
    for qid in query_ids:
        t0 = time.perf_counter()
        cur.execute(
            """
            SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM bench_proteins_external WHERE uniprot_id = %s)) AS dist
            FROM bench_proteins_external p
            ORDER BY p.embedding <=> (SELECT embedding FROM bench_proteins_external WHERE uniprot_id = %s) ASC
            LIMIT 51;
            """,
            (qid, qid),
        )
        cur.fetchall()
        ext_lats.append((time.perf_counter() - t0) * 1000.0)
    print(f"  • PRISM on STORAGE EXTERNAL (TOASTed): build={prism_ext_build_s:.2f}s, query p50={statistics.median(ext_lats):.2f} ms, mean={statistics.mean(ext_lats):.2f} ms")
    cur.execute("DROP TABLE IF EXISTS bench_proteins_external CASCADE;")

    print(f"\n[Phase 2] Computing Exact Brute-Force Ground Truth (Top-50 & Hybrid <~> Top-5) for {len(query_ids)} queries...")

    ground_truth_10 = {}
    ground_truth_50 = {}
    ground_truth_hybrid_5 = {}
    exact_latencies_ms = []

    cur.execute("SET enable_indexscan = off; SET enable_bitmapscan = off;")
    for idx, qid in enumerate(query_ids):
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH target AS (
                SELECT embedding AS emb FROM bench_proteins_plain WHERE uniprot_id = %s
            )
            SELECT p.uniprot_id, (p.embedding <=> t.emb) AS dist
            FROM bench_proteins_plain p, target t
            WHERE p.uniprot_id != %s
            ORDER BY p.embedding <=> t.emb ASC
            LIMIT 50;
            """,
            (qid, qid),
        )
        rows = cur.fetchall()
        dt_ms = (time.perf_counter() - t0) * 1000.0
        exact_latencies_ms.append(dt_ms)
        ground_truth_10[qid] = {
            "ids": set(r[0] for r in rows[:10]),
            "max_dist": float(rows[9][1]) if len(rows) >= 10 else 2.0,
        }
        ground_truth_50[qid] = {
            "ids": set(r[0] for r in rows[:50]),
            "max_dist": float(rows[49][1]) if len(rows) >= 50 else 2.0,
        }

        # For first 15 queries, also compute exact Hybrid <~> top-5 over top-100 true cosine candidates
        if idx < 15:
            cur.execute(
                """
                WITH target AS (
                    SELECT uniprot_id, embedding::vector AS emb, sequence AS seq, length(sequence)::float8 AS len
                    FROM bench_proteins_plain
                    WHERE uniprot_id = %s
                ),
                candidates AS (
                    SELECT p.uniprot_id, p.embedding::vector AS emb, p.sequence AS seq,
                           (p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s)) AS cos_dist
                    FROM bench_proteins_plain p
                    ORDER BY p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s) ASC
                    LIMIT 100
                )
                SELECT c.uniprot_id,
                       (ROW(c.emb, c.seq)::bio_feature <~> ROW(t.emb, t.seq)::bio_feature) AS hybrid_score
                FROM candidates c, target t
                WHERE c.uniprot_id != t.uniprot_id
                  AND length(c.seq) BETWEEN t.len * 0.7 AND t.len * 1.3
                  AND sequence_entropy(c.seq) >= 2.0
                ORDER BY hybrid_score ASC
                LIMIT 5;
                """,
                (qid, qid, qid),
            )
            h_rows = cur.fetchall()
            ground_truth_hybrid_5[qid] = {
                "ids": set(r[0] for r in h_rows),
                "max_score": float(h_rows[-1][1]) if h_rows else 2.0,
                "count": len(h_rows),
            }

    cur.execute("SET enable_indexscan = on; SET enable_bitmapscan = on;")
    print(
        f"  • Brute-force Sequential Scan: mean={statistics.mean(exact_latencies_ms):.2f} ms, "
        f"p50={statistics.median(exact_latencies_ms):.2f} ms"
    )

    def evaluate_queries(guc_setup_sql: str, k: int = 50):
        if guc_setup_sql:
            cur.execute(guc_setup_sql)
        latencies = []
        recalls_10 = []
        recalls_50 = []
        eps = 1e-5

        for qid in query_ids[:3]:
            cur.execute(
                """
                SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s)) AS dist
                FROM bench_proteins_plain p
                ORDER BY p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s) ASC
                LIMIT %s;
                """,
                (qid, qid, k + 1),
            )
            cur.fetchall()

        for qid in query_ids:
            t0 = time.perf_counter()
            cur.execute(
                """
                SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s)) AS dist
                FROM bench_proteins_plain p
                ORDER BY p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s) ASC
                LIMIT %s;
                """,
                (qid, qid, k + 1),
            )
            rows = [r for r in cur.fetchall() if r[0] != qid][:k]
            dt_ms = (time.perf_counter() - t0) * 1000.0
            latencies.append(dt_ms)

            gt10 = ground_truth_10[qid]
            hits_10 = sum(
                1 for r in rows[:10] if (r[0] in gt10["ids"] or float(r[1]) <= gt10["max_dist"] + eps)
            )
            recalls_10.append(hits_10 / 10.0)

            gt50 = ground_truth_50[qid]
            hits_50 = sum(
                1 for r in rows[:50] if (r[0] in gt50["ids"] or float(r[1]) <= gt50["max_dist"] + eps)
            )
            recalls_50.append(hits_50 / 50.0)

        latencies_sorted = sorted(latencies)
        p95_idx = min(int(len(latencies_sorted) * 0.95), len(latencies_sorted) - 1)
        return {
            "mean_ms": round(statistics.mean(latencies), 2),
            "p50_ms": round(statistics.median(latencies), 2),
            "p95_ms": round(latencies_sorted[p95_idx], 2),
            "qps": round(1000.0 / statistics.mean(latencies), 1),
            "recall_at_10": round(statistics.mean(recalls_10) * 100.0, 2),
            "recall_at_50": round(statistics.mean(recalls_50) * 100.0, 2),
        }

    def evaluate_hybrid_pg_bio_pipeline(guc_setup_sql: str, oversample: int = 50):
        if guc_setup_sql:
            cur.execute(guc_setup_sql)
        latencies = []
        hybrid_recalls = []
        eps = 1e-5
        for qid in query_ids[:15]:
            t0 = time.perf_counter()
            cur.execute(
                """
                WITH target AS (
                    SELECT uniprot_id, embedding::vector AS emb, sequence AS seq, length(sequence)::float8 AS len
                    FROM bench_proteins_plain
                    WHERE uniprot_id = %s
                ),
                candidates AS (
                    SELECT p.uniprot_id, p.name, p.embedding::vector AS emb, p.sequence AS seq,
                           (p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s)) AS cos_dist
                    FROM bench_proteins_plain p
                    ORDER BY p.embedding <=> (SELECT embedding FROM bench_proteins_plain WHERE uniprot_id = %s) ASC
                    LIMIT %s
                )
                SELECT c.uniprot_id,
                       c.cos_dist,
                       (ROW(c.emb, c.seq)::bio_feature <~> ROW(t.emb, t.seq)::bio_feature) AS hybrid_score
                FROM candidates c, target t
                WHERE c.uniprot_id != t.uniprot_id
                  AND length(c.seq) BETWEEN t.len * 0.7 AND t.len * 1.3
                  AND sequence_entropy(c.seq) >= 2.0
                ORDER BY hybrid_score ASC
                LIMIT 5;
                """,
                (qid, qid, qid, oversample),
            )
            rows = cur.fetchall()
            dt_ms = (time.perf_counter() - t0) * 1000.0
            latencies.append(dt_ms)

            gt_h = ground_truth_hybrid_5.get(qid)
            if gt_h and gt_h["count"] > 0:
                hits = sum(
                    1 for r in rows if (r[0] in gt_h["ids"] or float(r[2]) <= gt_h["max_score"] + eps)
                )
                hybrid_recalls.append(hits / gt_h["count"])

        return {
            "mean_ms": round(statistics.mean(latencies), 2),
            "p50_ms": round(statistics.median(latencies), 2),
            "hybrid_recall_at_5": round(statistics.mean(hybrid_recalls) * 100.0, 2) if hybrid_recalls else 0.0,
        }

    results = {
        "sample_size": SAMPLE_SIZE,
        "num_queries": NUM_QUERIES,
        "storage_comparison": {
            "external_default": ext_storage,
            "plain_inline": plain_storage,
            "prism_external_build_s": round(prism_ext_build_s, 2),
            "prism_external_p50_ms": round(statistics.median(ext_lats), 2),
            "prism_external_mean_ms": round(statistics.mean(ext_lats), 2),
        },
        "brute_force": {
            "mean_ms": round(statistics.mean(exact_latencies_ms), 2),
            "p50_ms": round(statistics.median(exact_latencies_ms), 2),
        },
        "indexes": {},
        "query_sweeps": [],
        "hybrid_pg_bio": {},
    }

    # ------------------------------------------------------------------
    # Index 1: pg_vectorsearch PRISM (Default: SOAR + Boundary Replication)
    # ------------------------------------------------------------------
    print("\n[Phase 3A] Building pg_vectorsearch PRISM Index (Default: SOAR + Boundary + RaBitQ Fastscan)...")
    t0 = time.perf_counter()
    cur.execute(
        """
        CREATE INDEX idx_bench_prism_default
        ON bench_proteins_plain
        USING prism (embedding vec16_cosine_ops);
        """
    )
    prism_build_s = time.perf_counter() - t0
    prism_size_mb = round(get_size_bytes(cur, "idx_bench_prism_default") / (1024 * 1024), 2)

    cur.execute("SELECT name, setting, source FROM prism_index_settings('idx_bench_prism_default');")
    prism_settings = {r[0]: {"setting": r[1], "source": r[2]} for r in cur.fetchall()}
    print(f"  • Build Time: {prism_build_s:.2f} s | Index Size: {prism_size_mb} MB")
    print(f"  • Auto Settings: nlist={prism_settings.get('nlist', {}).get('setting')}, "
          f"nlevels={prism_settings.get('nlevels', {}).get('setting')}, "
          f"auto_nprobe={prism_settings.get('nprobe', {}).get('setting')}")

    results["indexes"]["prism_default"] = {
        "build_time_s": round(prism_build_s, 2),
        "size_mb": prism_size_mb,
        "settings": prism_settings,
    }

    for nprobe in [2, 5, 0, 20, 40, 60]:
        label = f"PRISM (default, nprobe={'auto=' + prism_settings.get('nprobe', {}).get('setting', '0') if nprobe == 0 else nprobe})"
        metrics = evaluate_queries(f"SET prism.nprobe = {nprobe};")
        metrics["index"] = "prism_default"
        metrics["param"] = f"nprobe={nprobe}"
        metrics["label"] = label
        results["query_sweeps"].append(metrics)
        print(f"    -> {label:36s}: p50={metrics['p50_ms']:5.2f} ms | mean={metrics['mean_ms']:5.2f} ms | QPS={metrics['qps']:6.1f} | Recall@10={metrics['recall_at_10']:6.2f}% | Recall@50={metrics['recall_at_50']:6.2f}%")

    for nprobe, key in [(5, "prism_nprobe5"), (0, "prism_auto"), (30, "prism_nprobe30")]:
        res_h = evaluate_hybrid_pg_bio_pipeline(f"SET prism.nprobe = {nprobe};")
        results["hybrid_pg_bio"][key] = res_h
        print(f"    -> Combined pg_bio <~> Hybrid Re-rank ({key:14s}): p50={res_h['p50_ms']:5.2f} ms | Hybrid Recall@5={res_h['hybrid_recall_at_5']:6.2f}%")

    cur.execute("DROP INDEX idx_bench_prism_default;")

    # ------------------------------------------------------------------
    # Index 2: pg_vectorsearch PRISM (Compact: No Replication, 1x postings)
    # ------------------------------------------------------------------
    print("\n[Phase 3B] Building pg_vectorsearch PRISM Index (Compact: soar_lambda=0, boundary_epsilon=0)...")
    t0 = time.perf_counter()
    cur.execute(
        """
        CREATE INDEX idx_bench_prism_compact
        ON bench_proteins_plain
        USING prism (embedding vec16_cosine_ops)
        WITH (soar_lambda = 0, boundary_epsilon = 0);
        """
    )
    prism_compact_build_s = time.perf_counter() - t0
    prism_compact_size_mb = round(get_size_bytes(cur, "idx_bench_prism_compact") / (1024 * 1024), 2)
    print(f"  • Build Time: {prism_compact_build_s:.2f} s | Index Size: {prism_compact_size_mb} MB")

    results["indexes"]["prism_compact"] = {
        "build_time_s": round(prism_compact_build_s, 2),
        "size_mb": prism_compact_size_mb,
    }

    for nprobe in [5, 0, 20, 40]:
        label = f"PRISM (compact, nprobe={'auto' if nprobe == 0 else nprobe})"
        metrics = evaluate_queries(f"SET prism.nprobe = {nprobe};")
        metrics["index"] = "prism_compact"
        metrics["param"] = f"nprobe={nprobe}"
        metrics["label"] = label
        results["query_sweeps"].append(metrics)
        print(f"    -> {label:36s}: p50={metrics['p50_ms']:5.2f} ms | mean={metrics['mean_ms']:5.2f} ms | QPS={metrics['qps']:6.1f} | Recall@10={metrics['recall_at_10']:6.2f}% | Recall@50={metrics['recall_at_50']:6.2f}%")

    cur.execute("DROP INDEX idx_bench_prism_compact;")

    # ------------------------------------------------------------------
    # Index 3: pgvector IVFFlat (lists = 316)
    # ------------------------------------------------------------------
    print("\n[Phase 3C] Building pgvector IVFFlat Index (lists = 316)...")
    t0 = time.perf_counter()
    cur.execute(
        """
        CREATE INDEX idx_bench_ivfflat
        ON bench_proteins_plain
        USING ivfflat (embedding halfvec_cosine_ops)
        WITH (lists = 316);
        """
    )
    ivf_build_s = time.perf_counter() - t0
    ivf_size_mb = round(get_size_bytes(cur, "idx_bench_ivfflat") / (1024 * 1024), 2)
    print(f"  • Build Time: {ivf_build_s:.2f} s | Index Size: {ivf_size_mb} MB")

    results["indexes"]["ivfflat"] = {
        "build_time_s": round(ivf_build_s, 2),
        "size_mb": ivf_size_mb,
    }

    for probes in [1, 5, 10, 20, 40]:
        label = f"pgvector IVFFlat (probes={probes})"
        metrics = evaluate_queries(f"SET ivfflat.probes = {probes};")
        metrics["index"] = "ivfflat"
        metrics["param"] = f"probes={probes}"
        metrics["label"] = label
        results["query_sweeps"].append(metrics)
        print(f"    -> {label:36s}: p50={metrics['p50_ms']:5.2f} ms | mean={metrics['mean_ms']:5.2f} ms | QPS={metrics['qps']:6.1f} | Recall@10={metrics['recall_at_10']:6.2f}% | Recall@50={metrics['recall_at_50']:6.2f}%")

    res_h = evaluate_hybrid_pg_bio_pipeline("SET ivfflat.probes = 10;")
    results["hybrid_pg_bio"]["ivfflat_probes10"] = res_h
    print(f"    -> Combined pg_bio <~> Hybrid Re-rank (ivfflat_probes10): p50={res_h['p50_ms']:5.2f} ms | Hybrid Recall@5={res_h['hybrid_recall_at_5']:6.2f}%")

    cur.execute("DROP INDEX idx_bench_ivfflat;")

    # ------------------------------------------------------------------
    # Index 4: pgvector HNSW (m = 16, ef_construction = 64)
    # ------------------------------------------------------------------
    print("\n[Phase 3D] Building pgvector HNSW Index (m = 16, ef_construction = 64)...")
    t0 = time.perf_counter()
    cur.execute(
        """
        CREATE INDEX idx_bench_hnsw
        ON bench_proteins_plain
        USING hnsw (embedding halfvec_cosine_ops)
        WITH (m = 16, ef_construction = 64);
        """
    )
    hnsw_build_s = time.perf_counter() - t0
    hnsw_size_mb = round(get_size_bytes(cur, "idx_bench_hnsw") / (1024 * 1024), 2)
    print(f"  • Build Time: {hnsw_build_s:.2f} s | Index Size: {hnsw_size_mb} MB")

    results["indexes"]["hnsw"] = {
        "build_time_s": round(hnsw_build_s, 2),
        "size_mb": hnsw_size_mb,
    }

    for ef_search in [10, 40, 100, 200]:
        label = f"pgvector HNSW (ef_search={ef_search})"
        metrics = evaluate_queries(f"SET hnsw.ef_search = {ef_search};")
        metrics["index"] = "hnsw"
        metrics["param"] = f"ef_search={ef_search}"
        metrics["label"] = label
        results["query_sweeps"].append(metrics)
        print(f"    -> {label:36s}: p50={metrics['p50_ms']:5.2f} ms | mean={metrics['mean_ms']:5.2f} ms | QPS={metrics['qps']:6.1f} | Recall@10={metrics['recall_at_10']:6.2f}% | Recall@50={metrics['recall_at_50']:6.2f}%")

    for ef, key in [(40, "hnsw_ef40"), (100, "hnsw_ef100"), (200, "hnsw_ef200")]:
        res_h = evaluate_hybrid_pg_bio_pipeline(f"SET hnsw.ef_search = {ef};")
        results["hybrid_pg_bio"][key] = res_h
        print(f"    -> Combined pg_bio <~> Hybrid Re-rank ({key:14s}): p50={res_h['p50_ms']:5.2f} ms | Hybrid Recall@5={res_h['hybrid_recall_at_5']:6.2f}%")

    # Clean up benchmark table
    cur.execute("DROP TABLE IF EXISTS bench_proteins_plain CASCADE;")

    out_path = os.path.join(os.path.dirname(__file__), "benchmark_pg_vectorsearch_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✅ Saved full benchmark results to {out_path}")

    cur.close()
    conn.close()


if __name__ == "__main__":
    run_benchmark()
