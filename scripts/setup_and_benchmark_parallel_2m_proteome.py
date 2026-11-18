# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "numpy",
# ]
# ///
"""
Parallel Full-Proteome Setup & Head-to-Head Benchmark (2,340,766 UniProt Proteins):
- Table 1 (Original, untouched): `proteins` (STORAGE EXTERNAL + pgvector HNSW `idx_protein_embedding`)
- Table 2 (Parallel copy):       `proteins_prism` (STORAGE PLAIN + pg_vectorsearch PRISM `idx_proteins_prism_embedding`)
"""

import json
import os
import statistics
import time
import psycopg

DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
NUM_QUERIES = int(os.getenv("BENCH_NUM_QUERIES", "30"))


def get_table_storage_breakdown(cur, table_name: str, index_name: str) -> dict:
    cur.execute(
        """
        SELECT
            pg_relation_size(c.oid, 'main') AS main_bytes,
            COALESCE(pg_relation_size(c.reltoastrelid), 0) AS toast_bytes,
            pg_relation_size(%s::regclass) AS vec_index_bytes,
            pg_indexes_size(c.oid) AS all_indexes_bytes,
            pg_total_relation_size(c.oid) AS total_bytes,
            c.reltuples::bigint AS est_rows
        FROM pg_class c
        WHERE c.relname = %s;
        """,
        (index_name, table_name),
    )
    row = cur.fetchone()
    mb = 1024 * 1024
    return {
        "table_name": table_name,
        "index_name": index_name,
        "main_heap_mb": round(row[0] / mb, 2),
        "toast_heap_mb": round(row[1] / mb, 2),
        "vector_index_mb": round(row[2] / mb, 2),
        "all_indexes_mb": round(row[3] / mb, 2),
        "total_footprint_mb": round(row[4] / mb, 2),
        "est_rows": row[5],
    }


def extract_buffers_from_explain(explain_rows: list[str]) -> dict:
    shared_hit = 0
    shared_read = 0
    exec_ms = 0.0
    for line in explain_rows:
        if "Buffers: shared" in line:
            parts = line.split("Buffers: shared")[1].strip().split()
            for p in parts:
                if p.startswith("hit="):
                    shared_hit += int(p.split("=")[1])
                elif p.startswith("read="):
                    shared_read += int(p.split("=")[1])
        if "Execution Time:" in line:
            exec_ms = float(line.split("Execution Time:")[1].split("ms")[0].strip())
    return {"shared_hit": shared_hit, "shared_read": shared_read, "exec_ms": exec_ms}


def main():
    print("=============================================================================")
    print("🧬 Full 2.34M-Proteome Parallel Copy & Benchmark: `proteins` vs `proteins_prism`")
    print("=============================================================================\n")

    conn = psycopg.connect(DB_URL, autocommit=True)
    cur = conn.cursor()

    cur.execute("CREATE EXTENSION IF NOT EXISTS vector;")
    cur.execute("CREATE EXTENSION IF NOT EXISTS pg_vectorsearch;")
    cur.execute("CREATE EXTENSION IF NOT EXISTS pg_bio;")
    cur.execute("SET maintenance_work_mem = '4GB';")
    cur.execute("SET max_parallel_maintenance_workers = 6;")

    # Check if proteins_prism already exists with full row count
    cur.execute("SELECT to_regclass('public.proteins_prism');")
    exists = cur.fetchone()[0] is not None
    copy_time_s = 0.0
    pkey_time_s = 0.0

    cur.execute("SELECT count(*) FROM proteins;")
    orig_count = cur.fetchone()[0]
    print(f"• Original `proteins` table row count: {orig_count:,} (kept 100% untouched)")

    need_copy = True
    if exists:
        cur.execute("SELECT count(*) FROM proteins_prism;")
        prism_count = cur.fetchone()[0]
        if prism_count == orig_count:
            print(f"• Parallel table `proteins_prism` already has {prism_count:,} rows.")
            need_copy = False

    if need_copy:
        print("\n[Step 1] Creating parallel table `proteins_prism` (with `embedding halfvec(1280) STORAGE PLAIN`)...")
        cur.execute("DROP TABLE IF EXISTS proteins_prism CASCADE;")
        cur.execute(
            """
            CREATE TABLE proteins_prism (
                uniprot_id character varying(20) NOT NULL,
                name text,
                sequence text,
                embedding halfvec(1280) STORAGE PLAIN
            );
            """
        )
        t0 = time.perf_counter()
        cur.execute(
            """
            INSERT INTO proteins_prism (uniprot_id, name, sequence, embedding)
            SELECT uniprot_id, name, sequence, embedding
            FROM proteins;
            """
        )
        copy_time_s = time.perf_counter() - t0
        print(f"  -> Copied {orig_count:,} proteins in {copy_time_s:.2f} s")

        print("[Step 2] Building PRIMARY KEY (uniprot_id) and running ANALYZE on `proteins_prism`...")
        t0 = time.perf_counter()
        cur.execute("ALTER TABLE proteins_prism ADD CONSTRAINT proteins_prism_pkey PRIMARY KEY (uniprot_id);")
        cur.execute("ANALYZE proteins_prism;")
        pkey_time_s = time.perf_counter() - t0
        print(f"  -> Built primary key & analyzed in {pkey_time_s:.2f} s")

    # Check if PRISM index exists
    cur.execute("SELECT to_regclass('public.idx_proteins_prism_embedding');")
    idx_exists = cur.fetchone()[0] is not None
    prism_build_s = 0.0

    if not idx_exists:
        print("\n[Step 3] Building `idx_proteins_prism_embedding` USING prism (embedding vec16_cosine_ops) on 2.34M rows...")
        t0 = time.perf_counter()
        cur.execute(
            """
            CREATE INDEX idx_proteins_prism_embedding
            ON proteins_prism
            USING prism (embedding vec16_cosine_ops);
            """
        )
        prism_build_s = time.perf_counter() - t0
        print(f"  -> Built PRISM index over {orig_count:,} 1280-D proteins in {prism_build_s:.2f} s!")

    # Inspect resolved PRISM settings
    cur.execute("SELECT name, setting, source FROM prism_index_settings('idx_proteins_prism_embedding');")
    prism_settings = {r[0]: {"setting": r[1], "source": r[2]} for r in cur.fetchall()}
    print(
        f"  -> PRISM Index Settings: nlist={prism_settings.get('nlist', {}).get('setting')}, "
        f"fan_out={prism_settings.get('fan_out', {}).get('setting')}, "
        f"nlevels={prism_settings.get('nlevels', {}).get('setting')}, "
        f"auto_nprobe={prism_settings.get('nprobe', {}).get('setting')}"
    )

    # Compare full-scale storage breakdown
    orig_storage = get_table_storage_breakdown(cur, "proteins", "idx_protein_embedding")
    prism_storage = get_table_storage_breakdown(cur, "proteins_prism", "idx_proteins_prism_embedding")

    print("\n[Step 4] Full 2.34M-Proteome Storage Comparison:")
    print(
        f"  • `proteins` (HNSW + EXTERNAL):      Main={orig_storage['main_heap_mb']:7.2f} MB | "
        f"TOAST={orig_storage['toast_heap_mb']:7.2f} MB | HNSW Index={orig_storage['vector_index_mb']:7.2f} MB | "
        f"Total={orig_storage['total_footprint_mb']:7.2f} MB"
    )
    print(
        f"  • `proteins_prism` (PRISM + PLAIN):  Main={prism_storage['main_heap_mb']:7.2f} MB | "
        f"TOAST={prism_storage['toast_heap_mb']:7.2f} MB | PRISM Index={prism_storage['vector_index_mb']:7.2f} MB | "
        f"Total={prism_storage['total_footprint_mb']:7.2f} MB"
    )
    idx_ratio = orig_storage["vector_index_mb"] / max(prism_storage["vector_index_mb"], 0.01)
    total_saved_mb = orig_storage["total_footprint_mb"] - prism_storage["total_footprint_mb"]
    print(f"  => Vector Index Reduction: {idx_ratio:.2f}x smaller ({orig_storage['vector_index_mb']} MB -> {prism_storage['vector_index_mb']} MB)")
    print(f"  => Total Table+Index Disk Saved: {total_saved_mb:.2f} MB ({total_saved_mb/1024:.2f} GB saved!)")

    # Select deterministic sample of query proteins
    print(f"\n[Step 5] Selecting {NUM_QUERIES} query proteins for parallel head-to-head benchmarking...")
    cur.execute(
        f"""
        SELECT uniprot_id
        FROM proteins_prism
        WHERE length(sequence) BETWEEN 150 AND 450
        LIMIT {NUM_QUERIES * 20};
        """
    )
    pool = [r[0] for r in cur.fetchall()]
    query_ids = pool[::20][:NUM_QUERIES]

    def bench_table(table_name: str, guc_sql: str, k: int = 10):
        if guc_sql:
            cur.execute(guc_sql)
        # Warmup
        for qid in query_ids[:3]:
            cur.execute(
                f"""
                SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s)) AS dist
                FROM {table_name} p
                ORDER BY p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s) ASC
                LIMIT %s;
                """,
                (qid, qid, k + 1),
            )
            cur.fetchall()

        latencies = []
        results_by_qid = {}
        total_mean_top10_dist = []
        shared_hits = []
        shared_reads = []

        for qid in query_ids:
            t0 = time.perf_counter()
            cur.execute(
                f"""
                SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s)) AS dist
                FROM {table_name} p
                ORDER BY p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s) ASC
                LIMIT %s;
                """,
                (qid, qid, k + 1),
            )
            rows = [r for r in cur.fetchall() if r[0] != qid][:k]
            dt_ms = (time.perf_counter() - t0) * 1000.0
            latencies.append(dt_ms)
            results_by_qid[qid] = rows
            if rows:
                total_mean_top10_dist.append(statistics.mean(float(r[1]) for r in rows[:10]))

        # Also collect EXPLAIN (ANALYZE, BUFFERS) on first 5 queries
        for qid in query_ids[:5]:
            cur.execute(
                f"""
                EXPLAIN (ANALYZE, BUFFERS)
                SELECT p.uniprot_id, (p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s)) AS dist
                FROM {table_name} p
                ORDER BY p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s) ASC
                LIMIT %s;
                """,
                (qid, qid, k + 1),
            )
            exp_lines = [r[0] for r in cur.fetchall()]
            binfo = extract_buffers_from_explain(exp_lines)
            shared_hits.append(binfo["shared_hit"])
            shared_reads.append(binfo["shared_read"])

        lats_sorted = sorted(latencies)
        p95_idx = min(int(len(lats_sorted) * 0.95), len(lats_sorted) - 1)
        return {
            "p50_ms": round(statistics.median(latencies), 2),
            "mean_ms": round(statistics.mean(latencies), 2),
            "p95_ms": round(lats_sorted[p95_idx], 2),
            "qps": round(1000.0 / statistics.mean(latencies), 1),
            "avg_top10_cosine_dist": round(statistics.mean(total_mean_top10_dist), 6),
            "avg_buffer_blocks_hit_or_read": round(statistics.mean(h + r for h, r in zip(shared_hits, shared_reads)), 1),
            "results_by_qid": results_by_qid,
        }

    print("\n[Step 6] Running Parallel Queries Across 2.34M Proteins (`proteins` [HNSW] vs `proteins_prism` [PRISM])...")

    # Benchmark HNSW at ef_search = 40 (default), 100, 200
    hnsw_runs = {}
    for ef in [40, 100, 200]:
        res = bench_table("proteins", f"SET hnsw.ef_search = {ef};", k=10)
        hnsw_runs[f"ef_{ef}"] = res
        print(
            f"  • HNSW  (ef_search={ef:3d})  on `proteins`:       p50={res['p50_ms']:6.2f} ms | "
            f"mean={res['mean_ms']:6.2f} ms | p95={res['p95_ms']:6.2f} ms | QPS={res['qps']:6.1f} | "
            f"avg_dist={res['avg_top10_cosine_dist']:.5f} | buffers/query={res['avg_buffer_blocks_hit_or_read']}"
        )

    # Benchmark PRISM at nprobe = 10, 0 (auto), 40, 80
    prism_runs = {}
    for nprobe in [10, 0, 40, 80]:
        label = f"auto={prism_settings.get('nprobe', {}).get('setting', '0')}" if nprobe == 0 else str(nprobe)
        res = bench_table("proteins_prism", f"SET prism.nprobe = {nprobe};", k=10)
        prism_runs[f"nprobe_{nprobe}"] = res
        print(
            f"  • PRISM (nprobe={label:7s}) on `proteins_prism`: p50={res['p50_ms']:6.2f} ms | "
            f"mean={res['mean_ms']:6.2f} ms | p95={res['p95_ms']:6.2f} ms | QPS={res['qps']:6.1f} | "
            f"avg_dist={res['avg_top10_cosine_dist']:.5f} | buffers/query={res['avg_buffer_blocks_hit_or_read']}"
        )

    # Compute agreement / distance-tie recall between HNSW (ef=200) and PRISM (nprobe=80)
    eps = 1e-5
    overlaps = []
    prism_better_or_equal = 0
    total_pairs = 0
    for qid in query_ids:
        h_rows = hnsw_runs["ef_200"]["results_by_qid"][qid]
        p_rows = prism_runs["nprobe_80"]["results_by_qid"][qid]
        if not h_rows or not p_rows:
            continue
        h_ids = set(r[0] for r in h_rows)
        h_max_dist = float(h_rows[-1][1])
        hits = sum(1 for r in p_rows if (r[0] in h_ids or float(r[1]) <= h_max_dist + eps))
        overlaps.append(hits / len(h_rows))
        for hr, pr in zip(h_rows, p_rows):
            total_pairs += 1
            if float(pr[1]) <= float(hr[1]) + eps:
                prism_better_or_equal += 1

    agreement_pct = round(statistics.mean(overlaps) * 100.0, 2)
    print(f"\n  => Top-10 Distance-Tie Agreement (PRISM nprobe=80 vs HNSW ef_search=200): {agreement_pct:.2f}%")

    # Benchmark Combined pg_bio Hybrid <~> Smith-Waterman + Cosine Re-ranking on 2.34M rows
    print("\n[Step 7] Benchmarking Combined `pg_bio` `<~>` Hybrid Smith-Waterman Re-Ranking on 2.34M Proteins...")

    def bench_hybrid_2m(table_name: str, guc_sql: str, oversample: int = 50):
        if guc_sql:
            cur.execute(guc_sql)
        lats = []
        scores = []
        for qid in query_ids[:15]:
            t0 = time.perf_counter()
            cur.execute(
                f"""
                WITH target AS (
                    SELECT uniprot_id, embedding AS emb_half, embedding::vector AS emb_vec, sequence AS seq, length(sequence)::float8 AS len
                    FROM {table_name}
                    WHERE uniprot_id = %s
                ),
                candidates AS (
                    SELECT p.uniprot_id, p.name, p.embedding::vector AS emb_vec, p.sequence AS seq,
                           (p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s)) AS cos_dist
                    FROM {table_name} p
                    ORDER BY p.embedding <=> (SELECT embedding FROM {table_name} WHERE uniprot_id = %s) ASC
                    LIMIT %s
                )
                SELECT c.uniprot_id,
                       c.cos_dist,
                       (ROW(c.emb_vec, c.seq)::bio_feature <~> ROW(t.emb_vec, t.seq)::bio_feature) AS hybrid_score
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
            lats.append((time.perf_counter() - t0) * 1000.0)
            if rows:
                scores.append(float(rows[0][2]))
        return {
            "p50_ms": round(statistics.median(lats), 2),
            "mean_ms": round(statistics.mean(lats), 2),
            "avg_top1_hybrid_score": round(statistics.mean(scores), 5) if scores else 0.0,
        }

    hybrid_hnsw_40 = bench_hybrid_2m("proteins", "SET hnsw.ef_search = 40;")
    hybrid_hnsw_100 = bench_hybrid_2m("proteins", "SET hnsw.ef_search = 100;")
    hybrid_prism_auto = bench_hybrid_2m("proteins_prism", "SET prism.nprobe = 0;")
    hybrid_prism_40 = bench_hybrid_2m("proteins_prism", "SET prism.nprobe = 40;")

    print(f"  • `proteins`       (HNSW ef=40   + pg_bio <~>): p50={hybrid_hnsw_40['p50_ms']:6.2f} ms | mean={hybrid_hnsw_40['mean_ms']:6.2f} ms | top1_hybrid={hybrid_hnsw_40['avg_top1_hybrid_score']:.5f}")
    print(f"  • `proteins`       (HNSW ef=100  + pg_bio <~>): p50={hybrid_hnsw_100['p50_ms']:6.2f} ms | mean={hybrid_hnsw_100['mean_ms']:6.2f} ms | top1_hybrid={hybrid_hnsw_100['avg_top1_hybrid_score']:.5f}")
    print(f"  • `proteins_prism` (PRISM auto   + pg_bio <~>): p50={hybrid_prism_auto['p50_ms']:6.2f} ms | mean={hybrid_prism_auto['mean_ms']:6.2f} ms | top1_hybrid={hybrid_prism_auto['avg_top1_hybrid_score']:.5f}")
    print(f"  • `proteins_prism` (PRISM np=40  + pg_bio <~>): p50={hybrid_prism_40['p50_ms']:6.2f} ms | mean={hybrid_prism_40['mean_ms']:6.2f} ms | top1_hybrid={hybrid_prism_40['avg_top1_hybrid_score']:.5f}")

    # Strip raw per-query rows before saving JSON summary
    for d in hnsw_runs.values():
        d.pop("results_by_qid", None)
    for d in prism_runs.values():
        d.pop("results_by_qid", None)

    summary = {
        "row_count": orig_count,
        "copy_time_s": round(copy_time_s, 2),
        "pkey_analyze_time_s": round(pkey_time_s, 2),
        "prism_build_time_s": round(prism_build_s, 2),
        "prism_settings": prism_settings,
        "storage": {
            "proteins_hnsw_external": orig_storage,
            "proteins_prism_plain": prism_storage,
            "index_size_reduction_factor": round(idx_ratio, 2),
            "total_disk_saved_mb": round(total_saved_mb, 2),
        },
        "ann_top10_runs": {
            "proteins_hnsw": hnsw_runs,
            "proteins_prism": prism_runs,
            "top10_agreement_pct": agreement_pct,
        },
        "hybrid_pg_bio_2m": {
            "hnsw_ef40": hybrid_hnsw_40,
            "hnsw_ef100": hybrid_hnsw_100,
            "prism_auto": hybrid_prism_auto,
            "prism_nprobe40": hybrid_prism_40,
        },
    }

    out_path = os.path.join(os.path.dirname(__file__), "benchmark_parallel_2m_proteome_results.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n✅ Saved 2.34M-proteome parallel benchmark summary to {out_path}")

    cur.close()
    conn.close()


if __name__ == "__main__":
    main()
