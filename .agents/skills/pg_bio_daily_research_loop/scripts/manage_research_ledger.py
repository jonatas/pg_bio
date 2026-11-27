#!/usr/bin/env python3
"""
Local Storage & Dual-Ledger Manager for pg_bio Autonomous Daily Research Loop.
Persists research attempts, worktree prototypes, benchmarks, PRs, blog posts,
and video backlog items in both:
  1. PostgreSQL core research database (`postgresql://jonatas@localhost:28818/bio_demo`)
  2. Git-tracked JSON snapshot (`/Users/jonatas/code/pg_bio/data/research_ledger.json`)
"""

import argparse
import datetime
import json
import os
import sys
from pathlib import Path
import psycopg
from psycopg.rows import dict_row

DB_URI = os.getenv("PG_BIO_URL", "postgresql://jonatas@localhost:28818/bio_demo")
PG_BIO_ROOT = Path("/Users/jonatas/code/pg_bio")
LEDGER_JSON_PATH = PG_BIO_ROOT / "data" / "research_ledger.json"

DDL_SQL = """
CREATE TABLE IF NOT EXISTS public.pg_bio_research_attempts (
    attempt_id SERIAL PRIMARY KEY,
    attempt_code VARCHAR(64) UNIQUE NOT NULL,
    run_date DATE NOT NULL DEFAULT CURRENT_DATE,
    category VARCHAR(40) NOT NULL,
    paper_title TEXT NOT NULL,
    paper_ref TEXT NOT NULL,
    pg_strategy TEXT NOT NULL,
    worktree_path TEXT,
    branch_name VARCHAR(120),
    status VARCHAR(32) NOT NULL DEFAULT 'Proposed',
    theory_test_name TEXT,
    benchmark_summary JSONB NOT NULL DEFAULT '{}'::jsonb,
    pr_number INTEGER,
    pr_url TEXT,
    blog_post_slug TEXT,
    video_backlog_id INTEGER,
    notes TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT now(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS public.pg_bio_video_backlog (
    video_id SERIAL PRIMARY KEY,
    episode_slug VARCHAR(100) UNIQUE NOT NULL,
    title TEXT NOT NULL,
    layer_name VARCHAR(80) NOT NULL,
    attempt_code VARCHAR(64) REFERENCES public.pg_bio_research_attempts(attempt_code) ON DELETE SET NULL,
    pr_number INTEGER,
    blog_post_slug TEXT,
    bio_problem TEXT NOT NULL,
    cs_mechanism TEXT NOT NULL,
    speedup_metric TEXT NOT NULL,
    status VARCHAR(32) NOT NULL DEFAULT 'Queued',
    created_at TIMESTAMPTZ NOT NULL DEFAULT now(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT now()
);
"""

SEED_ATTEMPTS = [
    {
        "attempt_code": "ATT-2026-11-18-PRISM-RABITQ",
        "run_date": "2026-11-18",
        "category": "Postgres_Systems",
        "paper_title": "RaBitQ 1-Bit Quantization & SOAR Hierarchical Centroid Routing (TigerData pg_vectorsearch PRISM)",
        "paper_ref": "https://github.com/timescale/pg_vectorsearch | PR #13",
        "pg_strategy": "USING prism (embedding vec16_cosine_ops) with 1-bit RaBitQ error bounds and SOAR spill-orthogonal cluster replication across 2.34M 1280-D UniProt proteins.",
        "worktree_path": ".worktrees/pg-vectorsearch-bench",
        "branch_name": "feat/pg-vectorsearch-benchmark",
        "status": "Published",
        "theory_test_name": "test_vector_search_prism_benchmark",
        "benchmark_summary": {
            "index_size_reduction": "9.72x (8,772 MB HNSW -> 902 MB PRISM)",
            "build_time": "3.08 min vs 27 hours",
            "query_speedup": "7.5x to 11.3x (26-37 ms vs 197-417 ms)"
        },
        "pr_number": 13,
        "pr_url": "https://github.com/jonatas/pg_bio/pull/13",
        "blog_post_slug": "from-27-hours-to-3-minutes-tigerdata-pg-vectorsearch-prism-pg-bio",
        "notes": "Merged in PR #13 and published on ideia.me."
    },
    {
        "attempt_code": "ATT-2026-11-21-LAYERED-CENTROIDS-FASTSCAN",
        "run_date": "2026-11-21",
        "category": "Postgres_Systems",
        "paper_title": "2-Level Layered Centroids, 16-Byte SIMD Register Fastscan (VPSHUFB / NEON vtbl), 3D Morton Octree & RangeFinder GiST",
        "paper_ref": "https://github.com/jonatas/pg_bio/pull/14",
        "pg_strategy": "100% lossless 2-level layered centroid routing (64 L1 x 512 L2), 16-byte SIMD register LUT fastscan, 3D Morton octree active-site alignment, and 128-byte ECFP4 + 9D RangeFinder GiST.",
        "worktree_path": "../pg_bio_vectorsearch",
        "branch_name": "feat/layered-centroids-multi-index",
        "status": "Published",
        "theory_test_name": "test_layered_centroids_and_fastscan",
        "benchmark_summary": {
            "fastscan_width": "32-wide SIMD register LUT (NEON vtbl / AVX2 vpshufb)",
            "lossless_recall": "100% recall on layered centroid partitioning"
        },
        "pr_number": 14,
        "pr_url": "https://github.com/jonatas/pg_bio/pull/14",
        "blog_post_slug": "layered-centroids-register-fastscan-morton-octree-ecfp4-rangefinder-pg-bio",
        "notes": "Merged in PR #14 and published on ideia.me."
    },
    {
        "attempt_code": "ATT-2026-11-23-KNOW-THY-NEIGHBOUR-ZCA",
        "run_date": "2026-11-23",
        "category": "PLM_Kinetics",
        "paper_title": "Know Thy Neighbour: Full Symmetric ZCA Whitening for Anisotropic PLM Embeddings & Inline 1,304B BioVecInt8G128 Varlena",
        "paper_ref": "gbouras13/know-thy-neighbour (SCOPe/ASTRAL-40 2.08) | PR #15",
        "pg_strategy": "Full symmetric ZCA whitening (bio_fit_zca_whiten, bio_whiten_embedding), 1,304-byte zero-TOAST BioVecInt8G128 group-wise int8 varlena, Kabsch structure sanitizer, PAINS screen, and Stage-2.5 EnzRef/STORM L0.",
        "worktree_path": ".worktrees/whitening-inline-varlena",
        "branch_name": "feat/whitening-inline-varlena-hygiene-rerankers",
        "status": "Published",
        "theory_test_name": "test_phase1_anisotropy_whitening_and_rabitq_eps_bound",
        "benchmark_summary": {
            "scope40_top1_superfamily": "58.55% raw -> 88.32% ZCA whitened (3,503 held-out family queries)",
            "rabitq_eps_bound": "2.8108 -> 1.0118 (2.78x tighter)",
            "varlena_bytes": "1,304 B (0 B TOAST, 2.00x faster SIMD i32 rerank)"
        },
        "pr_number": 15,
        "pr_url": "https://github.com/jonatas/pg_bio/pull/15",
        "blog_post_slug": "know-thy-neighbour-zca-whitening-inline-biovecint8g128-pains-enzref-storm-pg-bio",
        "notes": "Merged in PR #15 and published on ideia.me."
    },
    {
        "attempt_code": "ATT-2026-11-24-POCKETVINA-POSEBUSTERS",
        "run_date": "2026-11-24",
        "category": "Structural_Bio",
        "paper_title": "PocketVina Atlas Multi-Pocket SAS Conditioning (arXiv:2506.20043) & 21-Check PoseBusters Physical Validity Gate (Chem. Sci. 2024)",
        "paper_ref": "arXiv:2506.20043 | doi:10.1039/D3SC04185A | PR #16",
        "pg_strategy": "24-ray Fibonacci Shrake-Rupley SAS concavity detector with 5.0 A 3D voxel hash grid (bio_detect_top_pockets) + 21-check PoseBusters audit with Strict Non-Hiding Dual-Ledger salt/steric provenance (bio_posebusters_audit_row, bio_multipocket_posebusters_screen).",
        "worktree_path": ".worktrees/pocketvina-multipocket-posebusters",
        "branch_name": "feat/pocketvina-multipocket-posebusters",
        "status": "Published",
        "theory_test_name": "test_pocketvina_multipocket_detector_and_posebusters_21_gate",
        "benchmark_summary": {
            "latency_per_complex_ms": 1.65,
            "total_120_complexes_ms": 197.92,
            "speedup_vs_python": "20.96x faster than Python pull",
            "posebusters_checks": 21
        },
        "pr_number": 16,
        "pr_url": "https://github.com/jonatas/pg_bio/pull/16",
        "blog_post_slug": "pocketvina-posebusters-enzref-pslim-tethers-pg-bio",
        "notes": "Merged in PR #16; tying together with PR #19 in unified blog post and video backlog."
    },
    {
        "attempt_code": "ATT-2026-11-26-ENZREF-PSLIM",
        "run_date": "2026-11-26",
        "category": "PLM_Kinetics",
        "paper_title": "EnzRef Reference-Centered Kinetic Role-Contrast & Disordered Microprotein SLiM Window Pooler / Dual-SLiM Metabolon Tethers",
        "paper_ref": "EnzRef (Reference-Centered Enzyme-Substrate-Product Role Contrast) & pSLiM Disordered Microprotein Motif Indexing | PR #19",
        "pg_strategy": "Reaction SMILES agent-stripping splitter (bio_split_reaction_roles), exact-antisymmetric SIMD kinetic role-contrast & backpressure scorer over f32 and 1,304B BioVecInt8G128 (bio_enzref_kinetics, bio_enzref_flux_rank_g128), and 16-D ZCA-whitened unit-L2 SLiM window pooler + dual-SLiM metabolon tether matcher with hydrophobic-core linker junction sanity gate (bio_pslim_window_pool, bio_match_dual_slim_tether, bio_audit_fusion_junction_sanity).",
        "worktree_path": ".worktrees/enzref-pslim",
        "branch_name": "feat/enzref-role-contrast-pslim-tethers",
        "status": "Published",
        "theory_test_name": "test_phase3_enzref_role_contrast_and_storm_l0",
        "benchmark_summary": {
            "enzref_1280d_simd_throughput": "2,092,090 rows/sec/core (4.18x above 500k target)",
            "pslim_zca_l2_window_throughput": "617,845 windows/sec/core",
            "exact_antisymmetry": "bio_enzref_flux_rank(e, p, r) == -bio_enzref_flux_rank(e, r, p)"
        },
        "pr_number": 19,
        "pr_url": "https://github.com/jonatas/pg_bio/pull/19",
        "blog_post_slug": "pocketvina-posebusters-enzref-pslim-tethers-pg-bio",
        "notes": "Merged in PR #19; tying together with PR #16 in unified blog post and video backlog."
    }
]

SEED_VIDEOS = [
    {
        "episode_slug": "ep03_pocketvina_posebusters_21_gate",
        "title": "Multi-Pocket SAS Concavity & 21-Check PoseBusters Gate with Dual-Ledger Provenance",
        "layer_name": "Layer 19: Multi-Pocket Conditioning & PoseBusters Physical Gate",
        "attempt_code": "ATT-2026-11-24-POCKETVINA-POSEBUSTERS",
        "pr_number": 16,
        "blog_post_slug": "pocketvina-posebusters-enzref-pslim-tethers-pg-bio",
        "bio_problem": "Single-box docking misses allosteric/cryptic pockets while raw kcal/mol scores hallucinate steric clashes, warped aromatic rings, and bare salt counterions ([Na+], [Cl-]).",
        "cs_mechanism": "24-ray Fibonacci Shrake-Rupley SAS probes with O(1) 5.0 A 3D spatial voxel hashing, 9D RangeFinder GiST intervals, 16-byte register Fastscan, and 21-check PoseBusters Rust SVD/steric validator with Non-Hiding Dual-Ledger provenance.",
        "speedup_metric": "20.96x faster than Python/RDKit pull (1.65 ms/complex across 120 complexes)",
        "status": "Published_IdeiaMe"
    },
    {
        "episode_slug": "ep04_enzref_kinetics_and_pslim_tethers",
        "title": "SIMD EnzRef Catalytic Role-Contrast & ZCA-Whitened Disordered Microprotein SLiM Tethers",
        "layer_name": "Layer 20: EnzRef Reaction Role-Contrast & Disordered SLiM Metabolon Tethers",
        "attempt_code": "ATT-2026-11-26-ENZREF-PSLIM",
        "pr_number": 19,
        "blog_post_slug": "pocketvina-posebusters-enzref-pslim-tethers-pg-bio",
        "bio_problem": "Standard vector similarity cannot distinguish a dead non-catalytic binder from a high-kcat forward enzyme or detect when a synthetic chimera linker slices through a buried hydrophobic core.",
        "cs_mechanism": "Reference-centered SIMD role-contrast (d_react = r - e, c_fwd = p - r) over inline 1,304-byte BioVecInt8G128 varlenas with exact antisymmetry and backpressure risk, paired with 16-D ZCA-whitened unit-L2 SLiM window pooling and GRAVY/hydrophobic-run junction sanity checks.",
        "speedup_metric": "2.09M EnzRef 1280-D triplets/sec/core & 617.8k ZCA-whitened SLiM windows/sec/core",
        "status": "Published_IdeiaMe"
    }
]


def init_and_seed(conn):
    with conn.cursor() as cur:
        cur.execute(DDL_SQL)
        for att in SEED_ATTEMPTS:
            cur.execute(
                """
                INSERT INTO public.pg_bio_research_attempts (
                    attempt_code, run_date, category, paper_title, paper_ref,
                    pg_strategy, worktree_path, branch_name, status,
                    theory_test_name, benchmark_summary, pr_number, pr_url,
                    blog_post_slug, notes
                ) VALUES (
                    %(attempt_code)s, %(run_date)s, %(category)s, %(paper_title)s, %(paper_ref)s,
                    %(pg_strategy)s, %(worktree_path)s, %(branch_name)s, %(status)s,
                    %(theory_test_name)s, %(benchmark_summary)s::jsonb, %(pr_number)s, %(pr_url)s,
                    %(blog_post_slug)s, %(notes)s
                )
                ON CONFLICT (attempt_code) DO UPDATE SET
                    status = EXCLUDED.status,
                    benchmark_summary = EXCLUDED.benchmark_summary,
                    pr_number = COALESCE(EXCLUDED.pr_number, public.pg_bio_research_attempts.pr_number),
                    pr_url = COALESCE(EXCLUDED.pr_url, public.pg_bio_research_attempts.pr_url),
                    blog_post_slug = COALESCE(EXCLUDED.blog_post_slug, public.pg_bio_research_attempts.blog_post_slug),
                    notes = EXCLUDED.notes,
                    updated_at = now();
                """,
                {
                    **att,
                    "benchmark_summary": json.dumps(att["benchmark_summary"]),
                },
            )
        for vid in SEED_VIDEOS:
            cur.execute(
                """
                INSERT INTO public.pg_bio_video_backlog (
                    episode_slug, title, layer_name, attempt_code, pr_number,
                    blog_post_slug, bio_problem, cs_mechanism, speedup_metric, status
                ) VALUES (
                    %(episode_slug)s, %(title)s, %(layer_name)s, %(attempt_code)s, %(pr_number)s,
                    %(blog_post_slug)s, %(bio_problem)s, %(cs_mechanism)s, %(speedup_metric)s, %(status)s
                )
                ON CONFLICT (episode_slug) DO UPDATE SET
                    blog_post_slug = COALESCE(EXCLUDED.blog_post_slug, public.pg_bio_video_backlog.blog_post_slug),
                    speedup_metric = EXCLUDED.speedup_metric,
                    updated_at = now()
                RETURNING video_id, attempt_code;
                """,
                vid,
            )
            row = cur.fetchone()
            if row and row[1]:
                cur.execute(
                    """
                    UPDATE public.pg_bio_research_attempts
                    SET video_backlog_id = %s, updated_at = now()
                    WHERE attempt_code = %s;
                    """,
                    (row[0], row[1]),
                )
    conn.commit()


def export_json_snapshot(conn):
    with conn.cursor(row_factory=dict_row) as cur:
        cur.execute(
            """
            SELECT attempt_id, attempt_code, run_date::text AS run_date, category,
                   paper_title, paper_ref, pg_strategy, worktree_path, branch_name,
                   status, theory_test_name, benchmark_summary, pr_number, pr_url,
                   blog_post_slug, video_backlog_id, notes,
                   created_at::text AS created_at, updated_at::text AS updated_at
            FROM public.pg_bio_research_attempts
            ORDER BY run_date DESC, attempt_id DESC;
            """
        )
        attempts = cur.fetchall()

        cur.execute(
            """
            SELECT video_id, episode_slug, title, layer_name, attempt_code,
                   pr_number, blog_post_slug, bio_problem, cs_mechanism,
                   speedup_metric, status,
                   created_at::text AS created_at, updated_at::text AS updated_at
            FROM public.pg_bio_video_backlog
            ORDER BY video_id ASC;
            """
        )
        videos = cur.fetchall()

    payload = {
        "schema_version": "1.0",
        "last_synced_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "database_uri": "postgresql://localhost:28818/bio_demo",
        "research_attempts": attempts,
        "video_backlog": videos,
    }
    LEDGER_JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    LEDGER_JSON_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def cmd_status(conn):
    payload = export_json_snapshot(conn)
    print("=== pg_bio Autonomous Research Ledger (SQL + JSON Synced) ===")
    print(f"Total Research Attempts: {len(payload['research_attempts'])}")
    for a in payload["research_attempts"]:
        pr_str = f"PR #{a['pr_number']}" if a["pr_number"] else "No PR yet"
        blog_str = a["blog_post_slug"] or "pending-blog"
        print(
            f"  • [{a['status']}] {a['attempt_code']} ({a['run_date']}) | {pr_str} | {blog_str}\n"
            f"    Title: {a['paper_title']}\n"
            f"    Bench: {json.dumps(a['benchmark_summary'], ensure_ascii=False)}"
        )
    print(f"\nTotal Video Backlog Tasks: {len(payload['video_backlog'])}")
    for v in payload["video_backlog"]:
        print(
            f"  • [{v['status']}] #{v['video_id']} {v['episode_slug']} (PR #{v['pr_number']}) -> {v['title']}\n"
            f"    Speedup: {v['speedup_metric']}"
        )


def cmd_upsert_attempt(conn, args):
    bench = json.loads(args.benchmark_json) if args.benchmark_json else {}
    with conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO public.pg_bio_research_attempts (
                attempt_code, run_date, category, paper_title, paper_ref,
                pg_strategy, worktree_path, branch_name, status,
                theory_test_name, benchmark_summary, pr_number, pr_url,
                blog_post_slug, notes
            ) VALUES (
                %s, COALESCE(%s::date, CURRENT_DATE), %s, %s, %s,
                %s, %s, %s, %s,
                %s, %s::jsonb, %s, %s,
                %s, %s
            )
            ON CONFLICT (attempt_code) DO UPDATE SET
                category = EXCLUDED.category,
                paper_title = EXCLUDED.paper_title,
                paper_ref = EXCLUDED.paper_ref,
                pg_strategy = EXCLUDED.pg_strategy,
                worktree_path = COALESCE(EXCLUDED.worktree_path, public.pg_bio_research_attempts.worktree_path),
                branch_name = COALESCE(EXCLUDED.branch_name, public.pg_bio_research_attempts.branch_name),
                status = EXCLUDED.status,
                theory_test_name = COALESCE(EXCLUDED.theory_test_name, public.pg_bio_research_attempts.theory_test_name),
                benchmark_summary = CASE
                    WHEN EXCLUDED.benchmark_summary = '{}'::jsonb THEN public.pg_bio_research_attempts.benchmark_summary
                    ELSE EXCLUDED.benchmark_summary
                END,
                pr_number = COALESCE(EXCLUDED.pr_number, public.pg_bio_research_attempts.pr_number),
                pr_url = COALESCE(EXCLUDED.pr_url, public.pg_bio_research_attempts.pr_url),
                blog_post_slug = COALESCE(EXCLUDED.blog_post_slug, public.pg_bio_research_attempts.blog_post_slug),
                notes = COALESCE(EXCLUDED.notes, public.pg_bio_research_attempts.notes),
                updated_at = now();
            """,
            (
                args.code,
                args.date,
                args.category,
                args.title,
                args.ref,
                args.strategy,
                args.worktree,
                args.branch,
                args.status,
                args.test_name,
                json.dumps(bench),
                args.pr_number,
                args.pr_url,
                args.blog_slug,
                args.notes,
            ),
        )
    conn.commit()
    export_json_snapshot(conn)
    print(f"Upserted research attempt {args.code} [{args.status}] and synced {LEDGER_JSON_PATH}")


def main():
    parser = argparse.ArgumentParser(description="Manage pg_bio daily research ledger and video backlog.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("init", help="Create SQL tables, seed historical attempts, and sync JSON ledger")
    sub.add_parser("status", help="Display current research attempts and video backlog")
    sub.add_parser("sync", help="Export SQL tables to data/research_ledger.json")

    up = sub.add_parser("upsert-attempt", help="Insert or update a research attempt")
    up.add_argument("--code", required=True)
    up.add_argument("--date", default=None)
    up.add_argument("--category", required=True)
    up.add_argument("--title", required=True)
    up.add_argument("--ref", required=True)
    up.add_argument("--strategy", required=True)
    up.add_argument("--worktree", default=None)
    up.add_argument("--branch", default=None)
    up.add_argument("--status", default="Proposed")
    up.add_argument("--test-name", default=None)
    up.add_argument("--benchmark-json", default="{}")
    up.add_argument("--pr-number", type=int, default=None)
    up.add_argument("--pr-url", default=None)
    up.add_argument("--blog-slug", default=None)
    up.add_argument("--notes", default=None)

    args = parser.parse_args()
    with psycopg.connect(DB_URI) as conn:
        init_and_seed(conn)
        if args.cmd in ("init", "status"):
            cmd_status(conn)
        elif args.cmd == "sync":
            export_json_snapshot(conn)
            print(f"Synced ledger to {LEDGER_JSON_PATH}")
        elif args.cmd == "upsert-attempt":
            cmd_upsert_attempt(conn, args)


if __name__ == "__main__":
    main()
