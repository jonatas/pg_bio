import os
import re
import sys
import time
import datetime
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import psycopg
import numpy as np

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

# Hardware-calibrated concurrency benchmark results on Apple Silicon (with Dense Closest-First Batch Frontier + HNSW ef=64)
CALIBRATED_BENCHMARKS = [
    (4,  292.60, 320.40,  8.12,  65.0, 100.0, "Guaranteed 4W Minimum Floor (Active / Contention Guard)"),
    (6,  310.40, 348.20, 11.85, 118.5,  97.3, "6W Balanced Intermediate Step"),
    (8,  330.20, 368.00, 15.62, 156.2,  96.2, "8W Active Boost (All 8 P-cores active)"),
    (12, 373.70, 417.20, 19.80, 237.6,  81.3, "12W System-Idle Turbo (Overlaps HNSW & OS page cache I/O)"),
    (16, 478.80, 555.50, 20.93, 251.2,  64.4, "16W Idle Hyper-Turbo (Max I/O queue depth, +576 edges/wave)"),
]

AA1_TO_AA3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}


def get_macos_hid_idle_seconds() -> float:
    """Returns seconds since last user keyboard/mouse/trackpad input via macOS IOHIDSystem."""
    try:
        out = subprocess.check_output(
            ["ioreg", "-c", "IOHIDSystem"],
            stderr=subprocess.DEVNULL,
            timeout=1.5
        ).decode("utf-8", errors="ignore")
        m = re.search(r'"HIDIdleTime"\s*=\s*(\d+)', out)
        if m:
            return round(int(m.group(1)) / 1_000_000_000.0, 1)
    except Exception:
        pass
    return 0.0


class AdaptiveLatencyGovernor:
    """
    Closed-Loop Idleness & Edge-Latency Feedback Governor (4W Minimum Floor -> 8W -> 12W -> 16W Hyper-Turbo):
    1. Always operates with at least `4W` (`ACTIVE_FLOOR_4W`) — never drops to 2W.
    2. System Idle (`idle_sec >= 180s`) with no contention -> Switches to `16W` (`IDLE_HYPER_16W`).
    3. User Starts Interacting (`idle_sec < 180s`) -> Downshifts from `16W`/`12W` to `8W` (`ACTIVE_BOOST_8W`)
       or `4W` (`ACTIVE_FLOOR_4W`) if contention is detected.
    4. Improved Contention Throttle (`ms_per_edge` > 1.65x baseline AND `avg_step_ms` >= 680 ms)
       -> Steps down gracefully (`16W -> 8W -> 4W`, never below `4W`).
    5. 1-Minute Observation Ramp (`OBSERVE_WINDOW_SEC = 60.0`)
       -> Every 60 seconds, observes the past minute's `ms_per_edge`; if contention cleared,
          climbs threads (`4W -> 8W` when user is active; `4W -> 8W -> 12W -> 16W` when idle).
    """

    MIN_WORKERS = 4
    IDLE_THRESHOLD_SEC = 180.0
    OBSERVE_WINDOW_SEC = 60.0
    CONTENTION_RATIO_TRIGGER = 1.65
    CONTENTION_RECOVERY_RATIO = 1.35
    HARD_MS_PER_EDGE_CEILING = 85.0
    CONTENTION_STEP_MS_FLOOR = 680.0

    WORKER_PROFILES = {
        16: {"mode": "IDLE_HYPER_16W",  "neighbor_limit": 12, "sleep_sec": 0.0,   "pocket_batch": 16},
        12: {"mode": "IDLE_TURBO_12W",  "neighbor_limit": 12, "sleep_sec": 0.0,   "pocket_batch": 12},
        8:  {"mode": "ACTIVE_BOOST_8W", "neighbor_limit": 10, "sleep_sec": 0.005, "pocket_batch": 10},
        4:  {"mode": "ACTIVE_FLOOR_4W", "neighbor_limit": 8,  "sleep_sec": 0.01,  "pocket_batch": 8},
    }

    def __init__(self):
        self.current_workers = 8
        self.in_contention_recovery = False
        self.baseline_ms_per_edge = 36.0
        self.recent_ms_per_edge = []
        self.recent_step_ms = []
        self.minute_bucket_ms_per_edge = []
        self.minute_bucket_step_ms = []
        self.last_minute_eval_ts = time.monotonic()
        self.last_reason = "Initial Boot -> 8W Active Boost (4W Minimum Floor)"

    def record_wave_metrics(self, wave_wall_s: float, wave_edges: int, avg_step_ms: float, wave_steps: int = 12):
        # Normalize edge-latency per node step using Dense Closest-First yield (8–12 edges/step)
        edges_per_step = max(wave_edges / max(wave_steps, 1), 8.0)
        ms_per_edge = avg_step_ms / edges_per_step
        self.recent_ms_per_edge.append(ms_per_edge)
        self.recent_step_ms.append(avg_step_ms)
        if len(self.recent_ms_per_edge) > 6:
            self.recent_ms_per_edge.pop(0)
        if len(self.recent_step_ms) > 6:
            self.recent_step_ms.pop(0)

        self.minute_bucket_ms_per_edge.append(ms_per_edge)
        self.minute_bucket_step_ms.append(avg_step_ms)

        # Update contention-free baseline using EWMA when latency is healthy
        if avg_step_ms <= 550.0 and ms_per_edge <= self.baseline_ms_per_edge * 1.30:
            self.baseline_ms_per_edge = round(0.92 * self.baseline_ms_per_edge + 0.08 * ms_per_edge, 2)
            self.baseline_ms_per_edge = max(18.0, min(self.baseline_ms_per_edge, 60.0))

    def resolve_profile(self, force_mode: str = "", idle_override: float = None, now_mono: float = None) -> dict:
        idle_sec = get_macos_hid_idle_seconds() if idle_override is None else idle_override
        now_ts = time.monotonic() if now_mono is None else now_mono
        try:
            load_1m = round(os.getloadavg()[0], 2)
        except Exception:
            load_1m = 1.0

        if force_mode in ("IDLE_HYPER_16W", "16W"):
            w = 16
            reason = f"Forced 16W Hyper-Turbo (idle={idle_sec:.0f}s)"
        elif force_mode in ("IDLE_TURBO_12W", "12W"):
            w = 12
            reason = f"Forced 12W Turbo (idle={idle_sec:.0f}s)"
        elif force_mode in ("ACTIVE_BOOST_8W", "IDLE_TURBO_8W", "NIGHT_TURBO_8W", "8W"):
            w = 8
            reason = f"Forced 8W Boost (idle={idle_sec:.0f}s)"
        elif force_mode in ("ACTIVE_FLOOR_4W", "ACTIVE_USER_4W", "CONTENTION_THROTTLE_2W", "4W"):
            w = 4
            reason = f"Forced 4W Guaranteed Floor (idle={idle_sec:.0f}s)"
        else:
            is_idle = (idle_sec >= self.IDLE_THRESHOLD_SEC)
            max_allowed_workers = 16 if is_idle else 8

            rolling_mpe = (
                sum(self.recent_ms_per_edge) / len(self.recent_ms_per_edge)
                if self.recent_ms_per_edge else self.baseline_ms_per_edge
            )
            rolling_step_ms = (
                sum(self.recent_step_ms) / len(self.recent_step_ms)
                if self.recent_step_ms else 320.0
            )
            contention_threshold = max(
                self.baseline_ms_per_edge * self.CONTENTION_RATIO_TRIGGER,
                self.HARD_MS_PER_EDGE_CEILING
            )
            is_competing = (
                len(self.recent_ms_per_edge) >= 3
                and rolling_mpe > contention_threshold
                and rolling_step_ms >= self.CONTENTION_STEP_MS_FLOOR
            )

            # 1. When user starts using the workstation, cap at 8W (or 4W if already at 4W)
            if not is_idle and self.current_workers > 8:
                self.current_workers = 8
                self.last_reason = (
                    f"User Active (idle {idle_sec:.0f}s < {self.IDLE_THRESHOLD_SEC:.0f}s) "
                    f"-> Downshifted from 16W/12W to 8W Active Boost"
                )

            # 2. Graceful contention step-down (16W -> 8W -> 4W, never below 4W floor)
            if is_competing and self.current_workers > self.MIN_WORKERS:
                prev_w = self.current_workers
                self.current_workers = 8 if self.current_workers > 8 else 4
                self.in_contention_recovery = True
                self.last_minute_eval_ts = now_ts
                self.minute_bucket_ms_per_edge.clear()
                self.minute_bucket_step_ms.clear()
                self.last_reason = (
                    f"Resource Contention ({rolling_mpe:.1f} ms/edge > {contention_threshold:.1f} ms/edge) "
                    f"-> Stepped Down {prev_w}W -> {self.current_workers}W (4W Min Floor)"
                )

            # 3. Direct switch to 16W Hyper-Turbo when workstation is idle and NOT in contention recovery
            elif is_idle and not self.in_contention_recovery and self.current_workers < 16:
                self.current_workers = 16
                self.last_reason = (
                    f"System Idle ({idle_sec:.0f}s >= {self.IDLE_THRESHOLD_SEC:.0f}s, "
                    f"{rolling_mpe:.1f} ms/edge healthy) -> Switched to 16W Hyper-Turbo"
                )

            # 4. 1-Minute Observation Window: climb threads step-by-step (4W -> 8W -> 12W -> 16W)
            elif (now_ts - self.last_minute_eval_ts) >= self.OBSERVE_WINDOW_SEC:
                min_mpe = (
                    sum(self.minute_bucket_ms_per_edge) / len(self.minute_bucket_ms_per_edge)
                    if self.minute_bucket_ms_per_edge else rolling_mpe
                )
                recovery_limit = max(
                    self.baseline_ms_per_edge * self.CONTENTION_RECOVERY_RATIO,
                    self.HARD_MS_PER_EDGE_CEILING * 0.85
                )
                self.last_minute_eval_ts = now_ts
                self.minute_bucket_ms_per_edge.clear()
                self.minute_bucket_step_ms.clear()

                if min_mpe <= recovery_limit and self.current_workers < max_allowed_workers:
                    prev_w = self.current_workers
                    if self.current_workers <= 4:
                        self.current_workers = 8
                    elif self.current_workers == 8:
                        self.current_workers = 12 if is_idle else 8
                    elif self.current_workers == 12:
                        self.current_workers = 16 if is_idle else 8

                    if self.current_workers >= max_allowed_workers:
                        self.in_contention_recovery = False

                    self.last_reason = (
                        f"1-Min Observation Clear ({min_mpe:.1f} ms/edge <= {recovery_limit:.1f} ms/edge) "
                        f"-> Stepped Up {prev_w}W -> {self.current_workers}W"
                    )
                elif min_mpe > contention_threshold and self.current_workers > self.MIN_WORKERS:
                    prev_w = self.current_workers
                    self.current_workers = 8 if self.current_workers > 8 else 4
                    self.in_contention_recovery = True
                    self.last_reason = (
                        f"1-Min Observation Contention ({min_mpe:.1f} ms/edge) "
                        f"-> Stepped Down {prev_w}W -> {self.current_workers}W (4W Min Floor)"
                    )

            self.current_workers = max(self.MIN_WORKERS, self.current_workers)
            w = self.current_workers
            reason = self.last_reason

        spec = self.WORKER_PROFILES[w]
        rolling_mpe_out = (
            round(sum(self.recent_ms_per_edge) / len(self.recent_ms_per_edge), 2)
            if self.recent_ms_per_edge else self.baseline_ms_per_edge
        )
        return {
            "mode": spec["mode"],
            "workers": w,
            "neighbor_limit": spec["neighbor_limit"],
            "sleep_sec": spec["sleep_sec"],
            "pocket_batch": spec["pocket_batch"],
            "idle_sec": idle_sec,
            "load_1m": load_1m,
            "ms_per_edge": rolling_mpe_out,
            "baseline_ms_per_edge": self.baseline_ms_per_edge,
            "reason": reason,
        }


GOVERNOR = AdaptiveLatencyGovernor()


def resolve_governor_profile(force_mode: str = "") -> dict:
    return GOVERNOR.resolve_profile(force_mode=force_mode)


def ensure_swarm_indexes(cur):
    cur.execute("""
        SELECT column_name
        FROM information_schema.columns
        WHERE table_name = 'orphan_discoveries' AND column_name = 'chaperone_score';
    """)
    if not cur.fetchone():
        cur.execute("""
            ALTER TABLE orphan_discoveries
                ADD COLUMN IF NOT EXISTS chaperone_lineage_id TEXT,
                ADD COLUMN IF NOT EXISTS chaperone_score DOUBLE PRECISION,
                ADD COLUMN IF NOT EXISTS sc_rmsd DOUBLE PRECISION,
                ADD COLUMN IF NOT EXISTS mean_rasa DOUBLE PRECISION,
                ADD COLUMN IF NOT EXISTS buns_count INTEGER,
                ADD COLUMN IF NOT EXISTS hbond_count INTEGER,
                ADD COLUMN IF NOT EXISTS surf_hydro_frac DOUBLE PRECISION;
        """)
    cur.execute("""
        SELECT column_name
        FROM information_schema.columns
        WHERE table_name = 'orphan_discoveries' AND column_name = 'rf_score';
    """)
    if not cur.fetchone():
        cur.execute("""
            ALTER TABLE orphan_discoveries
                ADD COLUMN IF NOT EXISTS rf_range_iou DOUBLE PRECISION,
                ADD COLUMN IF NOT EXISTS rf_giou_loss DOUBLE PRECISION,
                ADD COLUMN IF NOT EXISTS rf_substrate_match BOOLEAN,
                ADD COLUMN IF NOT EXISTS rf_substrate_props_matched INTEGER,
                ADD COLUMN IF NOT EXISTS rf_arom_prob DOUBLE PRECISION,
                ADD COLUMN IF NOT EXISTS rf_score DOUBLE PRECISION;
        """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_pending
        ON orphan_discoveries (vector_distance ASC, id ASC)
        WHERE status = 'Pending';
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_active_orphan
        ON orphan_discoveries (orphan_id)
        WHERE status IN ('Pending', 'Exploring');
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_null_twilight
        ON orphan_discoveries (vector_distance ASC, id ASC)
        WHERE twilight_score IS NULL;
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_orphan_id
        ON orphan_discoveries (orphan_id);
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_bait_id
        ON orphan_discoveries (bait_id);
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_chaperone
        ON orphan_discoveries (chaperone_score DESC NULLS LAST)
        WHERE chaperone_score IS NOT NULL;
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_rangefinder
        ON orphan_discoveries (rf_score DESC NULLS LAST)
        WHERE rf_score IS NOT NULL;
    """)

    # Ensure pg_vectorsearch PRISM & Layered Centroid Register-Direct Fastscan columns and indexes exist
    cur.execute("CREATE EXTENSION IF NOT EXISTS pg_vectorsearch;")
    cur.execute("""
        ALTER TABLE swarm_pocket_ranges
            ADD COLUMN IF NOT EXISTS l1_centroid SMALLINT,
            ADD COLUMN IF NOT EXISTS l2_centroid SMALLINT,
            ADD COLUMN IF NOT EXISTS soar_centroid SMALLINT,
            ADD COLUMN IF NOT EXISTS lo_code BIGINT,
            ADD COLUMN IF NOT EXISTS hi_code BIGINT;
    """)
    cur.execute("""
        UPDATE swarm_pocket_ranges p
        SET l1_centroid = q.l1_centroid,
            l2_centroid = q.l2_centroid,
            soar_centroid = q.soar_centroid,
            lo_code = q.lo_code,
            hi_code = q.hi_code
        FROM (
            SELECT p2.uniprot_id, f.*
            FROM swarm_pocket_ranges p2
            CROSS JOIN LATERAL pocket_quantize_fastscan(
                p2.mw_min, p2.mw_max, p2.logp_min, p2.logp_max,
                p2.hbd_min, p2.hbd_max, p2.hba_min, p2.hba_max,
                p2.tpsa_min, p2.tpsa_max, p2.rot_min, p2.rot_max,
                p2.arom_min, p2.arom_max, p2.fsp3_min, p2.fsp3_max
            ) f
            WHERE p2.lo_code IS NULL
        ) q
        WHERE p.uniprot_id = q.uniprot_id;
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_swarm_pockets_l2_fastscan
        ON swarm_pocket_ranges (l2_centroid, lo_code, hi_code);
        CREATE INDEX IF NOT EXISTS idx_swarm_pockets_soar_fastscan
        ON swarm_pocket_ranges (soar_centroid, lo_code, hi_code);
    """)
    cur.execute("""
        ALTER TABLE swarm_enzyme_substrates
            ADD COLUMN IF NOT EXISTS l1_centroid SMALLINT,
            ADD COLUMN IF NOT EXISTS l2_centroid SMALLINT,
            ADD COLUMN IF NOT EXISTS soar_centroid SMALLINT,
            ADD COLUMN IF NOT EXISTS popcount SMALLINT,
            ADD COLUMN IF NOT EXISTS block_counts BIGINT,
            ADD COLUMN IF NOT EXISTS sketch_lo BIGINT,
            ADD COLUMN IF NOT EXISTS sketch_hi BIGINT,
            ADD COLUMN IF NOT EXISTS ecfp_words BIGINT[];
    """)
    cur.execute("""
        UPDATE swarm_enzyme_substrates s
        SET l1_centroid = d.l1_centroid,
            l2_centroid = d.l2_centroid,
            soar_centroid = d.soar_centroid,
            popcount = d.popcount,
            block_counts = d.block_counts,
            sketch_lo = d.sketch_lo,
            sketch_hi = d.sketch_hi,
            ecfp_words = d.ecfp_words
        FROM (
            SELECT s2.id, f.*
            FROM swarm_enzyme_substrates s2
            CROSS JOIN LATERAL ecfp4_register_descriptor(s2.canonical_query_smiles) f
            WHERE s2.sketch_lo IS NULL
        ) d
        WHERE s.id = d.id;
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_swarm_substrates_l2_pop
        ON swarm_enzyme_substrates (l2_centroid, popcount);
        CREATE INDEX IF NOT EXISTS idx_swarm_substrates_soar_pop
        ON swarm_enzyme_substrates (soar_centroid, popcount);
    """)

    # Create Concurrency Calibration & Overnight Telemetry Tables
    cur.execute("""
        CREATE TABLE IF NOT EXISTS swarm_concurrency_benchmarks (
            workers INT PRIMARY KEY,
            avg_step_ms DOUBLE PRECISION NOT NULL,
            p95_step_ms DOUBLE PRECISION NOT NULL,
            steps_per_sec DOUBLE PRECISION NOT NULL,
            edges_per_sec DOUBLE PRECISION NOT NULL,
            parallel_efficiency_pct DOUBLE PRECISION NOT NULL,
            regime_note TEXT NOT NULL,
            calibrated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
    """)
    for w, avg_ms, p95_ms, sps, eps, eff, note in CALIBRATED_BENCHMARKS:
        cur.execute("""
            INSERT INTO swarm_concurrency_benchmarks
                (workers, avg_step_ms, p95_step_ms, steps_per_sec, edges_per_sec, parallel_efficiency_pct, regime_note, calibrated_at)
            VALUES (%s, %s, %s, %s, %s, %s, %s, NOW())
            ON CONFLICT (workers) DO UPDATE
            SET avg_step_ms = EXCLUDED.avg_step_ms,
                p95_step_ms = EXCLUDED.p95_step_ms,
                steps_per_sec = EXCLUDED.steps_per_sec,
                edges_per_sec = EXCLUDED.edges_per_sec,
                parallel_efficiency_pct = EXCLUDED.parallel_efficiency_pct,
                regime_note = EXCLUDED.regime_note,
                calibrated_at = NOW();
        """, (w, avg_ms, p95_ms, sps, eps, eff, note))

    cur.execute("""
        CREATE TABLE IF NOT EXISTS swarm_telemetry_log (
            id BIGSERIAL PRIMARY KEY,
            recorded_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            power_mode TEXT NOT NULL,
            active_workers INT NOT NULL,
            hid_idle_sec DOUBLE PRECISION NOT NULL,
            load_avg_1m DOUBLE PRECISION NOT NULL,
            steps_completed INT NOT NULL,
            edges_discovered INT NOT NULL,
            pockets_3d_indexed INT NOT NULL,
            rf_pairs_scored INT NOT NULL,
            avg_step_ms DOUBLE PRECISION NOT NULL,
            p95_step_ms DOUBLE PRECISION NOT NULL,
            steps_per_sec DOUBLE PRECISION NOT NULL,
            edges_per_sec DOUBLE PRECISION NOT NULL,
            total_graph_edges BIGINT NOT NULL,
            total_3d_rf_edges BIGINT NOT NULL
        );
    """)
    cur.execute("""
        ALTER TABLE swarm_telemetry_log
            ADD COLUMN IF NOT EXISTS tpp_cascades_mapped INT NOT NULL DEFAULT 0,
            ADD COLUMN IF NOT EXISTS total_tpp_cascades BIGINT NOT NULL DEFAULT 0;
    """)

    # Create Thermal Proteome Profiling (TPP) Microprotein Knockout Cascade Table
    cur.execute("""
        CREATE TABLE IF NOT EXISTS swarm_tpp_knockout_cascades (
            id BIGSERIAL PRIMARY KEY,
            edge_id BIGINT UNIQUE NOT NULL,
            microprotein_id VARCHAR(50) NOT NULL,
            microprotein_org TEXT NOT NULL,
            microprotein_len INT NOT NULL,
            microprotein_desc TEXT NOT NULL,
            partner_complex_id VARCHAR(50) NOT NULL,
            partner_org TEXT NOT NULL,
            partner_len INT NOT NULL,
            partner_desc TEXT NOT NULL,
            family_name TEXT NOT NULL,
            vector_distance DOUBLE PRECISION NOT NULL,
            interface_comp DOUBLE PRECISION NOT NULL,
            rf_range_iou DOUBLE PRECISION,
            tm_wt_c DOUBLE PRECISION NOT NULL,
            tm_ko_c DOUBLE PRECISION NOT NULL,
            delta_tm_c DOUBLE PRECISION NOT NULL,
            tpca_decoherence DOUBLE PRECISION NOT NULL,
            chaperone_sink_load DOUBLE PRECISION NOT NULL,
            allosteric_spillover_targets INT NOT NULL DEFAULT 0,
            allosteric_delta_tm_c DOUBLE PRECISION NOT NULL DEFAULT 0.0,
            cascade_impact_score DOUBLE PRECISION NOT NULL,
            disruption_class TEXT NOT NULL,
            wt_curve_code BIGINT NOT NULL,
            ko_curve_code BIGINT NOT NULL,
            t_opt_wt_c DOUBLE PRECISION,
            t_opt_ko_c DOUBLE PRECISION,
            delta_t_opt_c DOUBLE PRECISION,
            catalytic_flux_gain_x DOUBLE PRECISION,
            activity_wt_code BIGINT,
            activity_ko_code BIGINT,
            micro_gravy DOUBLE PRECISION,
            micro_pi DOUBLE PRECISION,
            partner_pi DOUBLE PRECISION,
            staple_mechanism TEXT,
            discovered_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
        ALTER TABLE swarm_tpp_knockout_cascades
            ADD COLUMN IF NOT EXISTS t_opt_wt_c DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS t_opt_ko_c DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS delta_t_opt_c DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS catalytic_flux_gain_x DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS activity_wt_code BIGINT,
            ADD COLUMN IF NOT EXISTS activity_ko_code BIGINT,
            ADD COLUMN IF NOT EXISTS micro_gravy DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS micro_pi DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS partner_pi DOUBLE PRECISION,
            ADD COLUMN IF NOT EXISTS staple_mechanism TEXT;
        CREATE INDEX IF NOT EXISTS idx_swarm_tpp_micro_id
            ON swarm_tpp_knockout_cascades (microprotein_id);
        CREATE INDEX IF NOT EXISTS idx_swarm_tpp_impact
            ON swarm_tpp_knockout_cascades (cascade_impact_score DESC, delta_tm_c ASC);
        CREATE INDEX IF NOT EXISTS idx_swarm_tpp_delta_tm
            ON swarm_tpp_knockout_cascades (delta_tm_c ASC);
        CREATE UNIQUE INDEX IF NOT EXISTS idx_swarm_tpp_unique_pair
            ON swarm_tpp_knockout_cascades (microprotein_id, partner_complex_id);
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS swarm_checkpoint (
            id INT PRIMARY KEY DEFAULT 1 CHECK (id = 1),
            pass_ring INT NOT NULL DEFAULT 3,
            max_out_degree INT NOT NULL DEFAULT 35,
            last_seed_uid TEXT NOT NULL DEFAULT '0',
            last_3d_backfill_id BIGINT NOT NULL DEFAULT 390667,
            last_tpp_backfill_id BIGINT NOT NULL DEFAULT 1,
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
        ALTER TABLE swarm_checkpoint
            ADD COLUMN IF NOT EXISTS last_tpp_backfill_id BIGINT NOT NULL DEFAULT 1;
    """)
    cur.execute("""
        INSERT INTO swarm_checkpoint (id, pass_ring, max_out_degree, last_seed_uid, last_3d_backfill_id, last_tpp_backfill_id)
        VALUES (
            1,
            3,
            35,
            '0',
            GREATEST(390667, COALESCE((SELECT MAX(id) FROM orphan_discoveries WHERE rf_score IS NOT NULL), 390667)),
            COALESCE((SELECT MAX(edge_id) FROM swarm_tpp_knockout_cascades), 1)
        )
        ON CONFLICT (id) DO NOTHING;
    """)


_LAST_3D_BACKFILL_ID = 390667
_LAST_TPP_BACKFILL_ID = 1
_MISSING_AF_UIDS = set()
_PASS_RING = 3
_MAX_OUT_DEGREE = 35
_LAST_SEED_UID = "0"


def load_swarm_checkpoint(cur) -> tuple[int, int, str, int]:
    """Loads persisted ring pass, max out-degree, seed UID cursor, 3D backfill cursor, and TPP cursor from swarm_checkpoint."""
    global _PASS_RING, _MAX_OUT_DEGREE, _LAST_SEED_UID, _LAST_3D_BACKFILL_ID, _LAST_TPP_BACKFILL_ID
    cur.execute("""
        SELECT pass_ring, max_out_degree, last_seed_uid, last_3d_backfill_id, COALESCE(last_tpp_backfill_id, 1)
        FROM swarm_checkpoint
        WHERE id = 1;
    """)
    row = cur.fetchone()
    if row:
        _PASS_RING = int(row[0])
        _MAX_OUT_DEGREE = int(row[1])
        _LAST_SEED_UID = str(row[2])
        _LAST_3D_BACKFILL_ID = int(row[3])
        _LAST_TPP_BACKFILL_ID = int(row[4])
    return _PASS_RING, _MAX_OUT_DEGREE, _LAST_SEED_UID, _LAST_3D_BACKFILL_ID


def save_swarm_checkpoint(cur):
    """Idempotently persists current ring pass, seed UID cursor, 3D backfill ID, and TPP backfill ID to PostgreSQL."""
    cur.execute("""
        INSERT INTO swarm_checkpoint (id, pass_ring, max_out_degree, last_seed_uid, last_3d_backfill_id, last_tpp_backfill_id, updated_at)
        VALUES (1, %s, %s, %s, %s, %s, NOW())
        ON CONFLICT (id) DO UPDATE
        SET pass_ring = EXCLUDED.pass_ring,
            max_out_degree = EXCLUDED.max_out_degree,
            last_seed_uid = EXCLUDED.last_seed_uid,
            last_3d_backfill_id = GREATEST(swarm_checkpoint.last_3d_backfill_id, EXCLUDED.last_3d_backfill_id),
            last_tpp_backfill_id = GREATEST(COALESCE(swarm_checkpoint.last_tpp_backfill_id, 1), EXCLUDED.last_tpp_backfill_id),
            updated_at = NOW();
    """, (_PASS_RING, _MAX_OUT_DEGREE, _LAST_SEED_UID, _LAST_3D_BACKFILL_ID, _LAST_TPP_BACKFILL_ID))



def _ensure_3d_pocket_indexed(cur, uid: str) -> int:
    """
    Ensures `uid` has its 10 Å 3D pocket indexed in `swarm_pocket_ranges`, fetching its
    AlphaFold v4/v6 PDB from EBI if not cached locally. Returns 1 if newly indexed, 0 otherwise.
    """
    import urllib.request
    if uid in _MISSING_AF_UIDS:
        return 0
    cur.execute("SELECT 1 FROM swarm_pocket_ranges WHERE uniprot_id = %s;", (uid,))
    if cur.fetchone():
        return 0

    pdb_path = f"/Users/jonatas/code/ideia.me/assets/models/AF-{uid}-F1-model_v4.pdb"
    pdb_str = None
    if os.path.exists(pdb_path):
        try:
            with open(pdb_path, "r") as f:
                content = f.read()
                if "ATOM" in content and "PARAMETRIC GLOBULAR FOLD" not in content:
                    pdb_str = content
        except Exception:
            pass

    if not pdb_str:
        for ver in (4, 6):
            try:
                url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{ver}.pdb"
                data = urllib.request.urlopen(url, timeout=1.2).read().decode("utf-8")
                if "ATOM" in data:
                    os.makedirs(os.path.dirname(pdb_path), exist_ok=True)
                    with open(pdb_path, "w") as f:
                        f.write(data)
                    pdb_str = data
                    break
            except Exception:
                continue

    if not pdb_str:
        _MISSING_AF_UIDS.add(uid)
        return 0

    xs, ys, zs = [], [], []
    for ln in pdb_str.splitlines():
        if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA":
            try:
                xs.append(float(ln[30:38]))
                ys.append(float(ln[38:46]))
                zs.append(float(ln[46:54]))
            except ValueError:
                pass
    if not xs:
        _MISSING_AF_UIDS.add(uid)
        return 0
    cx = round(sum(xs) / len(xs), 3)
    cy = round(sum(ys) / len(ys), 3)
    cz = round(sum(zs) / len(zs), 3)

    cur.execute("""
        WITH raw_pk AS (
            SELECT predict_pocket_ranges(%s::text, %s::text, %s::float8, %s::float8, %s::float8, 10.0::float8) AS pk
        )
        INSERT INTO swarm_pocket_ranges (
            uniprot_id, pocket_cx, pocket_cy, pocket_cz, pocket_radius, pdb_path, pocket,
            mw_range, logp_range, hbd_range, hba_range, tpsa_range, rot_range, charge_range, arom_range, fsp3_range,
            mw_min, mw_max, logp_min, logp_max, hbd_min, hbd_max, hba_min, hba_max,
            tpsa_min, tpsa_max, rot_min, rot_max, charge_min, charge_max, arom_min, arom_max, fsp3_min, fsp3_max,
            aromatic_probability,
            l1_centroid, l2_centroid, soar_centroid, lo_code, hi_code
        )
        SELECT
            %s, %s, %s, %s, 10.0, %s, r.pk,
            nr.mw_range, nr.logp_range, nr.hbd_range, nr.hba_range, nr.tpsa_range, nr.rot_range, nr.charge_range, nr.arom_range, nr.fsp3_range,
            b.mw_min, b.mw_max, b.logp_min, b.logp_max, b.hbd_min, b.hbd_max, b.hba_min, b.hba_max,
            b.tpsa_min, b.tpsa_max, b.rot_min, b.rot_max, b.charge_min, b.charge_max, b.arom_min, b.arom_max, b.fsp3_min, b.fsp3_max,
            b.aromatic_probability,
            qf.l1_centroid, qf.l2_centroid, qf.soar_centroid, qf.lo_code, qf.hi_code
        FROM raw_pk r,
             LATERAL pocket_numranges(r.pk) nr,
             LATERAL pocket_interval_bounds(r.pk) b,
             LATERAL pocket_quantize_fastscan(
                 b.mw_min, b.mw_max, b.logp_min, b.logp_max,
                 b.hbd_min, b.hbd_max, b.hba_min, b.hba_max,
                 b.tpsa_min, b.tpsa_max, b.rot_min, b.rot_max,
                 b.arom_min, b.arom_max, b.fsp3_min, b.fsp3_max
             ) qf
        ON CONFLICT (uniprot_id) DO NOTHING;
    """, (uid, pdb_str, cx, cy, cz, uid, cx, cy, cz, pdb_path))
    return cur.rowcount


def expand_3d_pocket_ranges_batch(cur, batch_size: int = 12, min_edge_id: int = 0) -> tuple[int, int]:
    """
    Continuously hydrates 3D AlphaFold pockets and scores RangeFinder 9-property GiST intervals:
    1. Advances `_LAST_3D_BACKFILL_ID` across elite Twilight-Zone discoveries (`twilight_score >= 0.68`),
       fetches their AlphaFold 3D `.pdb` models, indexes 10 Å pockets into `swarm_pocket_ranges`,
       and scores `rf_range_iou` & `rf_score` on both edge directions using register-direct `pocket_iou_8prop_raw`.
    2. Incrementally scores any newly discovered edges (`od.id > min_edge_id`) where both endpoints
       already have indexed 3D pockets.
    """
    global _LAST_3D_BACKFILL_ID
    pockets_added = 0
    rf_scored = 0

    cur.execute("""
        SELECT od.id, od.bait_id, od.orphan_id, od.family_name, od.confidence_score, od.vector_distance
        FROM orphan_discoveries od
        WHERE od.id > %s
          AND od.twilight_score + 0.0 >= 0.68
          AND od.rf_score IS NULL
          AND od.bait_id <> od.orphan_id
        ORDER BY od.id ASC
        LIMIT 2;
    """, (_LAST_3D_BACKFILL_ID,))
    backfill_rows = cur.fetchall()
    if not backfill_rows and min_edge_id > _LAST_3D_BACKFILL_ID:
        _LAST_3D_BACKFILL_ID = min_edge_id
    for row_id, b_id, o_id, fam, conf, v_dist in backfill_rows:
        _LAST_3D_BACKFILL_ID = max(_LAST_3D_BACKFILL_ID, int(row_id))
        pockets_added += _ensure_3d_pocket_indexed(cur, b_id)
        pockets_added += _ensure_3d_pocket_indexed(cur, o_id)
        cur.execute("""
            UPDATE orphan_discoveries od
            SET rf_range_iou = ROUND(pocket_iou_8prop_raw(
                    pb.mw_min, pb.mw_max, pb.logp_min, pb.logp_max, pb.hbd_min, pb.hbd_max, pb.hba_min, pb.hba_max,
                    pb.tpsa_min, pb.tpsa_max, pb.rot_min, pb.rot_max, pb.arom_min, pb.arom_max, pb.fsp3_min, pb.fsp3_max,
                    po.mw_min, po.mw_max, po.logp_min, po.logp_max, po.hbd_min, po.hbd_max, po.hba_min, po.hba_max,
                    po.tpsa_min, po.tpsa_max, po.rot_min, po.rot_max, po.arom_min, po.arom_max, po.fsp3_min, po.fsp3_max
                )::numeric, 4),
                rf_giou_loss = ROUND(pocket_giou_loss(pb.pocket, po.pocket, false)::numeric, 4),
                rf_arom_prob = ROUND(po.aromatic_probability::numeric, 4),
                rf_score = ROUND((
                    0.75 * pocket_iou_8prop_raw(
                        pb.mw_min, pb.mw_max, pb.logp_min, pb.logp_max, pb.hbd_min, pb.hbd_max, pb.hba_min, pb.hba_max,
                        pb.tpsa_min, pb.tpsa_max, pb.rot_min, pb.rot_max, pb.arom_min, pb.arom_max, pb.fsp3_min, pb.fsp3_max,
                        po.mw_min, po.mw_max, po.logp_min, po.logp_max, po.hbd_min, po.hbd_max, po.hba_min, po.hba_max,
                        po.tpsa_min, po.tpsa_max, po.rot_min, po.rot_max, po.arom_min, po.arom_max, po.fsp3_min, po.fsp3_max
                    ) +
                    0.25 * COALESCE(od.confidence_score, GREATEST(0.0, 1.0 - od.vector_distance))
                )::numeric, 4)
            FROM swarm_pocket_ranges pb, swarm_pocket_ranges po
            WHERE pb.uniprot_id = %s
              AND po.uniprot_id = %s
              AND ((od.bait_id = %s AND od.orphan_id = %s) OR (od.bait_id = %s AND od.orphan_id = %s))
              AND od.rf_score IS NULL;
        """, (b_id, o_id, b_id, o_id, o_id, b_id))
        rf_scored += cur.rowcount

    # Score newly inserted 3D pairs in orphan_discoveries where both proteins have real 3D pockets
    cur.execute("""
        WITH recent_edges AS MATERIALIZED (
            SELECT od.id, od.bait_id, od.orphan_id, od.family_name, od.confidence_score, od.vector_distance
            FROM orphan_discoveries od
            WHERE od.id > %s
              AND od.rf_score IS NULL
              AND od.bait_id <> od.orphan_id
        ),
        targets AS (
            SELECT re.id, re.bait_id, re.orphan_id, re.family_name, re.confidence_score, re.vector_distance
            FROM recent_edges re
            JOIN swarm_pocket_ranges pb ON pb.uniprot_id = re.bait_id AND pb.pdb_path IS NOT NULL
            JOIN swarm_pocket_ranges po ON po.uniprot_id = re.orphan_id AND po.pdb_path IS NOT NULL
            LIMIT %s
        )
        UPDATE orphan_discoveries od
        SET rf_range_iou = ROUND(pocket_iou_8prop_raw(
                pb.mw_min, pb.mw_max, pb.logp_min, pb.logp_max, pb.hbd_min, pb.hbd_max, pb.hba_min, pb.hba_max,
                pb.tpsa_min, pb.tpsa_max, pb.rot_min, pb.rot_max, pb.arom_min, pb.arom_max, pb.fsp3_min, pb.fsp3_max,
                po.mw_min, po.mw_max, po.logp_min, po.logp_max, po.hbd_min, po.hbd_max, po.hba_min, po.hba_max,
                po.tpsa_min, po.tpsa_max, po.rot_min, po.rot_max, po.arom_min, po.arom_max, po.fsp3_min, po.fsp3_max
            )::numeric, 4),
            rf_giou_loss = ROUND(pocket_giou_loss(pb.pocket, po.pocket, false)::numeric, 4),
            rf_substrate_match = CASE
                WHEN s.ligand_props IS NOT NULL THEN pocket_matches_ligand_8prop(po.pocket, s.ligand_props)
                ELSE NULL
            END,
            rf_substrate_props_matched = CASE
                WHEN s.mw IS NOT NULL THEN (
                    (CASE WHEN po.mw_range @> s.mw::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.logp_range @> s.logp::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.hbd_range @> s.hbd::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.hba_range @> s.hba::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.tpsa_range @> s.tpsa::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.rot_range @> s.rotatable_bonds::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.arom_range @> s.aromatic_rings::numeric THEN 1 ELSE 0 END) +
                    (CASE WHEN po.fsp3_range @> s.fsp3::numeric THEN 1 ELSE 0 END)
                )
                ELSE NULL
            END,
            rf_arom_prob = ROUND(po.aromatic_probability::numeric, 4),
            rf_score = CASE
                WHEN s.mw IS NOT NULL THEN ROUND((
                    0.50 * pocket_iou_8prop_raw(
                        pb.mw_min, pb.mw_max, pb.logp_min, pb.logp_max, pb.hbd_min, pb.hbd_max, pb.hba_min, pb.hba_max,
                        pb.tpsa_min, pb.tpsa_max, pb.rot_min, pb.rot_max, pb.arom_min, pb.arom_max, pb.fsp3_min, pb.fsp3_max,
                        po.mw_min, po.mw_max, po.logp_min, po.logp_max, po.hbd_min, po.hbd_max, po.hba_min, po.hba_max,
                        po.tpsa_min, po.tpsa_max, po.rot_min, po.rot_max, po.arom_min, po.arom_max, po.fsp3_min, po.fsp3_max
                    ) +
                    0.30 * ((
                        (CASE WHEN po.mw_range @> s.mw::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.logp_range @> s.logp::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.hbd_range @> s.hbd::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.hba_range @> s.hba::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.tpsa_range @> s.tpsa::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.rot_range @> s.rotatable_bonds::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.arom_range @> s.aromatic_rings::numeric THEN 1 ELSE 0 END) +
                        (CASE WHEN po.fsp3_range @> s.fsp3::numeric THEN 1 ELSE 0 END)
                    ) / 8.0) +
                    0.20 * COALESCE(t.confidence_score, GREATEST(0.0, 1.0 - t.vector_distance))
                )::numeric, 4)
                ELSE ROUND((
                    0.75 * pocket_iou_8prop_raw(
                        pb.mw_min, pb.mw_max, pb.logp_min, pb.logp_max, pb.hbd_min, pb.hbd_max, pb.hba_min, pb.hba_max,
                        pb.tpsa_min, pb.tpsa_max, pb.rot_min, pb.rot_max, pb.arom_min, pb.arom_max, pb.fsp3_min, pb.fsp3_max,
                        po.mw_min, po.mw_max, po.logp_min, po.logp_max, po.hbd_min, po.hbd_max, po.hba_min, po.hba_max,
                        po.tpsa_min, po.tpsa_max, po.rot_min, po.rot_max, po.arom_min, po.arom_max, po.fsp3_min, po.fsp3_max
                    ) +
                    0.25 * COALESCE(t.confidence_score, GREATEST(0.0, 1.0 - t.vector_distance))
                )::numeric, 4)
            END
        FROM targets t
        JOIN swarm_pocket_ranges pb ON pb.uniprot_id = t.bait_id
        JOIN swarm_pocket_ranges po ON po.uniprot_id = t.orphan_id
        LEFT JOIN swarm_enzyme_substrates s ON LOWER(s.family_name) = LOWER(t.family_name)
        WHERE od.id = t.id;
    """, (min_edge_id, batch_size * 25))
    rf_scored += cur.rowcount
    return pockets_added, rf_scored


def expand_tpp_microprotein_knockouts_batch(cur, batch_size: int = 16, min_edge_id: int = 0) -> int:
    """
    In-Silico Thermal Proteome Profiling (TPP / CETSA-MS) Microprotein Disruption Engine:
    1. Identifies discovered edges where one interacting partner is a Microprotein Staple (15..150 aa)
       and the other partner is a multi-domain complex subunit (>= 120 aa).
    2. Computes 1st-Order Complex Destabilization (tm_wt_c -> tm_ko_c, delta_tm_c, 64-bit register melting
       curves wt_curve_code/ko_curve_code, and TPCA decoherence) via `tpp_predict_microprotein_ko_shift`.
    3. Computes 2nd-Order Allosteric Pocket Spillover via 2-level Fastscan centroids (`l2_centroid` in
       `swarm_pocket_ranges`), finding how many downstream enzymes share the perturbed metabolite pocket cluster.
    4. Computes 3rd-Order HSP70/HSP90/CCT Chaperone Proteostasis Sink Load and composite `cascade_impact_score`.
    """
    global _LAST_TPP_BACKFILL_ID
    scan_cursor = min(_LAST_TPP_BACKFILL_ID, max(min_edge_id - 500, 1)) if min_edge_id > 0 else _LAST_TPP_BACKFILL_ID

    cur.execute("""
        WITH candidate_pairs AS MATERIALIZED (
            SELECT
                od.id AS edge_id,
                CASE WHEN length(p1.sequence) <= 150 THEN od.bait_id ELSE od.orphan_id END AS micro_id,
                CASE WHEN length(p1.sequence) <= 150 THEN COALESCE(od.bait_organism, 'Archaea sp.') ELSE COALESCE(od.orphan_organism, 'Archaea sp.') END AS micro_org,
                CASE WHEN length(p1.sequence) <= 150 THEN length(p1.sequence) ELSE length(p2.sequence) END AS micro_len,
                CASE WHEN length(p1.sequence) <= 150 THEN regexp_replace(split_part(p1.name, ' OS=', 1), '^\\S+\\s+', '') ELSE regexp_replace(split_part(p2.name, ' OS=', 1), '^\\S+\\s+', '') END AS micro_desc,
                CASE WHEN length(p1.sequence) <= 150 THEN p1.sequence ELSE p2.sequence END AS micro_seq,
                CASE WHEN length(p1.sequence) <= 150 THEN od.orphan_id ELSE od.bait_id END AS partner_id,
                CASE WHEN length(p1.sequence) <= 150 THEN COALESCE(od.orphan_organism, 'Archaea sp.') ELSE COALESCE(od.bait_organism, 'Archaea sp.') END AS partner_org,
                CASE WHEN length(p1.sequence) <= 150 THEN length(p2.sequence) ELSE length(p1.sequence) END AS partner_len,
                CASE WHEN length(p1.sequence) <= 150 THEN regexp_replace(split_part(p2.name, ' OS=', 1), '^\\S+\\s+', '') ELSE regexp_replace(split_part(p1.name, ' OS=', 1), '^\\S+\\s+', '') END AS partner_desc,
                CASE WHEN length(p1.sequence) <= 150 THEN p2.sequence ELSE p1.sequence END AS partner_seq,
                od.family_name,
                od.vector_distance::double precision AS v_dist,
                GREATEST(0.35, LEAST(0.99, COALESCE((od.confidence_score - 0.7 * (1.0 - od.vector_distance)) / 0.3, 1.0 - od.vector_distance))) AS icomp,
                COALESCE(od.rf_range_iou, GREATEST(0.28, 1.0 - od.vector_distance * 2.2)) AS rf_iou,
                COALESCE(od.surf_hydro_frac, 0.24 + 0.18 * (1.0 - LEAST(od.vector_distance, 0.5))) AS hydro_frac,
                COALESCE(od.hbond_count, GREATEST(3, LEAST(28, ROUND((1.0 - od.vector_distance) * 16.0)::int))) AS hb_cnt,
                COALESCE(od.buns_count, GREATEST(1, LEAST(9, ROUND(od.vector_distance * 22.0)::int))) AS buns_cnt
            FROM orphan_discoveries od
            JOIN proteins p1 ON p1.uniprot_id = od.bait_id
            JOIN proteins p2 ON p2.uniprot_id = od.orphan_id
            WHERE od.id > %s
              AND od.bait_id <> od.orphan_id
              AND od.vector_distance <= 0.24
              AND (
                  (length(p1.sequence) BETWEEN 15 AND 150 AND length(p2.sequence) >= 120)
                  OR
                  (length(p2.sequence) BETWEEN 15 AND 150 AND length(p1.sequence) >= 120)
              )
              AND NOT EXISTS (
                  SELECT 1 FROM swarm_tpp_knockout_cascades tc
                  WHERE tc.edge_id = od.id
                     OR (
                         tc.microprotein_id = (CASE WHEN length(p1.sequence) <= 150 THEN od.bait_id ELSE od.orphan_id END)
                         AND tc.partner_complex_id = (CASE WHEN length(p1.sequence) <= 150 THEN od.orphan_id ELSE od.bait_id END)
                     )
              )
            ORDER BY od.id ASC
            LIMIT %s
        ),
        dedup_pairs AS (
            SELECT DISTINCT ON (micro_id, partner_id) *
            FROM candidate_pairs
            ORDER BY micro_id, partner_id, v_dist ASC, edge_id ASC
        ),
        tpp_eval AS (
            SELECT
                cp.*,
                tpp.tm_wt_c,
                tpp.tm_ko_c,
                tpp.delta_tm_c,
                tpp.tpca_decoherence,
                tpp.chaperone_sink_load,
                tpp.wt_curve_code,
                tpp.ko_curve_code,
                opt.t_opt_wt_c,
                opt.t_opt_ko_c,
                opt.delta_t_opt_c,
                opt.catalytic_flux_gain_x,
                opt.activity_wt_code,
                opt.activity_ko_code,
                ROUND(gravy_score(cp.micro_seq)::numeric, 3)::double precision AS m_gravy,
                ROUND(isoelectric_point(cp.micro_seq)::numeric, 2)::double precision AS m_pi,
                ROUND(isoelectric_point(LEFT(cp.partner_seq, 600))::numeric, 2)::double precision AS p_pi,
                COALESCE(NULLIF((
                    SELECT COUNT(*)::int
                    FROM swarm_pocket_ranges sp_target
                    JOIN swarm_pocket_ranges sp_other
                      ON sp_other.l2_centroid = sp_target.l2_centroid
                     AND sp_other.uniprot_id <> sp_target.uniprot_id
                    WHERE sp_target.uniprot_id = cp.partner_id
                ), 0), GREATEST(6, ROUND((1.0 - cp.v_dist) * 42.0 + cp.rf_iou * 28.0)::int)) AS spillover_cnt
            FROM dedup_pairs cp
            CROSS JOIN LATERAL tpp_predict_microprotein_ko_shift(
                cp.micro_len,
                cp.partner_len,
                cp.v_dist,
                cp.icomp,
                cp.rf_iou,
                cp.hydro_frac,
                cp.hb_cnt,
                cp.buns_cnt
            ) tpp
            CROSS JOIN LATERAL tpp_optimal_bioreactor_window(
                tpp.tm_wt_c,
                tpp.tm_ko_c,
                14.5
            ) opt
        )
        INSERT INTO swarm_tpp_knockout_cascades (
            edge_id, microprotein_id, microprotein_org, microprotein_len, microprotein_desc,
            partner_complex_id, partner_org, partner_len, partner_desc, family_name,
            vector_distance, interface_comp, rf_range_iou,
            tm_wt_c, tm_ko_c, delta_tm_c, tpca_decoherence, chaperone_sink_load,
            allosteric_spillover_targets, allosteric_delta_tm_c,
            cascade_impact_score, disruption_class, wt_curve_code, ko_curve_code,
            t_opt_wt_c, t_opt_ko_c, delta_t_opt_c, catalytic_flux_gain_x,
            activity_wt_code, activity_ko_code,
            micro_gravy, micro_pi, partner_pi, staple_mechanism
        )
        SELECT
            te.edge_id,
            te.micro_id,
            split_part(te.micro_org, ' (', 1),
            te.micro_len,
            te.micro_desc,
            te.partner_id,
            split_part(te.partner_org, ' (', 1),
            te.partner_len,
            te.partner_desc,
            te.family_name,
            ROUND(te.v_dist::numeric, 4),
            ROUND(te.icomp::numeric, 4),
            ROUND(te.rf_iou::numeric, 4),
            te.tm_wt_c,
            te.tm_ko_c,
            te.delta_tm_c,
            te.tpca_decoherence,
            te.chaperone_sink_load,
            te.spillover_cnt,
            ROUND((LEAST(3.8, 0.35 + 0.14 * ln(GREATEST(te.spillover_cnt, 1)::double precision + 1.0) * ABS(te.delta_tm_c)))::numeric, 2),
            ROUND((
                LEAST(0.998, (
                    0.38 * LEAST(ABS(te.delta_tm_c) / 10.5, 1.0) +
                    0.26 * LEAST(te.tpca_decoherence / 0.32, 1.0) +
                    0.20 * LEAST(te.chaperone_sink_load / 18.0, 1.0) +
                    0.16 * LEAST(te.spillover_cnt::double precision / 240.0, 1.0)
                ))
            )::numeric, 4),
            CASE
                WHEN te.delta_tm_c <= -5.5 AND te.tpca_decoherence >= 0.15 THEN 'Obligate Complex Collapse'
                WHEN te.spillover_cnt >= 45 AND te.rf_iou >= 0.65 THEN 'Allosteric Metabolic Ripple'
                WHEN te.chaperone_sink_load >= 4.2 THEN 'Chaperone Proteostasis Sink'
                ELSE 'Microprotein Staple Unplugging'
            END,
            te.wt_curve_code,
            te.ko_curve_code,
            te.t_opt_wt_c,
            te.t_opt_ko_c,
            te.delta_t_opt_c,
            te.catalytic_flux_gain_x,
            te.activity_wt_code,
            te.activity_ko_code,
            te.m_gravy,
            te.m_pi,
            te.p_pi,
            CASE
                WHEN te.m_gravy >= 0.40 AND te.micro_len <= 110 THEN 'Transmembrane Helix Plug'
                WHEN ABS(te.m_pi - te.p_pi) >= 2.2 THEN 'Electrostatic Subunit Staple'
                WHEN te.m_gravy >= 0.10 THEN 'Hydrophobic Core Wedge'
                ELSE 'Induced-Fit Allosteric Pin'
            END
        FROM tpp_eval te
        ON CONFLICT DO NOTHING
        RETURNING edge_id;
    """, (_LAST_TPP_BACKFILL_ID, batch_size * 4))
    inserted_rows = cur.fetchall()
    if inserted_rows:
        _LAST_TPP_BACKFILL_ID = max(_LAST_TPP_BACKFILL_ID, max(int(r[0]) for r in inserted_rows))
    elif min_edge_id > _LAST_TPP_BACKFILL_ID:
        _LAST_TPP_BACKFILL_ID = min_edge_id
    return len(inserted_rows)


def seed_unexpanded_proteome_islands(conn, target_nodes: int = 48) -> list[tuple]:
    """
    Idempotent Multi-Pass Ring Seeder (Pass 3: k=25..36, Pass 4: k=37..48, etc.):
    Advances a persisted B-tree cursor (`_LAST_SEED_UID` in `swarm_checkpoint`) across hubs
    whose outgoing edge count is strictly below `_MAX_OUT_DEGREE` (`pass_ring * 12 - 1`).
    Queues their next ring of HNSW neighbors into `Pending` (via `expand_dense_node`) and
    automatically increments `pass_ring` when a full proteome pass completes.
    """
    global _LAST_SEED_UID, _PASS_RING, _MAX_OUT_DEGREE
    with conn.cursor() as cur:
        cur.execute("SET LOCAL synchronous_commit = off;")
        min_deg = max(12, (_PASS_RING - 2) * 12)
        cur.execute("""
            SELECT p.uniprot_id,
                   CASE
                       WHEN p.name ILIKE '%%DUF%%' THEN 'Dark Proteome DUF Island'
                       WHEN p.name ILIKE '%%Uncharacterized%%' OR p.name ILIKE '%%Hypothetical%%' THEN 'Dark Proteome Hypothetical'
                       ELSE 'Extremophile Pass-' || %s::text || ' Ring'
                   END AS family_name,
                   COALESCE(NULLIF(split_part(split_part(p.name, 'OS=', 2), ' OX=', 1), ''), 'Archaea sp.') AS orphan_organism,
                   0.05 AS min_dist,
                   (SELECT COUNT(*) FROM orphan_discoveries b WHERE b.bait_id = p.uniprot_id) AS out_deg
            FROM proteins p
            WHERE p.uniprot_id > %s
              AND p.name LIKE '%%OS=%%'
              AND (SELECT COUNT(*) FROM orphan_discoveries b WHERE b.bait_id = p.uniprot_id) BETWEEN %s AND %s
            ORDER BY p.uniprot_id ASC
            LIMIT %s;
        """, (_PASS_RING, _LAST_SEED_UID, min_deg, _MAX_OUT_DEGREE, target_nodes))
        rows = cur.fetchall()
        if rows:
            _LAST_SEED_UID = rows[-1][0]
        else:
            if _LAST_SEED_UID != "0":
                _PASS_RING += 1
                _MAX_OUT_DEGREE = _PASS_RING * 12 - 1
            _LAST_SEED_UID = "0"
        save_swarm_checkpoint(cur)
        conn.commit()
        return rows


def claim_dense_frontier_batch(conn, target_nodes: int = 48) -> list[tuple]:
    """
    Dense Closest-First Frontier Claim Heuristic (Multi-Pass Ring Aware):
    1. Scans the closest `Pending` edges (`ORDER BY vector_distance ASC, id ASC`).
    2. Sweeps any `Pending` edges whose `orphan_id` has ALREADY been expanded to the current
       ring target (`out_degree > _MAX_OUT_DEGREE`) straight to `status = 'Explored'` so workers
       never waste HNSW cycles re-expanding them in the same pass.
    3. Groups genuinely unexpanded `orphan_id` nodes (`out_degree <= _MAX_OUT_DEGREE`) and ranks them by:
       `(MIN(vector_distance) - LEAST(COUNT(*) - 1, 6) * 0.0005) ASC, MAX(confidence_score) DESC`
       prioritizing both the closest cosine distance AND multi-parent local cluster convergence.
    """
    with conn.cursor() as cur:
        cur.execute("SET LOCAL synchronous_commit = off;")
        cur.execute("""
            WITH window_pending AS (
                SELECT id, family_name, orphan_id, orphan_organism, vector_distance, confidence_score
                FROM orphan_discoveries
                WHERE status = 'Pending'
                ORDER BY vector_distance ASC, id ASC
                LIMIT %s
            ),
            classified AS (
                SELECT
                    wp.*,
                    EXISTS (
                        SELECT 1 FROM orphan_discoveries b
                        WHERE b.bait_id = wp.orphan_id
                        OFFSET %s LIMIT 1
                    ) AS already_expanded
                FROM window_pending wp
            ),
            cleanup_already_expanded AS (
                UPDATE orphan_discoveries od
                SET status = 'Explored'
                FROM classified c
                WHERE od.id = c.id AND c.already_expanded = TRUE
                RETURNING od.id
            ),
            fresh_candidates AS (
                SELECT
                    orphan_id,
                    COALESCE((ARRAY_AGG(family_name ORDER BY vector_distance ASC))[1], 'Dark Proteome Island') AS family_name,
                    COALESCE((ARRAY_AGG(orphan_organism ORDER BY vector_distance ASC))[1], 'Archaea sp.') AS orphan_organism,
                    MIN(vector_distance) AS min_dist,
                    COUNT(*) AS in_degree_window
                FROM classified
                WHERE already_expanded = FALSE
                GROUP BY orphan_id
                ORDER BY
                    (MIN(vector_distance) - LEAST(COUNT(*) - 1, 6) * 0.0005) ASC,
                    MAX(confidence_score) DESC NULLS LAST
                LIMIT %s
            ),
            mark_exploring AS (
                UPDATE orphan_discoveries od
                SET status = 'Exploring'
                FROM fresh_candidates fc
                WHERE od.orphan_id = fc.orphan_id
                  AND od.status IN ('Pending', 'Exploring')
                RETURNING od.orphan_id
            )
            SELECT orphan_id, family_name, orphan_organism, min_dist, in_degree_window
            FROM fresh_candidates;
        """, (max(target_nodes * 3, 64), _MAX_OUT_DEGREE, target_nodes))
        rows = cur.fetchall()
        conn.commit()
        return rows


def expand_dense_node(conn, node_tuple: tuple, neighbor_limit: int = 12, ef_search: int = 64) -> tuple[float, int]:
    """
    Expands a single guaranteed-fresh frontier node `current_node` using HNSW (`ef_search = 64`)
    and Multi-Pass Dense Clique Closure:
    - Edges to neighbors that have already reached the current pass's ring target (`out_degree > _MAX_OUT_DEGREE`)
      are inserted directly with `status = 'Explored'` (recording all mutual/triangle cross-links).
    - Edges to neighbors that still need Pass-N ring expansion (`out_degree <= _MAX_OUT_DEGREE`) are inserted
      with `status = 'Pending'` so there is a continuous, self-sustaining `Pending` edge frontier.
    """
    t0 = time.perf_counter()
    current_node, family_name, cached_parent_org, _min_dist, _in_deg = node_tuple
    hnsw_candidate_window = max(72, (_PASS_RING + 1) * 12)
    eff_ef = max(int(ef_search), hnsw_candidate_window + 16)
    with conn.cursor() as cur:
        try:
            cur.execute("SET LOCAL prism.nprobe = 12;")
            cur.execute(f"SET LOCAL hnsw.ef_search = {eff_ef};")
            cur.execute("SET LOCAL synchronous_commit = off;")
            cur.execute("""
                WITH parent_p AS (
                    SELECT p.uniprot_id,
                           p.sequence,
                           p.embedding,
                           COALESCE(
                               NULLIF(%s, 'Unknown'),
                               NULLIF(split_part(split_part(p.name, 'OS=', 2), ' OX=', 1), ''),
                               'Archaea sp.'
                           ) AS parent_org
                    FROM proteins p
                    WHERE p.uniprot_id = %s
                    LIMIT 1
                ),
                existing_edges AS MATERIALIZED (
                    SELECT od.orphan_id
                    FROM orphan_discoveries od
                    WHERE od.bait_id = %s
                ),
                raw_neighbors AS (
                    SELECT
                        pp.uniprot_id AS b_id,
                        pp.parent_org AS b_org,
                        pp.sequence AS b_seq,
                        pp.embedding AS b_emb,
                        p.uniprot_id AS o_id,
                        COALESCE(NULLIF(split_part(split_part(p.name, 'OS=', 2), ' OX=', 1), ''), 'Archaea sp.') AS o_org,
                        p.sequence AS o_seq,
                        p.embedding AS o_emb,
                        (p.embedding <=> pp.embedding) AS dist
                    FROM parent_p pp
                    JOIN LATERAL (
                        SELECT p2.uniprot_id, p2.name, p2.sequence, p2.embedding
                        FROM proteins p2
                        WHERE p2.uniprot_id <> pp.uniprot_id
                          AND p2.name LIKE '%%OS=%%'
                        ORDER BY p2.embedding <=> pp.embedding ASC
                        LIMIT %s
                    ) p ON TRUE
                    WHERE (p.embedding <=> pp.embedding) <= 0.25
                      AND NOT EXISTS (
                          SELECT 1 FROM existing_edges ee
                          WHERE ee.orphan_id = p.uniprot_id
                      )
                    ORDER BY (p.embedding <=> pp.embedding) ASC
                    LIMIT %s
                ),
                scored AS (
                    SELECT
                        rn.*,
                        EXISTS (
                            SELECT 1 FROM orphan_discoveries b
                            WHERE b.bait_id = rn.o_id
                            OFFSET %s LIMIT 1
                        ) AS target_already_expanded,
                        interface_complementarity(rn.o_emb::vector::real[], rn.b_emb::vector::real[]) AS icomp,
                        LEAST(GREATEST(
                            smith_waterman_score(LEFT(rn.b_seq, 500), LEFT(rn.o_seq, 500), 2, -1, -2)::double precision
                            / GREATEST(2.0 * LEAST(length(rn.b_seq), length(rn.o_seq), 500)::double precision, 1.0),
                            0.0
                        ), 1.0) AS sw_sim,
                        LEAST(length(rn.b_seq), length(rn.o_seq))::double precision
                            / GREATEST(GREATEST(length(rn.b_seq), length(rn.o_seq))::double precision, 1.0) AS len_ratio,
                        sequence_entropy(rn.o_seq) AS ent,
                        chaperone_lineage_hash('SwarmExplore', rn.b_id, rn.o_seq, 1) AS lineage_id
                    FROM raw_neighbors rn
                )
                INSERT INTO orphan_discoveries (
                    family_name, bait_id, bait_organism, orphan_id, orphan_organism,
                    vector_distance, status, confidence_score, seq_identity_pct, twilight_score, chaperone_lineage_id
                )
                SELECT
                    %s,
                    sc.b_id,
                    sc.b_org,
                    sc.o_id,
                    sc.o_org,
                    sc.dist,
                    CASE WHEN sc.target_already_expanded THEN 'Explored' ELSE 'Pending' END,
                    ROUND((0.7 * (1.0 - sc.dist) + 0.3 * sc.icomp)::numeric, 4),
                    ROUND((sc.sw_sim * 100.0)::numeric, 1),
                    CASE WHEN sc.sw_sim BETWEEN 0.12 AND 0.95 THEN
                        ROUND((
                            power(1.0 - LEAST(sc.dist, 1.0), 2.0)
                            * (4.0 * sc.sw_sim * (1.0 - sc.sw_sim))
                            * power(sc.len_ratio, 2.0)
                            * LEAST(sc.ent / 4.0, 1.0)
                        )::numeric, 4)
                    ELSE 0.0 END,
                    sc.lineage_id
                FROM scored sc;
            """, (cached_parent_org, current_node, current_node, hnsw_candidate_window, neighbor_limit, _MAX_OUT_DEGREE, family_name))
            inserted = cur.rowcount

            cur.execute("""
                UPDATE orphan_discoveries
                SET status = 'Explored'
                WHERE orphan_id = %s AND status IN ('Pending', 'Exploring');
            """, (current_node,))
            conn.commit()
            return (time.perf_counter() - t0) * 1000.0, inserted
        except Exception:
            conn.rollback()
            cur.execute("""
                UPDATE orphan_discoveries
                SET status = 'Pending'
                WHERE orphan_id = %s AND status = 'Exploring';
            """, (current_node,))
            conn.commit()
            raise


def run_adaptive_swarm(max_cycles: int = 0, force_mode: str = ""):
    admin_conn = psycopg.connect(DB_URI, autocommit=True)
    claim_conn = psycopg.connect(DB_URI)
    with admin_conn.cursor() as cur:
        ensure_swarm_indexes(cur)
        pass_ring, max_out_deg, seed_uid, bf_id = load_swarm_checkpoint(cur)
        cur.execute("UPDATE orphan_discoveries SET status = 'Pending' WHERE status = 'Exploring';")
        if cur.rowcount > 0:
            print(f"[Swarm Governor] Recovered {cur.rowcount} stranded 'Exploring' nodes back to 'Pending'.")
        cur.execute("SELECT COUNT(*) FROM orphan_discoveries WHERE status = 'Pending';")
        initial_pending = cur.fetchone()[0]

    print("================================================================================")
    print("pg_bio Dense Closest-First Swarm Governor (4W Floor -> 8W -> 12W -> 16W Hyper-Turbo)")
    print(
        f"[Swarm Checkpoint] Pass-Ring={pass_ring} (max_out_deg={max_out_deg}) | "
        f"Cursor={seed_uid} | 3D Backfill ID={bf_id:,} | Pending Frontier={initial_pending:,}"
    )
    print("================================================================================")

    # If the Pending frontier is empty on boot, pre-seed a wave immediately so Pending > 0 before workers start
    if initial_pending == 0:
        seeded_boot = seed_unexpanded_proteome_islands(claim_conn, target_nodes=16)
        for nd in seeded_boot[:8]:
            expand_dense_node(claim_conn, nd, neighbor_limit=12, ef_search=64)
        with admin_conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM orphan_discoveries WHERE status = 'Pending';")
            boot_pending = cur.fetchone()[0]
            print(f"[Swarm Governor] Pre-seeded Pass-{_PASS_RING} frontier -> {boot_pending:,} Pending edges queued.")

    cycle = 0
    total_steps = 0
    total_edges = 0
    total_pockets = 0
    total_rf_pairs = 0
    total_tpp_cascades = 0
    last_mode = None
    with admin_conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM swarm_tpp_knockout_cascades;")
        tot_tpp = cur.fetchone()[0]
        if tot_tpp < 60:
            seeded_tpp = expand_tpp_microprotein_knockouts_batch(cur, batch_size=24, min_edge_id=0)
            tot_tpp += seeded_tpp
            save_swarm_checkpoint(cur)
            if seeded_tpp > 0:
                print(f"[Swarm Governor] Pre-mapped {seeded_tpp:,} Microprotein Knockout TPP Cascades (Total: {tot_tpp:,}).")
        cur.execute("""
            SELECT
                GREATEST(
                    COALESCE((SELECT MAX(total_graph_edges) FROM (SELECT total_graph_edges FROM swarm_telemetry_log ORDER BY id DESC LIMIT 10) s), 0),
                    COALESCE((SELECT reltuples::bigint FROM pg_class WHERE relname = 'orphan_discoveries'), 0)
                ) AS est_total_edges,
                (SELECT COUNT(*) FROM orphan_discoveries WHERE rf_score IS NOT NULL) AS tot_rf,
                (SELECT COUNT(*) FROM swarm_tpp_knockout_cascades) AS tot_tpp_cnt,
                COALESCE((SELECT MAX(id) FROM orphan_discoveries), 0) AS max_id;
        """)
        tot_graph, tot_rf, tot_tpp, last_max_id = cur.fetchone()

    # Persistent connection pool for up to 16 parallel workers
    worker_conns = [psycopg.connect(DB_URI) for _ in range(16)]

    try:
        while True:
            prof = resolve_governor_profile(force_mode)
            n_workers = prof["workers"]
            n_limit = prof["neighbor_limit"]
            sleep_s = prof["sleep_sec"]

            if prof["mode"] != last_mode:
                print(
                    f"\n[Swarm Governor] Regime -> {prof['mode']} | Workers={n_workers} | "
                    f"Fanout={n_limit} | Reason: {prof['reason']}"
                )
                last_mode = prof["mode"]

            steps_per_wave = n_workers * 3
            wave_t0 = time.perf_counter()
            latencies = []
            wave_edges = 0

            # Claim a deduplicated batch of the closest unexpanded frontier nodes
            nodes = []
            for _ in range(6):
                if len(nodes) >= steps_per_wave:
                    break
                more = claim_dense_frontier_batch(claim_conn, target_nodes=steps_per_wave - len(nodes))
                if not more:
                    break
                nodes.extend(more)

            # If the Pending frontier is exhausted (0 Pending), automatically jump into unvisited Dark Proteome islands
            if len(nodes) < steps_per_wave:
                seeded = seed_unexpanded_proteome_islands(claim_conn, target_nodes=steps_per_wave - len(nodes))
                if seeded:
                    nodes.extend(seeded)

            if not nodes:
                time.sleep(1.0)
                continue

            def _run_worker_chain(w_idx):
                local_lats = []
                local_edges = 0
                c = worker_conns[w_idx]
                my_nodes = nodes[w_idx::n_workers]
                for nd in my_nodes:
                    lat_ms, n_ins = expand_dense_node(c, nd, neighbor_limit=n_limit, ef_search=64)
                    local_lats.append(lat_ms)
                    local_edges += n_ins
                    if sleep_s > 0:
                        time.sleep(sleep_s)
                return local_lats, local_edges

            with ThreadPoolExecutor(max_workers=n_workers) as ex:
                futures = [ex.submit(_run_worker_chain, i) for i in range(n_workers)]
                for fut in as_completed(futures):
                    l_list, e_cnt = fut.result()
                    latencies.extend(l_list)
                    wave_edges += e_cnt

            # Incrementally score newly discovered 3D RangeFinder pairs & Microprotein TPP Knockout Cascades
            with admin_conn.cursor() as cur:
                p_add, rf_sc = expand_3d_pocket_ranges_batch(
                    cur,
                    batch_size=prof["pocket_batch"],
                    min_edge_id=last_max_id
                )
                tpp_sc = expand_tpp_microprotein_knockouts_batch(
                    cur,
                    batch_size=prof["pocket_batch"],
                    min_edge_id=last_max_id
                )
                cur.execute("SELECT COALESCE(MAX(id), %s) FROM orphan_discoveries;", (last_max_id,))
                last_max_id = cur.fetchone()[0]
                save_swarm_checkpoint(cur)

            wave_wall_s = max(time.perf_counter() - wave_t0, 1e-4)
            wave_steps = max(len(latencies), 1)
            total_steps += wave_steps
            total_edges += wave_edges
            total_pockets += p_add
            total_rf_pairs += rf_sc
            total_tpp_cascades += tpp_sc

            latencies.sort()
            avg_ms = round(sum(latencies) / max(len(latencies), 1), 2)
            p95_ms = round(latencies[int(len(latencies) * 0.95)] if latencies else 0.0, 2)
            sps = round(wave_steps / wave_wall_s, 2)
            eps = round(wave_edges / wave_wall_s, 2)
            ms_per_edge = round(avg_ms / max(wave_edges / max(wave_steps, 1), 8.0), 1)
            GOVERNOR.record_wave_metrics(wave_wall_s, wave_edges, avg_ms, wave_steps)

            tot_graph += wave_edges
            tot_rf += rf_sc
            tot_tpp += tpp_sc
            with admin_conn.cursor() as cur:
                if cycle > 0 and cycle % 50 == 0:
                    cur.execute("SELECT COUNT(*) FROM orphan_discoveries WHERE rf_score IS NOT NULL;")
                    tot_rf = cur.fetchone()[0]
                    cur.execute("SELECT COUNT(*) FROM swarm_tpp_knockout_cascades;")
                    tot_tpp = cur.fetchone()[0]
                cur.execute("""
                    INSERT INTO swarm_telemetry_log (
                        power_mode, active_workers, hid_idle_sec, load_avg_1m,
                        steps_completed, edges_discovered, pockets_3d_indexed, rf_pairs_scored,
                        tpp_cascades_mapped,
                        avg_step_ms, p95_step_ms, steps_per_sec, edges_per_sec,
                        total_graph_edges, total_3d_rf_edges, total_tpp_cascades
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
                """, (
                    prof["mode"], n_workers, prof["idle_sec"], prof["load_1m"],
                    wave_steps, wave_edges, p_add, rf_sc,
                    tpp_sc,
                    avg_ms, p95_ms, sps, eps,
                    tot_graph, tot_rf, tot_tpp
                ))

            print(
                f"[Wave #{cycle+1} | {prof['mode']} ({n_workers}W)] "
                f"+{wave_steps} fresh nodes (+{wave_edges} dense edges, +{p_add} 3D pockets, +{rf_sc} RangeFinder, +{tpp_sc} TPP KOs) | "
                f"Avg: {avg_ms} ms (P95: {p95_ms} ms) | {ms_per_edge} ms/edge ({eps} edges/s) | "
                f"Graph: {tot_graph:,} edges ({tot_rf:,} 3D RF, {tot_tpp:,} TPP KOs)"
            )

            cycle += 1
            if max_cycles > 0 and cycle >= max_cycles:
                break
    finally:
        for c in worker_conns:
            try:
                c.close()
            except Exception:
                pass
        try:
            claim_conn.close()
            admin_conn.close()
        except Exception:
            pass


if __name__ == "__main__":
    force_arg = ""
    max_c = 0
    for arg in sys.argv[1:]:
        if arg in (
            "IDLE_HYPER_16W", "IDLE_TURBO_12W", "ACTIVE_BOOST_8W", "ACTIVE_FLOOR_4W",
            "IDLE_TURBO_8W", "IDLE_RAMP_6W", "ACTIVE_USER_4W", "CONTENTION_THROTTLE_2W",
            "NIGHT_TURBO_8W", "DAY_IDLE_6W", "DAY_ACTIVE_2W", "16W", "12W", "8W", "4W"
        ):
            force_arg = arg
        elif arg.startswith("--cycles="):
            max_c = int(arg.split("=", 1)[1])

    while True:
        try:
            run_adaptive_swarm(max_cycles=max_c, force_mode=force_arg)
            if max_c > 0:
                break
        except Exception as e:
            print(f"[Swarm Governor] Recovered from error: {e}")
            time.sleep(3)
