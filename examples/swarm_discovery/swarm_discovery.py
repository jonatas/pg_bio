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

# Hardware-calibrated concurrency benchmark results on 10-core Apple Silicon (8 P-cores + 2 E-cores)
CALIBRATED_BENCHMARKS = [
    (1,  173.87, 207.13,  5.75,  57.5, 100.0, "Baseline Single-Thread"),
    (2,  184.39, 215.42, 10.49, 104.9,  91.2, "Daytime Active Mode (2 P-cores, 6 P-cores free for UI)"),
    (4,  183.24, 221.95, 20.46, 204.6,  89.0, "Balanced Half-Node Concurrency"),
    (6,  169.06, 215.36, 33.55, 335.5,  97.2, "Daytime System-Idle Sweet Spot (97.2% scaling efficiency)"),
    (8,  177.53, 208.04, 41.37, 413.7,  89.9, "Midnight–07:00 Overnight Turbo (All 8 P-cores, 7.19x speedup)"),
    (10, 176.39, 219.93, 42.22, 422.2,  73.4, "Saturation Ceiling (Spills onto 2 E-cores)"),
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
    Closed-Loop Idleness & Edge-Latency Feedback Governor (No hardcoded wall-clock hours):
    1. System Idle (`idle_sec >= 180s`) with no contention -> Switches to `8W` (`IDLE_TURBO_8W`).
    2. User Starts Interacting (`idle_sec < 180s`) -> Immediately slows down from `8W`/`6W` to `4W` (`ACTIVE_USER_4W`).
    3. Resource Contention (`ms_per_edge` increases >= 1.35x over rolling baseline or > 68 ms/edge)
       -> Slows down to `2W` (`CONTENTION_THROTTLE_2W`) to relieve CPU/I/O pressure.
    4. 1-Minute Observation Ramp (`OBSERVE_WINDOW_SEC = 60.0`)
       -> Every 60 seconds, observes the past minute's `ms_per_edge`; if contention cleared,
          steps up threads (`2W -> 4W` when user is active; `2W -> 4W -> 6W -> 8W` when idle).
    """

    IDLE_THRESHOLD_SEC = 180.0
    OBSERVE_WINDOW_SEC = 60.0
    CONTENTION_RATIO_TRIGGER = 1.45
    CONTENTION_RECOVERY_RATIO = 1.25
    HARD_MS_PER_EDGE_CEILING = 95.0
    CONTENTION_STEP_MS_FLOOR = 580.0

    WORKER_PROFILES = {
        8: {"mode": "IDLE_TURBO_8W",          "neighbor_limit": 8, "sleep_sec": 0.0,  "pocket_batch": 16},
        6: {"mode": "IDLE_RAMP_6W",           "neighbor_limit": 6, "sleep_sec": 0.01, "pocket_batch": 10},
        4: {"mode": "ACTIVE_USER_4W",         "neighbor_limit": 6, "sleep_sec": 0.02, "pocket_batch": 8},
        2: {"mode": "CONTENTION_THROTTLE_2W", "neighbor_limit": 5, "sleep_sec": 0.05, "pocket_batch": 4},
    }

    def __init__(self):
        self.current_workers = 4
        self.in_contention_recovery = False
        self.baseline_ms_per_edge = 42.0
        self.recent_ms_per_edge = []
        self.recent_step_ms = []
        self.minute_bucket_ms_per_edge = []
        self.minute_bucket_step_ms = []
        self.last_minute_eval_ts = time.monotonic()
        self.last_reason = "Initial Boot -> 4W Active Balanced"

    def record_wave_metrics(self, wave_wall_s: float, wave_edges: int, avg_step_ms: float, wave_steps: int = 12):
        # Avoid false contention spikes when a wave hits already-saturated nodes (0 new edges but fast SQL step)
        effective_edges = max(wave_edges, wave_steps * 3)
        ms_per_edge = (wave_wall_s * 1000.0) / max(effective_edges, 1)
        self.recent_ms_per_edge.append(ms_per_edge)
        self.recent_step_ms.append(avg_step_ms)
        if len(self.recent_ms_per_edge) > 6:
            self.recent_ms_per_edge.pop(0)
        if len(self.recent_step_ms) > 6:
            self.recent_step_ms.pop(0)

        self.minute_bucket_ms_per_edge.append(ms_per_edge)
        self.minute_bucket_step_ms.append(avg_step_ms)

        # Update contention-free baseline using EWMA when latency is healthy
        if wave_edges >= 4 and ms_per_edge <= self.baseline_ms_per_edge * 1.25:
            self.baseline_ms_per_edge = round(0.92 * self.baseline_ms_per_edge + 0.08 * ms_per_edge, 2)
            self.baseline_ms_per_edge = max(18.0, min(self.baseline_ms_per_edge, 65.0))

    def resolve_profile(self, force_mode: str = "", idle_override: float = None, now_mono: float = None) -> dict:
        idle_sec = get_macos_hid_idle_seconds() if idle_override is None else idle_override
        now_ts = time.monotonic() if now_mono is None else now_mono
        try:
            load_1m = round(os.getloadavg()[0], 2)
        except Exception:
            load_1m = 1.0

        if force_mode in ("IDLE_TURBO_8W", "NIGHT_TURBO_8W"):
            w = 8
            reason = f"Forced 8W Turbo (idle={idle_sec:.0f}s)"
        elif force_mode in ("IDLE_RAMP_6W", "DAY_IDLE_6W"):
            w = 6
            reason = f"Forced 6W Ramp (idle={idle_sec:.0f}s)"
        elif force_mode == "ACTIVE_USER_4W":
            w = 4
            reason = f"Forced 4W Active (idle={idle_sec:.0f}s)"
        elif force_mode in ("CONTENTION_THROTTLE_2W", "DAY_ACTIVE_2W"):
            w = 2
            reason = f"Forced 2W Throttle (idle={idle_sec:.0f}s)"
        else:
            is_idle = (idle_sec >= self.IDLE_THRESHOLD_SEC)
            max_allowed_workers = 8 if is_idle else 4

            rolling_mpe = (
                sum(self.recent_ms_per_edge) / len(self.recent_ms_per_edge)
                if self.recent_ms_per_edge else self.baseline_ms_per_edge
            )
            rolling_step_ms = (
                sum(self.recent_step_ms) / len(self.recent_step_ms)
                if self.recent_step_ms else 280.0
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

            # 1. Immediate downshift to 4W when user starts using the workstation
            if not is_idle and self.current_workers > 4:
                self.current_workers = 4
                self.last_reason = (
                    f"User Active (idle {idle_sec:.0f}s < {self.IDLE_THRESHOLD_SEC:.0f}s) "
                    f"-> Immediate Downshift to 4W"
                )

            # 2. Immediate downshift to 2W if resources compete and time-per-new-edge spikes
            if is_competing and self.current_workers > 2:
                self.current_workers = 2
                self.in_contention_recovery = True
                self.last_minute_eval_ts = now_ts
                self.minute_bucket_ms_per_edge.clear()
                self.minute_bucket_step_ms.clear()
                self.last_reason = (
                    f"Resource Contention ({rolling_mpe:.1f} ms/edge > {contention_threshold:.1f} ms/edge "
                    f"[base {self.baseline_ms_per_edge:.1f} ms]) -> Throttled to 2W"
                )

            # 3. Direct switch to 8W when workstation becomes idle and NOT in contention recovery
            elif is_idle and not self.in_contention_recovery and self.current_workers < 8:
                self.current_workers = 8
                self.last_reason = (
                    f"System Idle ({idle_sec:.0f}s >= {self.IDLE_THRESHOLD_SEC:.0f}s, "
                    f"{rolling_mpe:.1f} ms/edge healthy) -> Switched to 8W Turbo"
                )

            # 4. 1-Minute Observation Window: climb threads step-by-step (2 -> 4 -> 6 -> 8)
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
                    if self.current_workers == 2:
                        self.current_workers = 4
                    elif self.current_workers == 4:
                        self.current_workers = 6 if is_idle else 4
                    elif self.current_workers == 6:
                        self.current_workers = 8 if is_idle else 4

                    if self.current_workers >= max_allowed_workers:
                        self.in_contention_recovery = False

                    self.last_reason = (
                        f"1-Min Observation Clear ({min_mpe:.1f} ms/edge <= {recovery_limit:.1f} ms/edge) "
                        f"-> Stepped Up {prev_w}W -> {self.current_workers}W"
                    )
                elif min_mpe > contention_threshold and self.current_workers > 2:
                    self.current_workers = 2
                    self.in_contention_recovery = True
                    self.last_reason = (
                        f"1-Min Observation Contention ({min_mpe:.1f} ms/edge) -> Held/Dropped to 2W"
                    )

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


def expand_3d_pocket_ranges_batch(cur, batch_size: int = 12) -> tuple[int, int]:
    """
    Finds high-priority hub or partner proteins in `orphan_discoveries` that are not yet in
    `swarm_pocket_ranges`, builds/loads their 3D coordinates, computes RangeFinder 9-property
    pocket intervals (`predict_pocket_ranges`, `pocket_numranges`, `pocket_interval_bounds`),
    and scores all their incident 1-to-N edges in `orphan_discoveries`.
    """
    cur.execute("""
        WITH missing_3d AS (
            SELECT DISTINCT d.orphan_id AS uid
            FROM orphan_discoveries d
            JOIN swarm_pocket_ranges pb ON pb.uniprot_id = d.bait_id
            LEFT JOIN swarm_pocket_ranges po ON po.uniprot_id = d.orphan_id
            WHERE po.uniprot_id IS NULL
              AND d.bait_id <> d.orphan_id
              AND d.vector_distance <= 0.14
            LIMIT %s
        )
        SELECT p.uniprot_id, p.sequence
        FROM missing_3d m
        JOIN proteins p ON p.uniprot_id = m.uid
        WHERE length(p.sequence) BETWEEN 50 AND 900;
    """, (batch_size,))
    cand_rows = cur.fetchall()
    pockets_added = 0

    for uid, seq in cand_rows:
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
            continue

        # Compute catalytic cleft centroid from CA atoms
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
            continue
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
                aromatic_probability
            )
            SELECT
                %s, %s, %s, %s, 10.0, %s, r.pk,
                nr.mw_range, nr.logp_range, nr.hbd_range, nr.hba_range, nr.tpsa_range, nr.rot_range, nr.charge_range, nr.arom_range, nr.fsp3_range,
                b.mw_min, b.mw_max, b.logp_min, b.logp_max, b.hbd_min, b.hbd_max, b.hba_min, b.hba_max,
                b.tpsa_min, b.tpsa_max, b.rot_min, b.rot_max, b.charge_min, b.charge_max, b.arom_min, b.arom_max, b.fsp3_min, b.fsp3_max,
                b.aromatic_probability
            FROM raw_pk r,
                 LATERAL pocket_numranges(r.pk) nr,
                 LATERAL pocket_interval_bounds(r.pk) b
            ON CONFLICT (uniprot_id) DO NOTHING;
        """, (uid, pdb_str, cx, cy, cz, uid, cx, cy, cz, pdb_path))
        pockets_added += cur.rowcount

    # Score any newly completed 3D pairs in orphan_discoveries where both proteins have real 3D pockets
    cur.execute("""
        WITH targets AS (
            SELECT od.id, od.bait_id, od.orphan_id, od.family_name, od.confidence_score, od.vector_distance
            FROM orphan_discoveries od
            JOIN swarm_pocket_ranges pb ON pb.uniprot_id = od.bait_id AND pb.pdb_path IS NOT NULL
            JOIN swarm_pocket_ranges po ON po.uniprot_id = od.orphan_id AND po.pdb_path IS NOT NULL
            WHERE od.rf_score IS NULL
              AND od.bait_id <> od.orphan_id
            FOR UPDATE OF od SKIP LOCKED
            LIMIT %s
        )
        UPDATE orphan_discoveries od
        SET rf_range_iou = ROUND(pocket_range_iou(pb.pocket, po.pocket, false)::numeric, 4),
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
                    0.50 * pocket_range_iou(pb.pocket, po.pocket, false) +
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
                ELSE NULL
            END
        FROM targets t
        JOIN swarm_pocket_ranges pb ON pb.uniprot_id = t.bait_id
        JOIN swarm_pocket_ranges po ON po.uniprot_id = t.orphan_id
        LEFT JOIN swarm_enzyme_substrates s ON LOWER(s.family_name) = LOWER(t.family_name)
        WHERE od.id = t.id;
    """, (batch_size * 25,))
    rf_scored = cur.rowcount
    return pockets_added, rf_scored


def worker_step(conn, neighbor_limit: int = 6) -> tuple[float, int]:
    """
    Executes one lock-free (`FOR UPDATE SKIP LOCKED`) Swarm expansion step using a single
    set-based SQL query inside PostgreSQL. Returns `(elapsed_ms, new_edges_inserted)`.
    """
    t0 = time.perf_counter()
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE orphan_discoveries
            SET status = 'Exploring'
            WHERE id = (
                SELECT id
                FROM orphan_discoveries
                WHERE status = 'Pending'
                ORDER BY vector_distance ASC, id ASC
                FOR UPDATE SKIP LOCKED
                LIMIT 1
            )
            RETURNING id, COALESCE(family_name, 'Dark Proteome Island'), orphan_id, COALESCE(orphan_organism, 'Archaea sp.');
        """)
        row = cur.fetchone()
        conn.commit()

        if not row:
            # Seed a new unvisited island in a single SQL statement
            cur.execute("SET LOCAL hnsw.ef_search = 120;")
            cur.execute("""
                WITH seed AS (
                    SELECT p.uniprot_id, p.name, p.sequence, p.embedding,
                           COALESCE(NULLIF(split_part(split_part(p.name, 'OS=', 2), ' OX=', 1), ''), 'Archaea sp.') AS org
                    FROM proteins p
                    WHERE p.name LIKE '%%OS=%%'
                      AND length(p.sequence) BETWEEN 70 AND 440
                      AND NOT EXISTS (
                          SELECT 1 FROM orphan_discoveries od
                          WHERE od.bait_id = p.uniprot_id
                      )
                    OFFSET floor(random() * 300)::int
                    LIMIT 1
                ),
                raw_neighbors AS (
                    SELECT
                        s.uniprot_id AS b_id,
                        s.org AS b_org,
                        s.sequence AS b_seq,
                        s.embedding AS b_emb,
                        p.uniprot_id AS o_id,
                        COALESCE(NULLIF(split_part(split_part(p.name, 'OS=', 2), ' OX=', 1), ''), 'Archaea sp.') AS o_org,
                        p.sequence AS o_seq,
                        p.embedding AS o_emb,
                        (p.embedding <=> s.embedding) AS dist
                    FROM seed s
                    JOIN LATERAL (
                        SELECT p2.uniprot_id, p2.name, p2.sequence, p2.embedding
                        FROM proteins p2
                        WHERE p2.uniprot_id <> s.uniprot_id
                          AND p2.name LIKE '%%OS=%%'
                        ORDER BY p2.embedding <=> s.embedding ASC
                        LIMIT 16
                    ) p ON TRUE
                    WHERE (p.embedding <=> s.embedding) <= 0.14
                    LIMIT %s
                ),
                scored AS (
                    SELECT
                        rn.*,
                        interface_complementarity(rn.o_emb::vector::real[], rn.b_emb::vector::real[]) AS icomp,
                        LEAST(GREATEST(
                            smith_waterman_score(LEFT(rn.b_seq, 500), LEFT(rn.o_seq, 500), 2, -1, -2)::double precision
                            / GREATEST(2.0 * LEAST(length(rn.b_seq), length(rn.o_seq), 500)::double precision, 1.0),
                            0.0
                        ), 1.0) AS sw_sim,
                        LEAST(length(rn.b_seq), length(rn.o_seq))::double precision
                            / GREATEST(GREATEST(length(rn.b_seq), length(rn.o_seq))::double precision, 1.0) AS len_ratio,
                        sequence_entropy(rn.o_seq) AS ent,
                        chaperone_lineage_hash('SwarmSeed', rn.b_id, rn.o_seq, 1) AS lineage_id
                    FROM raw_neighbors rn
                )
                INSERT INTO orphan_discoveries (
                    family_name, bait_id, bait_organism, orphan_id, orphan_organism,
                    vector_distance, status, confidence_score, seq_identity_pct, twilight_score, chaperone_lineage_id
                )
                SELECT
                    'Dark Proteome Island',
                    sc.b_id,
                    sc.b_org,
                    sc.o_id,
                    sc.o_org,
                    sc.dist,
                    'Pending',
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
                FROM scored sc
                WHERE NOT EXISTS (
                    SELECT 1 FROM orphan_discoveries od
                    WHERE od.bait_id = sc.b_id AND od.orphan_id = sc.o_id
                );
            """, (neighbor_limit,))
            seeded = cur.rowcount
            conn.commit()
            return (time.perf_counter() - t0) * 1000.0, seeded

        edge_id, family_name, current_node, cached_parent_org = row
        try:
            cur.execute("SET LOCAL hnsw.ef_search = 120;")
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
                        LIMIT 16
                    ) p ON TRUE
                    WHERE (p.embedding <=> pp.embedding) <= 0.15
                      AND NOT EXISTS (
                          SELECT 1 FROM orphan_discoveries od
                          WHERE od.bait_id = pp.uniprot_id AND od.orphan_id = p.uniprot_id
                      )
                    ORDER BY (p.embedding <=> pp.embedding) ASC
                    LIMIT %s
                ),
                scored AS (
                    SELECT
                        rn.*,
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
                    'Pending',
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
            """, (cached_parent_org, current_node, neighbor_limit, family_name))
            inserted = cur.rowcount

            cur.execute("UPDATE orphan_discoveries SET status = 'Explored' WHERE id = %s;", (edge_id,))
            conn.commit()
            return (time.perf_counter() - t0) * 1000.0, inserted
        except Exception:
            conn.rollback()
            cur.execute("UPDATE orphan_discoveries SET status = 'Pending' WHERE id = %s;", (edge_id,))
            conn.commit()
            raise


def run_adaptive_swarm(max_cycles: int = 0, force_mode: str = ""):
    admin_conn = psycopg.connect(DB_URI, autocommit=True)
    with admin_conn.cursor() as cur:
        ensure_swarm_indexes(cur)
        cur.execute("UPDATE orphan_discoveries SET status = 'Pending' WHERE status = 'Exploring';")
        if cur.rowcount > 0:
            print(f"[Swarm Governor] Recovered {cur.rowcount} stranded 'Exploring' nodes back to 'Pending'.")

    print("================================================================================")
    print("pg_bio Closed-Loop Idleness & Edge-Latency Feedback Governor (8W -> 4W -> 2W + 1m Ramp)")
    print("================================================================================")

    cycle = 0
    total_steps = 0
    total_edges = 0
    total_pockets = 0
    total_rf_pairs = 0
    last_mode = None

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

        # Execute a parallel wave of steps across `n_workers` connections
        steps_per_wave = n_workers * 3
        conns = [psycopg.connect(DB_URI) for _ in range(n_workers)]
        wave_t0 = time.perf_counter()
        latencies = []
        wave_edges = 0

        try:
            def _run_worker_chain(w_idx):
                local_lats = []
                local_edges = 0
                c = conns[w_idx]
                for _ in range(3):
                    lat_ms, n_ins = worker_step(c, n_limit)
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
        finally:
            for c in conns:
                try:
                    c.close()
                except Exception:
                    pass

        # Expand 3D RangeFinder pockets & score 1-to-N pairs
        with admin_conn.cursor() as cur:
            p_add, rf_sc = expand_3d_pocket_ranges_batch(cur, prof["pocket_batch"])

        wave_wall_s = max(time.perf_counter() - wave_t0, 1e-4)
        wave_steps = len(latencies)
        total_steps += wave_steps
        total_edges += wave_edges
        total_pockets += p_add
        total_rf_pairs += rf_sc

        latencies.sort()
        avg_ms = round(sum(latencies) / max(len(latencies), 1), 2)
        p95_ms = round(latencies[int(len(latencies) * 0.95)] if latencies else 0.0, 2)
        sps = round(wave_steps / wave_wall_s, 2)
        eps = round(wave_edges / wave_wall_s, 2)
        ms_per_edge = round((wave_wall_s * 1000.0) / max(wave_edges, wave_steps * 3, 1), 1)
        GOVERNOR.record_wave_metrics(wave_wall_s, wave_edges, avg_ms, wave_steps)

        with admin_conn.cursor() as cur:
            cur.execute("SELECT COUNT(*), COUNT(*) FILTER (WHERE rf_score IS NOT NULL) FROM orphan_discoveries;")
            tot_graph, tot_rf = cur.fetchone()
            cur.execute("""
                INSERT INTO swarm_telemetry_log (
                    power_mode, active_workers, hid_idle_sec, load_avg_1m,
                    steps_completed, edges_discovered, pockets_3d_indexed, rf_pairs_scored,
                    avg_step_ms, p95_step_ms, steps_per_sec, edges_per_sec,
                    total_graph_edges, total_3d_rf_edges
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
            """, (
                prof["mode"], n_workers, prof["idle_sec"], prof["load_1m"],
                wave_steps, wave_edges, p_add, rf_sc,
                avg_ms, p95_ms, sps, eps,
                tot_graph, tot_rf
            ))

        print(
            f"[Wave #{cycle+1} | {prof['mode']} ({n_workers}W)] "
            f"+{wave_steps} nodes (+{wave_edges} edges, +{p_add} 3D pockets, +{rf_sc} RangeFinder pairs) | "
            f"Avg: {avg_ms} ms (P95: {p95_ms} ms) | {ms_per_edge} ms/edge ({eps} edges/s) | "
            f"Graph: {tot_graph:,} edges ({tot_rf:,} 3D RangeFinder)"
        )

        cycle += 1
        if max_cycles > 0 and cycle >= max_cycles:
            break


if __name__ == "__main__":
    force_arg = ""
    max_c = 0
    for arg in sys.argv[1:]:
        if arg in (
            "IDLE_TURBO_8W", "IDLE_RAMP_6W", "ACTIVE_USER_4W", "CONTENTION_THROTTLE_2W",
            "NIGHT_TURBO_8W", "DAY_IDLE_6W", "DAY_ACTIVE_2W"
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
