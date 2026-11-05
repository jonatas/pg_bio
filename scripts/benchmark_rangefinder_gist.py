import json
import math
import time
import psycopg

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"


def build_6n9q_like_pocket_pdb() -> str:
    """
    Constructs a realistic 3D all-atom binding pocket (~140 atoms within 10 A of center)
    modeled after N-acyl-homoserine-lactone lactonase (PDB 6N9Q, Figure 1 of arXiv:2611.03379v1).
    """
    residues = [
        ("PHE", "C"), ("ASP", "O"), ("LEU", "C"), ("TRP", "C"), ("HIS", "N"),
        ("TYR", "O"), ("VAL", "C"), ("SER", "O"), ("ILE", "C"), ("ARG", "N"),
        ("GLU", "O"), ("ALA", "C"), ("PHE", "C"), ("THR", "O"), ("MET", "S"),
        ("ASN", "N"), ("LEU", "C"), ("HIS", "N"), ("TYR", "O"), ("VAL", "C"),
    ]
    lines = []
    serial = 1
    for idx, (res_name, sc_elem) in enumerate(residues):
        theta = idx * 0.65
        phi = 0.35 + (idx / len(residues)) * 2.4
        r = 5.2 + 2.8 * ((idx % 3) / 2.0)
        cx = r * math.cos(theta) * math.sin(phi)
        cy = r * math.sin(theta) * math.sin(phi)
        cz = r * math.cos(phi)

        atoms = [
            ("N", -0.75, 0.0, 0.0, "N"),
            ("CA", 0.0, 0.0, 0.0, "C"),
            ("C", 0.78, 0.0, 0.0, "C"),
            ("O", 1.35, 0.85, 0.0, "O"),
            ("CB", -0.35, -1.20, 0.45, "C"),
            ("CG", -0.85, -2.10, 0.90, "C"),
            (f"{sc_elem}D1", -1.35, -2.85, 1.25, sc_elem),
        ]
        for aname, dx, dy, dz, elem in atoms:
            lines.append(
                f"ATOM  {serial:5d}  {aname:<3s} {res_name:3s} A{idx + 10:4d}    "
                f"{cx + dx:8.3f}{cy + dy:8.3f}{cz + dz:8.3f}  1.00 92.40           {elem:>2s}  "
            )
            serial += 1
    return "\n".join(lines)


def run_benchmarks():
    print("=" * 88)
    print("pg_bio RangeFinder (arXiv:2611.03379v1) GiST vs. B-Tree Benchmark Suite")
    print("=" * 88)

    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            # ------------------------------------------------------------------
            # 1. Single-Pocket Inference Latency & SE(3) Invariance Verification
            # ------------------------------------------------------------------
            pdb_text = build_6n9q_like_pocket_pdb()
            # Warmup
            cur.execute(
                "SELECT predict_pocket_ranges('6N9Q', %s, 0.0, 0.0, 0.0, 10.0)::text;",
                (pdb_text,),
            )
            pocket_json = json.loads(cur.fetchone()[0])

            n_iters = 200
            t0 = time.perf_counter()
            for _ in range(n_iters):
                cur.execute(
                    "SELECT predict_pocket_ranges('6N9Q', %s, 0.0, 0.0, 0.0, 10.0);",
                    (pdb_text,),
                )
                cur.fetchone()
            pocket_ms = ((time.perf_counter() - t0) * 1000.0) / n_iters

            print(f"\n[1] Pocket-to-Interval Prediction Latency (10 A pocket, {pocket_json['atom_count']} atoms, {pocket_json['edge_count_5a']} 5A edges):")
            print(f"    - pg_bio Native Rust (CPU single-core): {pocket_ms:.3f} ms / pocket")
            print(f"    - RangeFinder SE(3)-Transformer (A10G GPU, Table 1): 12.200 ms / pocket")
            print(f"    - Foldseek Baseline (Table 1):                      127.000 ms / pocket")
            print(f"    - APoc Baseline (Table 1):                       26,300.000 ms / pocket")

            # ------------------------------------------------------------------
            # 2. Forward Virtual Screening: 1 Pocket -> 250,000 Molecules
            # ------------------------------------------------------------------
            print("\n[2] Populating 250,000-Compound Drug-Like Library & Building Indexes...")
            cur.execute("DROP TABLE IF EXISTS rf_bench_molecules CASCADE;")
            t_pop0 = time.perf_counter()
            cur.execute(
                """
                CREATE UNLOGGED TABLE rf_bench_molecules AS
                SELECT
                    g AS mol_id,
                    'MOL_' || g AS chembl_id,
                    (120.0 + (hashtext('mw_' || g)::bigint & 65535) / 65535.0 * 580.0)::float8 AS mw,
                    (-2.0 + (hashtext('lp_' || g)::bigint & 65535) / 65535.0 * 8.5)::float8 AS logp,
                    ((hashtext('hbd_' || g)::bigint & 7))::float8 AS hbd,
                    (1 + (hashtext('hba_' || g)::bigint & 11))::float8 AS hba,
                    (10.0 + (hashtext('tpsa_' || g)::bigint & 65535) / 65535.0 * 170.0)::float8 AS tpsa,
                    ((hashtext('rot_' || g)::bigint & 11))::float8 AS rotatable_bonds,
                    (CASE WHEN (g % 10) = 0 THEN 1.0 WHEN (g % 15) = 0 THEN -1.0 ELSE 0.0 END)::float8 AS formal_charge,
                    ((hashtext('arom_' || g)::bigint & 3))::float8 AS aromatic_rings,
                    ((hashtext('fsp3_' || g)::bigint & 65535) / 65535.0)::float8 AS fsp3
                FROM generate_series(1, 250000) AS g;
                """
            )
            pop_ms = (time.perf_counter() - t_pop0) * 1000.0

            t_bt0 = time.perf_counter()
            cur.execute(
                """
                CREATE INDEX idx_rf_mol_btree
                ON rf_bench_molecules (mw, logp, tpsa, hba);
                """
            )
            btree_build_ms = (time.perf_counter() - t_bt0) * 1000.0

            cur.execute("ANALYZE rf_bench_molecules;")
            cur.execute("SELECT pg_relation_size('idx_rf_mol_btree') / (1024.0 * 1024.0);")
            btree_mb = float(cur.fetchone()[0])
            print(f"    - 250k Table Population: {pop_ms:.1f} ms | B-Tree Index Build: {btree_build_ms:.1f} ms ({btree_mb:.2f} MB)")

            # Target pocket intervals (typical kinase/NSAID pocket capturing ~2.5% of chemical space -> ~40x reduction)
            q_params = {
                "mw_min": 280.0, "mw_max": 420.0,
                "logp_min": 1.5, "logp_max": 4.2,
                "hbd_min": 0.0, "hbd_max": 3.0,
                "hba_min": 2.0, "hba_max": 7.0,
                "tpsa_min": 40.0, "tpsa_max": 105.0,
                "rot_min": 1.0, "rot_max": 6.0,
                "arom_min": 1.0, "arom_max": 3.0,
                "fsp3_min": 0.15, "fsp3_max": 0.65,
            }

            # 2A. Python Client-Side Loop over 250,000 molecules
            t_py0 = time.perf_counter()
            cur.execute(
                "SELECT mol_id, mw, logp, hbd, hba, tpsa, rotatable_bonds, aromatic_rings, fsp3 FROM rf_bench_molecules;"
            )
            rows = cur.fetchall()
            py_hits = [
                r[0] for r in rows
                if (q_params["mw_min"] <= r[1] <= q_params["mw_max"]
                    and q_params["logp_min"] <= r[2] <= q_params["logp_max"]
                    and q_params["hbd_min"] <= r[3] <= q_params["hbd_max"]
                    and q_params["hba_min"] <= r[4] <= q_params["hba_max"]
                    and q_params["tpsa_min"] <= r[5] <= q_params["tpsa_max"]
                    and q_params["rot_min"] <= r[6] <= q_params["rot_max"]
                    and q_params["arom_min"] <= r[7] <= q_params["arom_max"]
                    and q_params["fsp3_min"] <= r[8] <= q_params["fsp3_max"])
            ]
            py_screen_ms = (time.perf_counter() - t_py0) * 1000.0

            # 2B. PostgreSQL Sequential Scan with numrange @> numeric (BCD conversion overhead)
            cur.execute("SET enable_indexscan = off; SET enable_bitmapscan = off;")
            t_nr0 = time.perf_counter()
            for _ in range(5):
                cur.execute(
                    """
                    SELECT count(*) FROM rf_bench_molecules
                    WHERE numrange(%(mw_min)s::numeric, %(mw_max)s::numeric, '[]') @> mw::numeric
                      AND numrange(%(logp_min)s::numeric, %(logp_max)s::numeric, '[]') @> logp::numeric
                      AND numrange(%(hbd_min)s::numeric, %(hbd_max)s::numeric, '[]') @> hbd::numeric
                      AND numrange(%(hba_min)s::numeric, %(hba_max)s::numeric, '[]') @> hba::numeric
                      AND numrange(%(tpsa_min)s::numeric, %(tpsa_max)s::numeric, '[]') @> tpsa::numeric
                      AND numrange(%(rot_min)s::numeric, %(rot_max)s::numeric, '[]') @> rotatable_bonds::numeric
                      AND numrange(%(arom_min)s::numeric, %(arom_max)s::numeric, '[]') @> aromatic_rings::numeric
                      AND numrange(%(fsp3_min)s::numeric, %(fsp3_max)s::numeric, '[]') @> fsp3::numeric;
                    """,
                    q_params,
                )
                nr_count = cur.fetchone()[0]
            seq_numrange_ms = ((time.perf_counter() - t_nr0) * 1000.0) / 5.0

            # 2C. PostgreSQL Sequential Scan with IEEE-754 float8 BETWEEN
            t_fl0 = time.perf_counter()
            for _ in range(5):
                cur.execute(
                    """
                    SELECT count(*) FROM rf_bench_molecules
                    WHERE mw BETWEEN %(mw_min)s AND %(mw_max)s
                      AND logp BETWEEN %(logp_min)s AND %(logp_max)s
                      AND hbd BETWEEN %(hbd_min)s AND %(hbd_max)s
                      AND hba BETWEEN %(hba_min)s AND %(hba_max)s
                      AND tpsa BETWEEN %(tpsa_min)s AND %(tpsa_max)s
                      AND rotatable_bonds BETWEEN %(rot_min)s AND %(rot_max)s
                      AND aromatic_rings BETWEEN %(arom_min)s AND %(arom_max)s
                      AND fsp3 BETWEEN %(fsp3_min)s AND %(fsp3_max)s;
                    """,
                    q_params,
                )
                fl_count = cur.fetchone()[0]
            seq_float_ms = ((time.perf_counter() - t_fl0) * 1000.0) / 5.0

            # 2D. PostgreSQL Composite B-Tree Index Scan
            cur.execute("SET enable_indexscan = on; SET enable_bitmapscan = on;")
            t_idx0 = time.perf_counter()
            for _ in range(20):
                cur.execute(
                    """
                    SELECT count(*) FROM rf_bench_molecules
                    WHERE mw BETWEEN %(mw_min)s AND %(mw_max)s
                      AND logp BETWEEN %(logp_min)s AND %(logp_max)s
                      AND hbd BETWEEN %(hbd_min)s AND %(hbd_max)s
                      AND hba BETWEEN %(hba_min)s AND %(hba_max)s
                      AND tpsa BETWEEN %(tpsa_min)s AND %(tpsa_max)s
                      AND rotatable_bonds BETWEEN %(rot_min)s AND %(rot_max)s
                      AND aromatic_rings BETWEEN %(arom_min)s AND %(arom_max)s
                      AND fsp3 BETWEEN %(fsp3_min)s AND %(fsp3_max)s;
                    """,
                    q_params,
                )
                bt_count = cur.fetchone()[0]
            btree_scan_ms = ((time.perf_counter() - t_idx0) * 1000.0) / 20.0

            assert len(py_hits) == nr_count == fl_count == bt_count
            reduction_250k = 250000.0 / max(1, bt_count)
            print(f"\n[2] Forward Virtual Screening Results (1 Pocket -> 250,000 Molecules | {bt_count} hits, {reduction_250k:.1f}x reduction):")
            print(f"    - Python Client-Side Pull + Filter:            {py_screen_ms:8.2f} ms (1.0x baseline)")
            print(f"    - PG Seq Scan (`numrange @> numeric` BCD):     {seq_numrange_ms:8.2f} ms ({py_screen_ms / seq_numrange_ms:5.1f}x speedup)")
            print(f"    - PG Seq Scan (`float8 BETWEEN` IEEE-754):     {seq_float_ms:8.2f} ms ({py_screen_ms / seq_float_ms:5.1f}x speedup)")
            print(f"    - PG Composite B-Tree Index Scan:              {btree_scan_ms:8.2f} ms ({py_screen_ms / btree_scan_ms:5.1f}x speedup)")

            # ------------------------------------------------------------------
            # 3. Reverse Target Fishing (Polypharmacology): 1 Ligand -> 20,000 Pockets
            #    Here each row is a POCKET with 9 `numrange` intervals [mu - 2|sigma|, mu + 2|sigma|]!
            #    This is where Multi-Column GiST on `numrange` dominates B-Tree!
            # ------------------------------------------------------------------
            print("\n[3] Populating 20,000 Protein Binding Pockets with 9 RangeFinder `numrange` Intervals...")
            cur.execute("DROP TABLE IF EXISTS rf_bench_pockets CASCADE;")
            cur.execute(
                """
                CREATE UNLOGGED TABLE rf_bench_pockets AS
                SELECT
                    p AS pocket_id,
                    'UNIPROT_P' || p AS uniprot_id,
                    mw_c - mw_w AS mw_min, mw_c + mw_w AS mw_max,
                    numrange((mw_c - mw_w)::numeric, (mw_c + mw_w)::numeric, '[]') AS mw_range,
                    numrange((lp_c - lp_w)::numeric, (lp_c + lp_w)::numeric, '[]') AS logp_range,
                    numrange((hbd_c - 1.5)::numeric, (hbd_c + 1.5)::numeric, '[]') AS hbd_range,
                    numrange((hba_c - 2.0)::numeric, (hba_c + 2.0)::numeric, '[]') AS hba_range,
                    numrange((tpsa_c - 28.0)::numeric, (tpsa_c + 28.0)::numeric, '[]') AS tpsa_range,
                    numrange((rot_c - 2.0)::numeric, (rot_c + 2.0)::numeric, '[]') AS rot_range,
                    numrange((-0.08)::numeric, (0.08)::numeric, '[]') AS charge_range,
                    numrange((arom_c - 1.0)::numeric, (arom_c + 1.0)::numeric, '[]') AS arom_range,
                    numrange(GREATEST(0.0, fsp3_c - 0.18)::numeric, LEAST(1.0, fsp3_c + 0.18)::numeric, '[]') AS fsp3_range
                FROM (
                    SELECT
                        p,
                        180.0 + (hashtext('pmw_' || p)::bigint & 65535) / 65535.0 * 480.0 AS mw_c,
                        45.0 + (hashtext('pmww_' || p)::bigint & 65535) / 65535.0 * 45.0 AS mw_w,
                        -1.0 + (hashtext('plp_' || p)::bigint & 65535) / 65535.0 * 6.5 AS lp_c,
                        0.9 + (hashtext('plpw_' || p)::bigint & 65535) / 65535.0 * 0.8 AS lp_w,
                        1.0 + (hashtext('phbd_' || p)::bigint & 3) AS hbd_c,
                        3.0 + (hashtext('phba_' || p)::bigint & 7) AS hba_c,
                        35.0 + (hashtext('ptpsa_' || p)::bigint & 65535) / 65535.0 * 130.0 AS tpsa_c,
                        2.0 + (hashtext('prot_' || p)::bigint & 7) AS rot_c,
                        1.0 + (hashtext('parom_' || p)::bigint & 3) AS arom_c,
                        0.15 + (hashtext('pfsp3_' || p)::bigint & 65535) / 65535.0 * 0.70 AS fsp3_c
                    FROM generate_series(1, 20000) AS p
                ) sub;
                """
            )

            t_g0 = time.perf_counter()
            cur.execute(
                """
                CREATE INDEX idx_rf_pockets_gist
                ON rf_bench_pockets
                USING gist (mw_range, logp_range, tpsa_range, hba_range, hbd_range, rot_range, arom_range, fsp3_range);
                """
            )
            gist_build_ms = (time.perf_counter() - t_g0) * 1000.0

            t_pbt0 = time.perf_counter()
            cur.execute(
                """
                CREATE INDEX idx_rf_pockets_btree
                ON rf_bench_pockets (mw_min, mw_max);
                """
            )
            pbtree_build_ms = (time.perf_counter() - t_pbt0) * 1000.0
            cur.execute("ANALYZE rf_bench_pockets;")
            cur.execute("SELECT pg_relation_size('idx_rf_pockets_gist') / (1024.0 * 1024.0);")
            gist_mb = float(cur.fetchone()[0])

            # Query ligand: Ibuprofen ('CC(C)Cc1ccc(C(C)C(=O)O)cc1')
            query_smiles = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"

            # 3A. Sequential Scan over 20,000 pockets
            cur.execute("SET enable_indexscan = off; SET enable_bitmapscan = off;")
            t_pseq0 = time.perf_counter()
            for _ in range(10):
                cur.execute(
                    """
                    WITH lig AS (SELECT * FROM compute_ligand_properties_row(%s))
                    SELECT count(*) FROM rf_bench_pockets p, lig
                    WHERE p.mw_range @> lig.mw::numeric
                      AND p.logp_range @> lig.logp::numeric
                      AND p.tpsa_range @> lig.tpsa::numeric
                      AND p.hba_range @> lig.hba::numeric
                      AND p.hbd_range @> lig.hbd::numeric
                      AND p.rot_range @> lig.rotatable_bonds::numeric
                      AND p.arom_range @> lig.aromatic_rings::numeric
                      AND p.fsp3_range @> lig.fsp3::numeric;
                    """,
                    (query_smiles,),
                )
                p_seq_hits = cur.fetchone()[0]
            pocket_seq_ms = ((time.perf_counter() - t_pseq0) * 1000.0) / 10.0

            # 3B. Multi-Column GiST Index Scan over 20,000 pockets
            cur.execute("SET enable_indexscan = on; SET enable_bitmapscan = on;")
            t_pgist0 = time.perf_counter()
            for _ in range(25):
                cur.execute(
                    """
                    WITH lig AS (SELECT * FROM compute_ligand_properties_row(%s))
                    SELECT count(*) FROM rf_bench_pockets p, lig
                    WHERE p.mw_range @> lig.mw::numeric
                      AND p.logp_range @> lig.logp::numeric
                      AND p.tpsa_range @> lig.tpsa::numeric
                      AND p.hba_range @> lig.hba::numeric
                      AND p.hbd_range @> lig.hbd::numeric
                      AND p.rot_range @> lig.rotatable_bonds::numeric
                      AND p.arom_range @> lig.aromatic_rings::numeric
                      AND p.fsp3_range @> lig.fsp3::numeric;
                    """,
                    (query_smiles,),
                )
                p_gist_hits = cur.fetchone()[0]
            pocket_gist_ms = ((time.perf_counter() - t_pgist0) * 1000.0) / 25.0

            assert p_seq_hits == p_gist_hits
            print(f"    - 8-Column GiST Index Build (20k pockets):     {gist_build_ms:8.1f} ms ({gist_mb:.2f} MB)")
            print(f"    - Sequential Scan (1 Ligand -> 20k Pockets):   {pocket_seq_ms:8.2f} ms ({p_seq_hits} compatible pockets)")
            print(f"    - Multi-Column GiST Index Scan (`@>`):         {pocket_gist_ms:8.2f} ms ({pocket_seq_ms / pocket_gist_ms:5.1f}x speedup over Seq Scan)")

            # ------------------------------------------------------------------
            # 4. Formal Charge Degeneracy Audit (8-Property vs. 9-Property Range-IoU & RRF)
            # ------------------------------------------------------------------
            print("\n[4] Auditing Appendix A.1 Formal Charge Degeneracy (8-Property vs. 9-Property Filtering):")
            cur.execute(
                """
                WITH sample_binders AS (
                    SELECT ARRAY[
                        'CC(=O)Oc1ccccc1C(=O)O',
                        'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
                        'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                        'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl'
                    ] AS binders,
                    ARRAY[
                        'CC(=O)Oc1ccccc1C(=O)O',
                        'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
                        'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                        'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl',
                        'C1CCCCC1',
                        'CCCCCCCCCCCCCCCCCC(=O)O',
                        'NCC(=O)O',
                        'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
                        'c1ccc2ccccc2c1',
                        'CC(=O)Nc1ccc(O)cc1'
                    ] AS bg_lib
                ),
                pred_pocket AS (
                    SELECT predict_pocket_ranges('6N9Q_active', %s, 0.0, 0.0, 0.0, 10.0) AS p
                )
                SELECT
                    evaluate_rangefinder_screening((SELECT p FROM pred_pocket), binders, bg_lib, true)::text AS eval_9,
                    evaluate_rangefinder_screening((SELECT p FROM pred_pocket), binders, bg_lib, false)::text AS eval_8
                FROM sample_binders;
                """,
                (pdb_text,),
            )
            ev9_str, ev8_str = cur.fetchone()
            ev9 = json.loads(ev9_str)
            ev8 = json.loads(ev8_str)
            print(f"    - 9-Property (with Formal Charge): Range-IoU = {ev9['range_iou']:.3f} | Recall = {ev9['recall']:.3f} | Reduction = {ev9['reduction_factor']:.2f}x | RRF = {ev9['rrf_score']:.3f}")
            print(f"    - 8-Property (no Formal Charge):   Range-IoU = {ev8['range_iou']:.3f} | Recall = {ev8['recall']:.3f} | Reduction = {ev8['reduction_factor']:.2f}x | RRF = {ev8['rrf_score']:.3f}")

            # Clean up benchmark tables
            cur.execute("DROP TABLE IF EXISTS rf_bench_molecules CASCADE;")
            cur.execute("DROP TABLE IF EXISTS rf_bench_pockets CASCADE;")

    print("\n" + "=" * 88)
    print("Benchmark Complete!")
    print("=" * 88)


if __name__ == "__main__":
    run_benchmarks()
