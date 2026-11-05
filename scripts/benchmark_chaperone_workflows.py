import math
import time
import psycopg

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"


def build_synthetic_complex_pdb(
    n_binder: int = 24,
    n_target: int = 12,
    plddt: float = 92.5,
    coord_jitter: float = 0.0,
    expose_hydrophobics: bool = False,
    bury_unsatisfied_polar: bool = False,
) -> str:
    """
    Generates a realistic 3D all-atom backbone + polar/apolar sidechain PDB complex
    (Chain A: de novo binder from RFdiffusion3/MPNN; Chain B: target protein/epitope).
    """
    lines = []
    serial = 1
    for i in range(n_binder):
        theta = i * 0.58
        phi = 0.45 + (i / max(1, n_binder)) * 2.1
        r = 7.2 + (0.25 if i % 3 == 0 else -0.15)
        cx = r * math.cos(theta) * math.sin(phi) + coord_jitter * math.sin(i * 1.7)
        cy = r * math.sin(theta) * math.sin(phi) + coord_jitter * math.cos(i * 1.3)
        cz = r * math.cos(phi) + coord_jitter * math.sin(i * 0.9)

        if expose_hydrophobics:
            res_name = "PHE" if i % 2 == 0 else "LEU"
        else:
            res_name = ("SER", "ALA", "GLU", "LYS", "ASN", "GLY")[i % 6]

        for aname, dx, dy, dz, elem in (
            ("N", -0.75, 0.0, 0.0, "N"),
            ("CA", 0.0, 0.0, 0.0, "C"),
            ("C", 0.78, 0.0, 0.0, "C"),
            ("O", 1.35, 0.85, 0.0, "O"),
        ):
            lines.append(
                f"ATOM  {serial:5d}  {aname:<3s} {res_name:3s} A{i + 1:4d}    "
                f"{cx + dx:8.3f}{cy + dy:8.3f}{cz + dz:8.3f}  1.00{plddt:6.2f}           {elem:>2s}  "
            )
            serial += 1

    # Add target Chain B across the binding interface (2.85 A H-bond register with binder)
    for j in range(n_target):
        theta = j * 0.58
        phi = 0.45 + (j / max(1, n_binder)) * 2.1
        r = 7.2 + (0.25 if j % 3 == 0 else -0.15)
        ox = r * math.cos(theta) * math.sin(phi) + 1.35
        oy = r * math.sin(theta) * math.sin(phi) + 0.85
        oz = r * math.cos(phi)

        res_name = "ASN" if j % 2 == 0 else "THR"
        # Place N at 2.85 A from binder O to form cross-interface H-bonds
        nx = ox + 2.85
        ny = oy
        nz = oz
        for aname, x, y, z, elem in (
            ("N", nx, ny, nz, "N"),
            ("CA", nx + 1.40, ny, nz, "C"),
            ("C", nx + 2.20, ny + 0.80, nz, "C"),
            ("O", nx + 2.90, ny + 1.55, nz, "O"),
        ):
            lines.append(
                f"ATOM  {serial:5d}  {aname:<3s} {res_name:3s} B{101 + j:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{plddt:6.2f}           {elem:>2s}  "
            )
            serial += 1

    if bury_unsatisfied_polar:
        # Place isolated buried polar N/O atoms encapsulated by carbon shell with no H-bond partner
        for k in range(5):
            bx, by, bz = 0.0, 0.0, k * 0.1
            lines.append(
                f"ATOM  {serial:5d}  ND2 ASN A{80 + k:4d}    "
                f"{bx:8.3f}{by:8.3f}{bz:8.3f}  1.00{plddt:6.2f}           N  "
            )
            serial += 1
            for s_idx, (sx, sy, sz) in enumerate(
                [(1.8, 0, 0), (-1.8, 0, 0), (0, 1.8, 0), (0, -1.8, 0), (0, 0, 1.8), (0, 0, -1.8)]
            ):
                lines.append(
                    f"ATOM  {serial:5d}  CB  ALA A{90 + s_idx:4d}    "
                    f"{bx + sx:8.3f}{by + sy:8.3f}{bz + sz:8.3f}  1.00{plddt:6.2f}           C  "
                )
                serial += 1

    return "\n".join(lines) + "\n"


def python_evaluate_chaperone_row(cand_pdb: str, ref_pdb: str) -> dict:
    """
    Reference pure-Python evaluation of Chaperone's post-oracle block
    (parsing PDB atoms, 64-pt Fibonacci sphere Shrake-Rupley SASA/RASA, H-bond scan, and CA RMSD).
    """
    atoms = []
    ca_cand = []
    for line in cand_pdb.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            aname = line[12:16].strip()
            rname = line[17:20].strip()
            chain = line[21:22].strip() or "A"
            rseq = int(line[22:26].strip() or 0)
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            bfac = float(line[60:66]) if len(line) >= 66 else 0.0
            elem = aname[0]
            vdw = 1.70 if elem == "C" else (1.55 if elem == "N" else 1.52)
            atoms.append((aname, rname, chain, rseq, x, y, z, vdw, elem in ("N", "O")))
            if aname == "CA":
                ca_cand.append((x, y, z, bfac))

    ca_ref = []
    for line in ref_pdb.splitlines():
        if line.startswith("ATOM  ") and len(line) >= 54 and line[12:16].strip() == "CA":
            ca_ref.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

    # 1. RMSD
    n_ca = min(len(ca_cand), len(ca_ref))
    rmsd = 0.0
    if n_ca > 0:
        rmsd = math.sqrt(
            sum(
                (ca_cand[i][0] - ca_ref[i][0]) ** 2
                + (ca_cand[i][1] - ca_ref[i][1]) ** 2
                + (ca_cand[i][2] - ca_ref[i][2]) ** 2
                for i in range(n_ca)
            )
            / n_ca
        )

    # 2. 64-point Fibonacci sphere SASA
    golden = math.pi * (3.0 - math.sqrt(5.0))
    pts = []
    for k in range(64):
        yk = 1.0 - (2.0 * k + 1.0) / 64.0
        rk = math.sqrt(max(0.0, 1.0 - yk * yk))
        th = golden * k
        pts.append((rk * math.cos(th), yk, rk * math.sin(th)))

    total_sasa = 0.0
    for i, a in enumerate(atoms):
        ri = a[7] + 1.4
        xi, yi, zi = a[4], a[5], a[6]
        nbrs = []
        for j, b in enumerate(atoms):
            if i == j:
                continue
            rj = b[7] + 1.4
            dx, dy, dz = xi - b[4], yi - b[5], zi - b[6]
            if dx * dx + dy * dy + dz * dz < (ri + rj) ** 2:
                nbrs.append((b[4], b[5], b[6], rj * rj))
        acc = 0
        for ux, uy, uz in pts:
            px, py, pz = xi + ri * ux, yi + ri * uy, zi + ri * uz
            buried = False
            for bx, by, bz, rj2 in nbrs:
                if (px - bx) ** 2 + (py - by) ** 2 + (pz - bz) ** 2 < rj2:
                    buried = True
                    break
            if not buried:
                acc += 1
        total_sasa += 4.0 * math.pi * ri * ri * (acc / 64.0)

    # 3. H-bond scan
    hbonds = 0
    polars = [a for a in atoms if a[8]]
    for i in range(len(polars)):
        for j in range(i + 1, len(polars)):
            if polars[i][3] == polars[j][3] and polars[i][2] == polars[j][2]:
                continue
            d2 = (
                (polars[i][4] - polars[j][4]) ** 2
                + (polars[i][5] - polars[j][5]) ** 2
                + (polars[i][6] - polars[j][6]) ** 2
            )
            if 2.2**2 <= d2 <= 3.5**2:
                hbonds += 1

    return {"rmsd": rmsd, "total_sasa": total_sasa, "hbonds": hbonds}


def main():
    print("=" * 88)
    print("pg_bio Benchmark: RosettaCommons Chaperone Workflows (RFD3 -> MPNN -> AF3/Boltz-2)")
    print("=" * 88)

    ref_pdb = build_synthetic_complex_pdb(n_binder=24, n_target=12, plddt=94.0, coord_jitter=0.0)
    good_cand = build_synthetic_complex_pdb(n_binder=24, n_target=12, plddt=93.2, coord_jitter=0.35)
    drifted_cand = build_synthetic_complex_pdb(n_binder=24, n_target=12, plddt=62.0, coord_jitter=4.20)
    hydro_cand = build_synthetic_complex_pdb(
        n_binder=24, n_target=12, plddt=89.0, coord_jitter=0.45, expose_hydrophobics=True
    )
    buns_cand = build_synthetic_complex_pdb(
        n_binder=24, n_target=12, plddt=88.5, coord_jitter=0.40, bury_unsatisfied_polar=True
    )

    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute("DROP TABLE IF EXISTS chaperone_benchmark_designs;")
            cur.execute(
                """
                CREATE TABLE chaperone_benchmark_designs (
                    design_id INT PRIMARY KEY,
                    archetype TEXT,
                    parent_rfd3_id TEXT,
                    parent_mpnn_id TEXT,
                    cycle_num INT,
                    sequence TEXT,
                    rfd3_backbone_pdb TEXT,
                    oracle_pdb TEXT
                );
                """
            )

            archetypes = [
                ("High-Affinity Satisfied Binder (Pass)", good_cand, "MSEKANGSEKANGSEKANGSEKAN"),
                ("High scRMSD Backbone Drift (Fail)", drifted_cand, "MSEKANGSEKANGSEKANGSEKAN"),
                ("Surface Hydrophobic Patch Aggregation (Fail)", hydro_cand, "FLFLFLFLFLFLFLFLFLFLFLFL"),
                ("Buried Unsatisfied Polar BUNS Trap (Fail)", buns_cand, "MSEKANGSEKANGSEKANGNNNNN"),
            ]

            n_rows = 400
            rows = []
            for idx in range(1, n_rows + 1):
                label, cand_pdb, seq = archetypes[(idx - 1) % len(archetypes)]
                cycle = idx % 3
                rows.append(
                    (
                        idx,
                        label,
                        f"chp_rfd3_c0_{idx:04d}",
                        f"chp_mpnn_c{cycle}_{idx:04d}",
                        cycle,
                        seq,
                        ref_pdb,
                        cand_pdb,
                    )
                )

            cur.executemany(
                """
                INSERT INTO chaperone_benchmark_designs
                (design_id, archetype, parent_rfd3_id, parent_mpnn_id, cycle_num, sequence, rfd3_backbone_pdb, oracle_pdb)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
                """,
                rows,
            )
            conn.commit()

            # 1. Verify Archetype Discrimination on the 4 canonical design states
            cur.execute(
                """
                SELECT
                    d.archetype,
                    ROUND(eval.sc_rmsd::numeric, 3) AS sc_rmsd,
                    ROUND(eval.mean_plddt::numeric, 1) AS plddt,
                    ROUND(eval.interface_delta_sasa::numeric, 1) AS delta_sasa,
                    eval.cross_interface_hbonds AS iface_hb,
                    eval.buried_unsatisfied_polars AS buns,
                    ROUND(eval.surface_hydrophobic_fraction::numeric, 3) AS surf_hydro,
                    ROUND(eval.design_composite_score::numeric, 3) AS score,
                    eval.passes_chaperone_gate AS pass_gate
                FROM chaperone_benchmark_designs d
                CROSS JOIN LATERAL evaluate_chaperone_block(d.oracle_pdb, d.rfd3_backbone_pdb, 'A', 'B') AS eval
                WHERE d.design_id <= 4
                ORDER BY d.design_id;
                """
            )
            print("\n[1] Chaperone Post-Oracle Archetype Audit (RFdiffusion3 -> Foundry-MPNN -> Boltz-2/AF3):")
            for r in cur.fetchall():
                print(
                    f"  - {r[0]:<44s} | scRMSD={r[1]} A | pLDDT={r[2]} | dSASA={r[3]} A^2 | "
                    f"iHB={r[4]} | BUNS={r[5]} | SurfHydro={r[6]} | Score={r[7]} | Gate={r[8]}"
                )

            # 2. Benchmark Default Python Pull + Evaluation vs. In-Database pg_bio Evaluation
            t0 = time.perf_counter()
            cur.execute(
                "SELECT design_id, oracle_pdb, rfd3_backbone_pdb FROM chaperone_benchmark_designs;"
            )
            fetched = cur.fetchall()
            py_survivors = []
            for did, c_pdb, r_pdb in fetched:
                res = python_evaluate_chaperone_row(c_pdb, r_pdb)
                if res["rmsd"] <= 2.0 and res["hbonds"] >= 5:
                    py_survivors.append(did)
            py_ms = (time.perf_counter() - t0) * 1000.0

            # Warmup & time native pg_bio SQL evaluation + lineage hash
            cur.execute(
                """
                SELECT
                    d.design_id,
                    chaperone_lineage_hash('boltz2', d.parent_mpnn_id, d.sequence, d.cycle_num) AS lineage_id,
                    eval.sc_rmsd,
                    eval.interface_delta_sasa,
                    eval.cross_interface_hbonds,
                    eval.buried_unsatisfied_polars,
                    eval.design_composite_score
                FROM chaperone_benchmark_designs d
                CROSS JOIN LATERAL evaluate_chaperone_block(d.oracle_pdb, d.rfd3_backbone_pdb, 'A', 'B') AS eval
                WHERE eval.passes_chaperone_gate = TRUE;
                """
            )
            _ = cur.fetchall()

            t1 = time.perf_counter()
            cur.execute(
                """
                SELECT
                    d.design_id,
                    chaperone_lineage_hash('boltz2', d.parent_mpnn_id, d.sequence, d.cycle_num) AS lineage_id,
                    eval.sc_rmsd,
                    eval.interface_delta_sasa,
                    eval.cross_interface_hbonds,
                    eval.buried_unsatisfied_polars,
                    eval.design_composite_score
                FROM chaperone_benchmark_designs d
                CROSS JOIN LATERAL evaluate_chaperone_block(d.oracle_pdb, d.rfd3_backbone_pdb, 'A', 'B') AS eval
                WHERE eval.passes_chaperone_gate = TRUE;
                """
            )
            pg_survivors = cur.fetchall()
            pg_ms = (time.perf_counter() - t1) * 1000.0

            speedup = py_ms / max(pg_ms, 0.001)
            print(f"\n[2] Throughput Benchmark ({n_rows} Multi-Chain Complexes, 144 Heavy Atoms Each):")
            print(f"  - Python Pull + Shrake-Rupley SASA + H-Bond + RMSD : {py_ms:8.2f} ms")
            print(f"  - Native pg_bio Fused Rust SQL Evaluation + Lineage : {pg_ms:8.2f} ms ({speedup:.1f}x faster)")
            print(f"  - Surviving Golden-Gate Ready Designs              : {len(pg_survivors)} / {n_rows}")

            # 3. Audit all 8 Official RosettaCommons Chaperone Core Examples
            t2 = time.perf_counter()
            cur.execute(
                """
                SELECT
                    ex.example_name,
                    ex.category,
                    ROUND(audit.sc_rmsd::numeric, 3) AS sc_rmsd,
                    ROUND(audit.mean_plddt::numeric, 1) AS plddt,
                    ROUND(audit.interface_delta_sasa::numeric, 1) AS delta_sasa,
                    audit.cross_interface_hbonds AS iface_hb,
                    audit.buried_unsatisfied_polars AS buns,
                    ROUND(audit.active_site_or_motif_rmsd::numeric, 3) AS motif_rmsd,
                    audit.ligand_or_cofactor_verified AS verified,
                    ROUND(audit.composite_score::numeric, 3) AS score,
                    audit.passes_example_gate AS pass_gate
                FROM chaperone_core_examples() AS ex
                CROSS JOIN LATERAL evaluate_chaperone_core_example_block(
                    ex.example_name, ex.candidate_pdb, ex.reference_pdb, ex.ligand_smiles
                ) AS audit;
                """
            )
            core_rows = cur.fetchall()
            core_ms = (time.perf_counter() - t2) * 1000.0
            print(f"\n[3] Official RosettaCommons Chaperone Core Examples Audit (8 Workflows in {core_ms:.2f} ms):")
            for r in core_rows:
                print(
                    f"  - {r[0]:<20s} ({r[1]:<46s}) | scRMSD={r[2]} A | pLDDT={r[3]} | "
                    f"dSASA={r[4]} A^2 | iHB={r[5]} | BUNS={r[6]} | MotifRMSD={r[7]} | "
                    f"Verified={r[8]} | Score={r[9]} | Gate={r[10]}"
                )


if __name__ == "__main__":
    main()

