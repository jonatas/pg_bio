"""
End-to-end runnable example for pg_bio Layer 15: RangeFinder (arXiv:2611.03379v1)
Pocket-Interval GiST & Composite B-Tree Virtual Screening.

Run with:
    python examples/rangefinder_pocket_screening/run_rangefinder_example.py
"""

import os
import sys
import psycopg

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, os.path.join(REPO_ROOT, "pgbio-py"))

from pgbio.client import PgBioClient  # noqa: E402

DB_URI = os.environ.get("PGBIO_DB_URI", "postgresql://jonatas@localhost:28818/bio_demo")
SQL_FILE = os.path.join(os.path.dirname(__file__), "schema_and_queries.sql")


def main():
    print("=" * 84)
    print("pg_bio Layer 15 — RangeFinder (arXiv:2611.03379v1) Pocket-Interval GiST Example")
    print("=" * 84)

    client = PgBioClient(DB_URI)

    # 1. Compute 9 RangeFinder physicochemical properties via PgBioClient
    ibuprofen_smiles = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
    props = client.compute_ligand_properties(ibuprofen_smiles)
    print(f"\n[1] Ligand Properties for Ibuprofen ({ibuprofen_smiles}):")
    print(
        f"    MW={props.mw:.2f} Da | logP={props.logp:.2f} | HBD={props.hbd:.0f} | "
        f"HBA={props.hba:.0f} | TPSA={props.tpsa:.1f} A^2 | Rot={props.rotatable_bonds:.0f} | "
        f"Charge={props.formal_charge:.0f} | Arom={props.aromatic_rings:.0f} | Fsp3={props.fsp3:.3f}"
    )

    # 2. Evaluate RangeFinder Screening (8-property non-degenerate vs 9-property)
    binders = [
        "CC(=O)Oc1ccccc1C(=O)O",          # Aspirin
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1",     # Ibuprofen
        "COc1ccc2cc(C(C)C(=O)O)ccc2c1",   # Naproxen
        "O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl", # Diclofenac
    ]
    background = binders + [
        "C1CCCCC1",                       # Cyclohexane
        "CCCCCCCCCCCCCCCCCC(=O)O",        # Stearic Acid
        "NCC(=O)O",                       # Glycine
        "CC1=CCC(C(=CCC=C(C)C)C)CC1",     # Bisabolene
    ]
    eval_8 = client.evaluate_rangefinder_screening(
        "COX2_NSAID_POCKET", binders, background, include_formal_charge=False
    )
    print("\n[2] RangeFinder 8-Property Screening Audit (COX2_NSAID_POCKET):")
    print(
        f"    Range-IoU={eval_8.range_iou:.3f} | Recall={eval_8.recall:.3f} | "
        f"Reduction={eval_8.reduction_factor:.2f}x | RRF={eval_8.rrf_score:.3f}"
    )

    # 3. Execute schema_and_queries.sql in PostgreSQL 18
    with open(SQL_FILE, "r") as f:
        sql_text = f.read()

    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            cur.execute(sql_text)

            print("\n[3] Forward Virtual Screen (COX2_NSAID_POCKET -> rf_example_ligands via Composite B-Tree):")
            cur.execute(
                """
                SELECT l.compound_name, ROUND(l.mw::numeric, 2), ROUND(l.logp::numeric, 2), ROUND(l.tpsa::numeric, 1)
                FROM pocket_interval_bounds(
                    (SELECT rf_profile FROM rf_example_pockets WHERE pocket_id = 'COX2_NSAID_POCKET')
                ) b
                JOIN rf_example_ligands l
                  ON l.mw              BETWEEN b.mw_min   AND b.mw_max
                 AND l.logp            BETWEEN b.logp_min AND b.logp_max
                 AND l.tpsa            BETWEEN b.tpsa_min AND b.tpsa_max
                 AND l.hba             BETWEEN b.hba_min  AND b.hba_max
                 AND l.hbd             BETWEEN b.hbd_min  AND b.hbd_max
                 AND l.rotatable_bonds BETWEEN b.rot_min  AND b.rot_max
                 AND l.aromatic_rings  BETWEEN b.arom_min AND b.arom_max
                 AND l.fsp3            BETWEEN b.fsp3_min AND b.fsp3_max
                ORDER BY l.mw;
                """
            )
            for name, mw, logp, tpsa in cur.fetchall():
                print(f"    - {name:<15s} | MW: {mw} Da | logP: {logp} | TPSA: {tpsa} A^2")

            print("\n[4] Reverse Target Fishing (Vanillin 'COc1cc(C=O)ccc1O' -> rf_example_pockets via Multi-Column GiST `@>`):")
            cur.execute(
                """
                WITH lig AS (
                    SELECT * FROM compute_ligand_properties_row('COc1cc(C=O)ccc1O')
                )
                SELECT p.pocket_id, p.uniprot_id, p.target_name, ROUND(p.aromatic_probability::numeric, 3)
                FROM rf_example_pockets p, lig
                WHERE p.mw_range   @> lig.mw::numeric
                  AND p.logp_range @> lig.logp::numeric
                  AND p.tpsa_range @> lig.tpsa::numeric
                  AND p.hba_range  @> lig.hba::numeric
                  AND p.hbd_range  @> lig.hbd::numeric
                  AND p.rot_range  @> lig.rotatable_bonds::numeric
                  AND p.arom_range @> lig.aromatic_rings::numeric
                  AND p.fsp3_range @> lig.fsp3::numeric;
                """
            )
            for pid, uniprot, target_name, arom_p in cur.fetchall():
                print(f"    - {pid} ({uniprot}) — {target_name} [Arom Prob: {arom_p}]")

            # Clean up example tables
            cur.execute("DROP TABLE IF EXISTS rf_example_pockets CASCADE;")
            cur.execute("DROP TABLE IF EXISTS rf_example_ligands CASCADE;")

    print("\n✅ RangeFinder example completed successfully!")


if __name__ == "__main__":
    main()
