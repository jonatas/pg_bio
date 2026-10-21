import argparse
import os
import time
import urllib.request
import numpy as np
import psycopg

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"
MODEL_CACHE_DIR = "/Users/jonatas/code/ideia.me/assets/models"

# Seed baits mapped directly to AmazonZyme's most-used organisms & published papers (UFAM / LAB-BIOSIN / PPGBIOTEC)
AMAZONZYME_SEEDS = [
    ("UFAM2 Cocktail (P. citrinum LMI01)", "B1B533", "Penicillium citrinum"),
    ("SSF Glucoamylase (A. awamori GA1/2)", "P69327", "Aspergillus awamori"),
    ("SSF Acid Amylase (A. niger)", "P56271", "Aspergillus niger"),
    ("Industrial Cellulase (T. reesei CBH2)", "P07987", "Trichoderma reesei"),
    ("Liquefaction Amylase (B. licheniformis)", "P06278", "Bacillus licheniformis"),
    ("TPI Clone FosLip2 (C. necator Lipase)", "Q0KD98", "Cupriavidus necator"),
    ("TPI Bioplastics (C. necator PhaC)", "P23608", "Cupriavidus necator"),
    ("TPI Clone Fos4 Ade1 (PDB 8V16)", "A0A8E4UR43", "Uncultured Bacteroidetes (TPI)"),
]


def fetch_or_cache_pdb(uniprot_id: str):
    """Fetch AlphaFold PDB structure (v6 or v4) or crystal PDB (8V16 for A0A8E4UR43)."""
    os.makedirs(MODEL_CACHE_DIR, exist_ok=True)
    if uniprot_id == "A0A8E4UR43" and os.path.exists("/Users/jonatas/code/pg_bio/examples/amazonzyme_poc/8V16.pdb"):
        with open("/Users/jonatas/code/pg_bio/examples/amazonzyme_poc/8V16.pdb", "r") as f:
            return f.read()
    path = os.path.join(MODEL_CACHE_DIR, f"AF-{uniprot_id}-F1-model_v4.pdb")
    if os.path.exists(path):
        with open(path, "r") as f:
            content = f.read()
            if "ATOM" in content:
                return content

    for ver in (6, 4):
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{ver}.pdb"
            pdb_data = urllib.request.urlopen(url, timeout=8).read().decode("utf-8")
            if "ATOM" in pdb_data:
                with open(path, "w") as f:
                    f.write(pdb_data)
                return pdb_data
        except Exception:
            pass
    return None


def evaluate_3d_globular_compactness(uniprot_id: str):
    """
    Computes 3D structural quality metrics from AlphaFold CA coordinates to reject
    coiled-coil rods and disordered spaghetti:
      1. Mean pLDDT confidence across CA atoms
      2. Radius of Gyration (Rg in Angstroms)
      3. Flory Globular Compactness Ratio: Rg / (2.25 * N^0.38)
         - Compact globular enzymes (TIM barrels, alpha/beta hydrolases): ratio ~ 0.85 - 1.28
         - Extended coiled-coils / single alpha-helical rods / disordered chains: ratio > 1.32
    Returns (is_globular, mean_plddt, rg_angstrom, compactness_ratio)
    """
    pdb_str = fetch_or_cache_pdb(uniprot_id)
    if not pdb_str:
        return False, 0.0, 99.9, 9.99

    coords = []
    plddts = []
    for line in pdb_str.splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                b = float(line[60:66])
                coords.append([x, y, z])
                plddts.append(b)
            except ValueError:
                pass

    n = len(coords)
    if n < 60:
        return False, 0.0, 99.9, 9.99

    ca = np.array(coords)
    mean_plddt = float(np.mean(plddts))
    centroid = np.mean(ca, axis=0)
    rg = float(np.sqrt(np.mean(np.sum((ca - centroid) ** 2, axis=1))))
    ideal_globular_rg = 2.25 * (n ** 0.38)
    compactness_ratio = rg / ideal_globular_rg

    # Also check principal axis anisotropy (eigenvalues of gyration tensor) to reject long rods
    centered = ca - centroid
    cov = np.dot(centered.T, centered) / n
    eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
    elongation = float(np.sqrt(eigvals[0] / max(eigvals[1], 1e-4)))

    is_globular = (
        mean_plddt >= 81.5
        and compactness_ratio <= 1.45
        and elongation <= 3.0
    )
    return is_globular, round(mean_plddt, 1), round(rg, 2), round(compactness_ratio, 3)


def init_swarm_schema(cur):
    cur.execute("""
        CREATE TABLE IF NOT EXISTS amazonzyme_swarm_discoveries (
            id SERIAL PRIMARY KEY,
            research_line TEXT NOT NULL,
            bait_id VARCHAR(50) NOT NULL,
            bait_organism TEXT,
            orphan_id VARCHAR(50) NOT NULL UNIQUE,
            orphan_organism TEXT,
            vector_distance REAL NOT NULL,
            hybrid_score REAL,
            mw_kda REAL,
            isoelectric_pi REAL,
            gravy_hydropathy REAL,
            shannon_entropy REAL,
            mean_plddt REAL DEFAULT 0.0,
            rg_angstrom REAL DEFAULT 0.0,
            compactness_ratio REAL DEFAULT 1.0,
            reactor_profile TEXT,
            wetlab_solubility TEXT,
            status VARCHAR(30) DEFAULT 'Pending',
            discovered_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
    """)
    # Ensure new 3D structural columns exist if table was created earlier
    cur.execute("ALTER TABLE amazonzyme_swarm_discoveries ADD COLUMN IF NOT EXISTS mean_plddt REAL DEFAULT 0.0;")
    cur.execute("ALTER TABLE amazonzyme_swarm_discoveries ADD COLUMN IF NOT EXISTS rg_angstrom REAL DEFAULT 0.0;")
    cur.execute("ALTER TABLE amazonzyme_swarm_discoveries ADD COLUMN IF NOT EXISTS compactness_ratio REAL DEFAULT 1.0;")

    # Seed initial nodes if table is empty
    cur.execute("SELECT count(*) FROM amazonzyme_swarm_discoveries;")
    if cur.fetchone()[0] == 0:
        print("[AmazonZyme Swarm] Seeding initial domain baits from bio_demo...")
        for line_name, preferred_id, fallback_org in AMAZONZYME_SEEDS:
            cur.execute("SELECT uniprot_id, sequence FROM proteins WHERE uniprot_id = %s", (preferred_id,))
            row = cur.fetchone()
            if not row:
                kw = line_name.split(" (")[0]
                cur.execute(
                    "SELECT uniprot_id, sequence FROM proteins WHERE name ILIKE %s AND length(sequence) BETWEEN 250 AND 550 OFFSET 1 LIMIT 1",
                    (f"%{kw}%",)
                )
                row = cur.fetchone()
            if row:
                b_id, b_seq = row
                _, plddt, rg, comp_ratio = evaluate_3d_globular_compactness(b_id)
                cur.execute("""
                    SELECT 
                        round((molecular_weight(%s)/1000.0)::numeric, 2),
                        round(isoelectric_point(%s)::numeric, 2),
                        round(gravy_score(%s)::numeric, 3),
                        round(sequence_entropy(%s)::numeric, 3)
                """, (b_seq, b_seq, b_seq, b_seq))
                mw, pi, gravy, ent = cur.fetchone()
                reactor = "Acidic_SSF (pH 3.5-5.5)" if pi < 5.8 else "Neutral_Alkaline (pH 6.0-8.5)"
                sol = "SOLUBLE_GLOBULAR_BARREL" if gravy < 0.0 else "GLOBULAR_CHAPERONE_ASSIST"
                cur.execute("""
                    INSERT INTO amazonzyme_swarm_discoveries
                    (research_line, bait_id, bait_organism, orphan_id, orphan_organism,
                     vector_distance, hybrid_score, mw_kda, isoelectric_pi, gravy_hydropathy,
                     shannon_entropy, mean_plddt, rg_angstrom, compactness_ratio,
                     reactor_profile, wetlab_solubility, status)
                    VALUES (%s, %s, %s, %s, %s, 0.0, 0.0, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'Pending')
                    ON CONFLICT (orphan_id) DO NOTHING;
                """, (
                    line_name, b_id, fallback_org, b_id, fallback_org,
                    mw, pi, gravy, ent, plddt, rg, comp_ratio, reactor, sol
                ))


def run_swarm_step(conn, cur):
    cur.execute("""
        SELECT id, research_line, orphan_id, orphan_organism
        FROM amazonzyme_swarm_discoveries
        WHERE status = 'Pending'
        ORDER BY vector_distance ASC, id ASC
        LIMIT 1;
    """)
    row = cur.fetchone()
    if not row:
        return False

    node_id, research_line, current_id, current_org = row
    print(f"\n[AmazonZyme Swarm] Exploring {current_id} ({current_org}) | Line: {research_line}")

    cur.execute("UPDATE amazonzyme_swarm_discoveries SET status = 'Exploring' WHERE id = %s", (node_id,))
    conn.commit()

    # Fetch parent embedding and sequence first so Postgres uses the HNSW index with ef_search=300
    cur.execute("SELECT embedding, sequence, length(sequence) FROM proteins WHERE uniprot_id = %s", (current_id,))
    prow = cur.fetchone()
    if not prow:
        return False
    parent_emb, parent_seq, parent_len = prow

    cur.execute("SET LOCAL hnsw.ef_search = 400;")
    cur.execute("""
        WITH nearest_candidates AS (
            SELECT 
                p.uniprot_id,
                p.name,
                p.sequence,
                (p.embedding <=> %s) AS v_dist
            FROM proteins p
            WHERE p.uniprot_id != %s
              AND p.name ILIKE '%%[AmazonZyme-Organism%%'
              AND length(p.sequence) BETWEEN GREATEST(180, (%s * 0.65)::int) AND (%s * 1.35)::int
            ORDER BY p.embedding <=> %s ASC
            LIMIT 25
        )
        SELECT 
            nc.uniprot_id,
            nc.name,
            nc.v_dist,
            (0.65 * nc.v_dist) + (0.35 * (150.0 / (150.0 + sequence_alignment_score(nc.sequence, %s)::float8))) AS h_score,
            round((molecular_weight(nc.sequence) / 1000.0)::numeric, 2) AS mw_kda,
            round(isoelectric_point(nc.sequence)::numeric, 2) AS pi,
            round(gravy_score(nc.sequence)::numeric, 3) AS gravy,
            round(sequence_entropy(nc.sequence)::numeric, 3) AS entropy
        FROM nearest_candidates nc
        WHERE nc.v_dist <= 0.72
          AND sequence_entropy(nc.sequence) >= 3.85
          AND gravy_score(nc.sequence) BETWEEN -0.55 AND 0.35
        ORDER BY h_score ASC
        LIMIT 12;
    """, (parent_emb, current_id, parent_len, parent_len, parent_emb, parent_seq))

    candidates = cur.fetchall()
    added = 0

    for cand_id, cand_name, v_dist, h_score, mw_kda, pi, gravy, entropy in candidates:
        if added >= 3:
            break

        # Gate 3D structure: reject coiled-coils, single-helix rods, and low-pLDDT disordered structures
        is_globular, plddt, rg, comp_ratio = evaluate_3d_globular_compactness(cand_id)
        if not is_globular:
            print(
                f"  - Rejected coil/disordered {cand_id}: pLDDT={plddt}, Rg={rg}A, compactness={comp_ratio} (>1.28 or elongated rod)"
            )
            continue

        org_name = "Unknown Amazonian/Extremophile Taxon"
        if "OS=" in cand_name and "OX=" in cand_name:
            org_name = cand_name.split("OS=")[1].split(" OX=")[0].strip()
        else:
            try:
                cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{cand_id}",))
                r = cur.fetchone()
                if r and r[0]:
                    org_name = r[0]
            except Exception:
                pass

        reactor_profile = (
            "Acidic_SSF (pH 3.5-5.5)" if float(pi) < 5.8 else "Neutral_Alkaline (pH 6.0-8.5)"
        )
        wetlab_solubility = (
            f"SOLUBLE_GLOBULAR (Rg={rg}A, Ratio={comp_ratio})"
            if float(gravy) < 0.0
            else f"GLOBULAR_HYDROPHOBIC (Rg={rg}A, Ratio={comp_ratio})"
        )

        cur.execute("""
            INSERT INTO amazonzyme_swarm_discoveries
            (research_line, bait_id, bait_organism, orphan_id, orphan_organism,
             vector_distance, hybrid_score, mw_kda, isoelectric_pi, gravy_hydropathy,
             shannon_entropy, mean_plddt, rg_angstrom, compactness_ratio,
             reactor_profile, wetlab_solubility, status)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'Pending')
            ON CONFLICT (orphan_id) DO NOTHING
            RETURNING id;
        """, (
            research_line, current_id, current_org, cand_id, org_name,
            float(v_dist), float(h_score), float(mw_kda), float(pi), float(gravy),
            float(entropy), plddt, rg, comp_ratio, reactor_profile, wetlab_solubility
        ))
        if cur.fetchone():
            added += 1
            print(
                f"  + Verified 3D Globular {cand_id} ({org_name}) | dist={v_dist:.4f} | pLDDT={plddt} "
                f"| Rg={rg}A (ratio={comp_ratio}) | MW={mw_kda}kDa | pI={pi}"
            )

    cur.execute("UPDATE amazonzyme_swarm_discoveries SET status = 'Explored' WHERE id = %s", (node_id,))
    conn.commit()
    print(f"  -> Node {current_id} explored ({added} verified globular enzymes queued).")
    return True


def main():
    parser = argparse.ArgumentParser(description="AmazonZyme Globular Enzyme Discovery Swarm on pg_bio")
    parser.add_argument("--cycles", type=int, default=0, help="Number of swarm nodes to explore (0 = infinite daemon)")
    parser.add_argument("--reset", action="store_true", help="Reset table and re-seed with strict 3D globular gate")
    args = parser.parse_args()

    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            if args.reset:
                cur.execute("DROP TABLE IF EXISTS amazonzyme_swarm_discoveries;")
                conn.commit()
            init_swarm_schema(cur)
            conn.commit()

            step = 0
            while True:
                ok = run_swarm_step(conn, cur)
                if not ok:
                    print("No Pending nodes left in amazonzyme_swarm_discoveries.")
                    break
                step += 1
                if args.cycles > 0 and step >= args.cycles:
                    break
                time.sleep(0.3)

            print("\n=== TOP AMAZONZYME GLOBULAR SWARM DISCOVERIES (COILS REJECTED) ===")
            cur.execute("""
                SELECT 
                    research_line,
                    bait_id,
                    orphan_id,
                    left(orphan_organism, 28) AS organism,
                    round(vector_distance::numeric, 4) AS v_dist,
                    round(hybrid_score::numeric, 4) AS h_score,
                    mw_kda,
                    isoelectric_pi AS pi,
                    gravy_hydropathy AS gravy,
                    mean_plddt,
                    rg_angstrom,
                    compactness_ratio
                FROM amazonzyme_swarm_discoveries
                WHERE vector_distance > 0.0
                ORDER BY hybrid_score ASC
                LIMIT 20;
            """)
            rows = cur.fetchall()
            for r in rows:
                print(
                    f"[{r[0]:28s}] {r[1]} -> {r[2]} ({r[3]:28s}) | "
                    f"dist={r[4]} | hyb={r[5]} | {r[6]:5.1f}kDa | pI={r[7]:4.2f} | GRAVY={r[8]:6.3f} | "
                    f"pLDDT={r[9]:4.1f} | Rg={r[10]:5.2f}A (comp={r[11]:4.3f})"
                )

if __name__ == "__main__":
    main()
