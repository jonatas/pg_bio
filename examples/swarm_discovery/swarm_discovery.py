import psycopg
import time
import numpy as np

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

def ensure_swarm_indexes(cur):
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_pending
        ON orphan_discoveries (vector_distance ASC, id ASC)
        WHERE status = 'Pending';
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_orphan_id
        ON orphan_discoveries (orphan_id);
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_orphan_disc_bait_id
        ON orphan_discoveries (bait_id);
    """)


def extract_os_from_header(header_name: str):
    if header_name and "OS=" in header_name and "OX=" in header_name:
        return header_name.split("OS=")[1].split(" OX=")[0].strip()
    return None


def run_swarm(max_steps: int = 0):
    conn = psycopg.connect(DB_URI)
    cur = conn.cursor()

    ensure_swarm_indexes(cur)
    conn.commit()

    print("Starting Optimized Swarm Discovery Agent (FOR UPDATE SKIP LOCKED + Layer 8 Interface Complementarity)...")

    steps = 0
    while True:
        t0 = time.perf_counter()
        # 1. Atomically claim the closest Pending node (lock-free multi-agent concurrency)
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
            RETURNING id, family_name, orphan_id, orphan_organism;
        """)
        row = cur.fetchone()
        conn.commit()

        if not row:
            print("No 'Pending' nodes in the graph! Sleeping...")
            if max_steps > 0:
                break
            time.sleep(10)
            continue

        edge_id, family_name, current_node, cached_parent_org = row
        print(f"\n[Swarm] Exploring node: {current_node} (Family: {family_name})")

        # Verify cognate substrate & catalytic warhead via pg_bio Hückel 4n+2 + WL hash + VF2 (@>) if available
        cur.execute("""
            SELECT substrate_name,
                   wl_hash,
                   molecules_are_identical(ingested_kekule_or_salt_smiles, canonical_query_smiles) AS kekule_arom_eq,
                   (mol @> parse_smiles(warhead_substructure_smiles)) AS warhead_vf2_ok,
                   warhead_substructure_smiles
            FROM swarm_enzyme_substrates
            WHERE LOWER(family_name) = LOWER(%s)
            LIMIT 1;
        """, (family_name,))
        sub_row = cur.fetchone()
        if sub_row:
            sub_name, wl_h, kek_eq, vf2_ok, warhead_s = sub_row
            print(
                f"  [Cheminfo] Cognate Substrate: {sub_name} | WL64={wl_h} "
                f"(Kekulé==Arom: {kek_eq}) | Warhead VF2 (@> {warhead_s}): {vf2_ok}"
            )

        # 2. Get the embedding & header name for this node
        cur.execute("SELECT embedding, name FROM proteins WHERE uniprot_id = %s", (current_node,))
        res = cur.fetchone()
        if not res:
            print(f"Could not find embedding for {current_node}. Marking as Error.")
            cur.execute("UPDATE orphan_discoveries SET status = 'Error' WHERE id = %s", (edge_id,))
            conn.commit()
            continue

        node_emb, parent_header = res
        parent_org = (
            cached_parent_org
            if cached_parent_org and cached_parent_org != "Unknown"
            else (extract_os_from_header(parent_header) or "Unknown")
        )

        # 3. Single-query HNSW search + anti-join duplicate filter + Layer 8 interface_complementarity
        cur.execute("SET LOCAL hnsw.ef_search = 200;")
        cur.execute("""
            WITH raw_neighbors AS (
                SELECT p.uniprot_id,
                       p.name,
                       p.embedding,
                       (p.embedding <=> %s) AS dist
                FROM proteins p
                WHERE p.name ILIKE '%%uncharacterized%%'
                  AND p.uniprot_id != %s
                ORDER BY p.embedding <=> %s ASC
                LIMIT 15
            )
            SELECT rn.uniprot_id,
                   rn.name,
                   rn.dist,
                   interface_complementarity(rn.embedding::vector::real[], %s::vector::real[]) AS icomp
            FROM raw_neighbors rn
            WHERE rn.dist <= 0.15
              AND NOT EXISTS (
                  SELECT 1 FROM orphan_discoveries od
                  WHERE od.orphan_id = rn.uniprot_id
              )
              AND NOT EXISTS (
                  SELECT 1 FROM orphan_discoveries od
                  WHERE od.bait_id = rn.uniprot_id
              )
            ORDER BY rn.dist ASC
            LIMIT 5;
        """, (node_emb, current_node, node_emb, node_emb))

        neighbors = cur.fetchall()
        new_discoveries = 0

        for n_id, n_name, n_dist, icomp in neighbors:
            n_org = extract_os_from_header(n_name)
            if not n_org:
                cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{n_id}",))
                org_res = cur.fetchone()
                n_org = org_res[0] if org_res else "Unknown"

            if parent_org == "Unknown":
                cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{current_node}",))
                parent_org_res = cur.fetchone()
                parent_org = parent_org_res[0] if parent_org_res else "Unknown"

            conf_score = round(0.7 * (1.0 - float(n_dist)) + 0.3 * float(icomp), 4)

            cur.execute("""
                INSERT INTO orphan_discoveries 
                (family_name, bait_id, bait_organism, orphan_id, orphan_organism, vector_distance, status, confidence_score)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            """, (family_name, current_node, parent_org, n_id, n_org, float(n_dist), 'Pending', conf_score))
            new_discoveries += 1

        # Mark as Explored
        cur.execute("UPDATE orphan_discoveries SET status = 'Explored' WHERE id = %s", (edge_id,))
        conn.commit()
        elapsed_ms = (time.perf_counter() - t0) * 1000.0

        print(f"  -> Found {new_discoveries} new branches in {elapsed_ms:.1f} ms. Node {current_node} is now Explored.")
        steps += 1
        if max_steps > 0 and steps >= max_steps:
            break
        time.sleep(0.2)

if __name__ == "__main__":
    run_loop = True
    while run_loop:
        try:
            run_swarm()
        except Exception as e:
            print(f"Swarm crashed: {e}")
            time.sleep(5)
