import psycopg
import time
import numpy as np

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

def run_swarm():
    conn = psycopg.connect(DB_URI)
    cur = conn.cursor()
    
    print("Starting Swarm Discovery Agent...")
    
    while True:
        # 1. Pick a Pending node
        cur.execute("""
            SELECT id, family_name, orphan_id 
            FROM orphan_discoveries 
            WHERE status = 'Pending' 
            ORDER BY vector_distance ASC 
            LIMIT 1;
        """)
        row = cur.fetchone()
        
        if not row:
            print("No 'Pending' nodes in the graph! Sleeping...")
            time.sleep(10)
            continue
            
        edge_id, family_name, current_node = row
        print(f"\n[Swarm] Exploring node: {current_node} (Family: {family_name})")
        
        # Mark as Exploring so we don't pick it again concurrently
        cur.execute("UPDATE orphan_discoveries SET status = 'Exploring' WHERE id = %s", (edge_id,))
        conn.commit()
        
        # 2. Get the embedding for this node
        cur.execute("SELECT embedding FROM proteins WHERE uniprot_id = %s", (current_node,))
        res = cur.fetchone()
        if not res:
            print(f"Could not find embedding for {current_node}. Marking as Error.")
            cur.execute("UPDATE orphan_discoveries SET status = 'Error' WHERE id = %s", (edge_id,))
            conn.commit()
            continue
            
        node_emb = res[0]
        
        # 3. Find top 5 closest uncharacterized neighbors
        cur.execute("""
            SELECT uniprot_id, (embedding <=> %s) as dist
            FROM proteins
            WHERE name ILIKE '%%uncharacterized%%'
              AND uniprot_id != %s
            ORDER BY embedding <=> %s ASC
            LIMIT 5;
        """, (node_emb, current_node, node_emb))
        
        neighbors = cur.fetchall()
        new_discoveries = 0
        
        for n_id, n_dist in neighbors:
            # Check if this neighbor is already in the graph
            cur.execute("""
                SELECT 1 FROM orphan_discoveries 
                WHERE orphan_id = %s OR bait_id = %s
            """, (n_id, n_id))
            if cur.fetchone():
                continue # Already in graph
                
            # If distance is too far, ignore
            if n_dist > 0.15:
                continue
                
            # Enrich organism
            cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{n_id}",))
            org_res = cur.fetchone()
            n_org = org_res[0] if org_res else "Unknown"
            
            # Get parent organism
            cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{current_node}",))
            parent_org_res = cur.fetchone()
            parent_org = parent_org_res[0] if parent_org_res else "Unknown"
            
            # Insert new edge
            cur.execute("""
                INSERT INTO orphan_discoveries 
                (family_name, bait_id, bait_organism, orphan_id, orphan_organism, vector_distance, status, confidence_score)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            """, (family_name, current_node, parent_org, n_id, n_org, n_dist, 'Pending', 1.0 - n_dist))
            new_discoveries += 1
            
        # Mark as Explored
        cur.execute("UPDATE orphan_discoveries SET status = 'Explored' WHERE id = %s", (edge_id,))
        conn.commit()
        
        print(f"  -> Found {new_discoveries} new branches. Node {current_node} is now Explored.")
        time.sleep(1) # Slight pause

if __name__ == "__main__":
    run_loop = True
    while run_loop:
        try:
            run_swarm()
        except Exception as e:
            print(f"Swarm crashed: {e}")
            time.sleep(5)
