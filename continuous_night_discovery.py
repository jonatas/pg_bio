import psycopg
import time
import urllib.request
import os
from datetime import datetime

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

FAMILIES = [
    # CRISPR & Biotechnology
    "Cas9", "Argonaute", "Taq polymerase", "Luciferase", "Topoisomerase", "Helicase", "Ligase", "Nuclease",
    
    # Industrial & Bioremediation
    "Amylase", "Lipase", "Protease", "Xylanase", "Chitinase", "Keratinase", "Pectinase", "Cellulase", "Cutinase",
    "Dehalogenase", "Amidase", "Methane monooxygenase", "Laccase", "Alkane hydroxylase", "PETase",
    
    # Extreme Environment & Antioxidants
    "Catalase", "Peroxidase", "Superoxide dismutase", "Manganese peroxidase", "Metallothionein",
    
    # Carbon / Nitrogen Cycle
    "Nitrogenase", "Nitrate reductase", "Nitrite reductase", "Carbonic anhydrase", "RuBisCO", "Hydrogenase", "PHA synthase",
    
    # General Enzymes
    "Kinase", "Phosphatase", "Oxidoreductase", "Desaturase", "Oxygenase", "Hydrolase", "Transferase", "Isomerase"
]

def fetch_pdb_or_fold(cur, uniprot_id, seq, output_path):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except:
        print(f"AlphaFold missing {uniprot_id}. Triggering ESMFold API via pg_bio...")
        try:
            cur.execute("SELECT bio_fold_sequence(%s)", (seq,))
            pdb_data = cur.fetchone()[0]
            if pdb_data and "ATOM" in pdb_data:
                with open(output_path, "w") as f:
                    f.write(pdb_data)
                return True
        except Exception as e:
            print(f"ESMFold failed: {e}")
        return False


def get_uniprot_details(uniprot_id):
    import json, urllib.request
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read())
            
            function = "No specific function described."
            name = "Unknown protein"
            
            if 'proteinDescription' in data and 'recommendedName' in data['proteinDescription']:
                name = data['proteinDescription']['recommendedName']['fullName']['value']
            elif 'proteinDescription' in data and 'submittedName' in data['proteinDescription']:
                name = data['proteinDescription']['submittedName'][0]['fullName']['value']
            
            for comment in data.get('comments', []):
                if comment['commentType'] == 'FUNCTION':
                    function = comment['texts'][0]['value']
                    break
                    
            return {"name": name, "function": function}
    except Exception as e:
        print(f"Failed to fetch UniProt details for {uniprot_id}: {e}")
        return {"name": "Unknown", "function": "Not available."}

def generate_blog_post(family, bait_id, bait_org, bait_seq, orphan_id, orphan_org, orphan_seq, dist):
    bait_details = get_uniprot_details(bait_id)
    orphan_details = get_uniprot_details(orphan_id)

    date_str = datetime.now().strftime("%Y-%m-%d")
    slug = f"mining-{family.lower().replace(' ', '-')}-{bait_id.lower()}-{orphan_id.lower()}-dark-proteome"
    filepath = f"/Users/jonatas/code/ideia.me/_posts/{date_str}-{slug}.md"
    
    # Match organism name by removing strains (usually after the first two words)
    org_base = " ".join(orphan_org.split()[:2])
    
    from organism_dict import enrichment_data
    rich_text = enrichment_data.get(org_base, None)
    
    if rich_text:
        discovery_sentence = f"Our search revealed an entirely uncharacterized protein (`{orphan_id}`) in *{orphan_org}*—{rich_text}."
    else:
        discovery_sentence = f"Our search revealed an entirely uncharacterized protein (`{orphan_id}`) in *{orphan_org}*."

    content = f"""---
layout: post
title: "Unearthing {family}: Exploring the Dark Proteome of Extreme Ecosystems!"
date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
categories: [bioinformatics, pgvector, machine-learning, structural-biology, pgbio]
---

As the `pg_bio` autonomous night pipeline continues its exciting sweep of the dark proteome, we set our sights on an incredible protein family: **{family}**! By bypassing months of wet-lab work, we are uncovering hidden secrets of nature using the immense power of native PostgreSQL multiomics engines scanning millions of vectors in milliseconds.

<!--more-->

Our SQL engine scanned the embedding space and found a high-confidence structural match that bridges two completely different biological worlds. We found an uncharacterized orphan protein that exhibits an almost identical 3D fold to a known, well-studied bait!

## The Bait: {bait_details['name']} ({bait_id})
To understand the magnitude of this discovery, we first must look at the known bait protein from *{bait_org}*. 
**What does it do?** 
{bait_details['function']}

This specific enzymatic function is crucial to its ecosystem. But what happens when we search the vast, uncharted territories of the database for something structurally similar?

## The Discovery: A Hidden Orphan in *{orphan_org}*
{discovery_sentence} Despite its label as "uncharacterized", its vector embeddings tell a different story!  

The structural similarity implies a massive evolutionary divergence or a conserved function adapted to a completely new environment. Could this extremophile or unique organism be harboring a more robust, efficient version of the enzyme? 

### Practical Applications & Impact
What does this mean for the real world? Proteins in the **{family}** family have massive potential in industrial biotechnology, bioremediation, medicine, and synthetic biology. By finding a novel version of this protein in *{orphan_org}*, we might have just discovered a variant that operates at extreme temperatures, pH levels, or with higher catalytic efficiency! This is the power of mining the dark proteome.

---

## The Math & The Pipeline
Using our newly built UniProt SQL Foreign Data Wrapper (`bio_search_uniprot`), we dynamically enriched the raw vector search directly inside the database:

| Category | Known Bait | Orphan Discovery |
| :--- | :--- | :--- |
| **UniProt ID** | `{bait_id}` | `{orphan_id}` |
| **Organism** | *{bait_org}* | *{orphan_org}* |
| **Status** | Characterized | Uncharacterized |
| **Cosine Distance** | - | **{dist:.4f}** |

*Note: A distance of {dist:.4f} means the 3D backbone is mathematically incredibly similar!*

### Interactive 3Dmol.js Preview
Dive into the structures below! *Tip: Double-click either 3D viewer to lock their cameras together for synchronized rotation, and click any fragment to automatically highlight the matching residue on the opposite protein!*

<div style="display: flex; justify-content: space-between; gap: 20px;">
  <div style="flex: 1;">
    <h4>Bait: {bait_id} ({bait_org})</h4>
    <div style="height: 400px; width: 100%; position: relative;" class="viewer_3Dmoljs" data-href="/assets/models/AF-{bait_id}-F1-model_v4.pdb" data-backgroundcolor="0xffffff" data-style="cartoon:color=cyan"></div>
  </div>
  <div style="flex: 1;">
    <h4>Discovery: {orphan_id} ({orphan_org})</h4>
    <div style="height: 400px; width: 100%; position: relative;" class="viewer_3Dmoljs" data-href="/assets/models/AF-{orphan_id}-F1-model_v4.pdb" data-backgroundcolor="0xffffff" data-style="cartoon:color=magenta"></div>
  </div>
</div>

### The SQL Query
This discovery was completely automated natively in PostgreSQL using our custom Z-Order indexing and the new UniProt SRF:

```sql
WITH closest AS (
    SELECT uniprot_id, name, embedding,
           (embedding <=> (SELECT embedding FROM proteins WHERE uniprot_id = '{bait_id}')) as dist
    FROM proteins
    WHERE name ILIKE '%uncharacterized%'
    ORDER BY dist ASC LIMIT 1
)
SELECT c.uniprot_id, c.dist, u.organism
FROM closest c
CROSS JOIN LATERAL bio_search_uniprot('accession:' || c.uniprot_id) u;
```


{{% include pg_bio_promo.md %}}
"""
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written {filepath}")


def run_loop():
    print("Starting Night Discoveries continuous loop indefinitely...")
    import random
    while True:
        random.shuffle(FAMILIES)
        for family in FAMILIES:
            now = datetime.now()
            print(f"[{now.strftime('%H:%M:%S')}] Scanning for {family}...")
            
            try:
                with psycopg.connect(DB_URI, autocommit=True) as conn:
                    with conn.cursor() as cur:
                        # 1. Get Bait
                        cur.execute("SELECT uniprot_id, name, sequence, embedding FROM proteins WHERE name ILIKE %s LIMIT 1", (f"%{family}%",))
                        bait = cur.fetchone()
                        if not bait:
                            continue
                        b_id, b_name, b_seq, b_emb = bait
                        
                        # 0. Get already discovered orphans
                        cur.execute("SELECT orphan_id FROM orphan_discoveries")
                        discovered_orphans = [row[0] for row in cur.fetchall()]
                        
                        # 2. Get Orphan
                        if discovered_orphans:
                            cur.execute("""
                                SELECT uniprot_id, name, sequence, (embedding <=> %s) as dist
                                FROM proteins
                                WHERE name ILIKE %s
                                  AND uniprot_id <> ALL(%s)
                                ORDER BY embedding <=> %s ASC
                                LIMIT 5;
                            """, (b_emb, '%uncharacterized%', discovered_orphans, b_emb))
                        else:
                            cur.execute("""
                                SELECT uniprot_id, name, sequence, (embedding <=> %s) as dist
                                FROM proteins
                                WHERE name ILIKE %s
                                ORDER BY embedding <=> %s ASC
                                LIMIT 5;
                            """, (b_emb, '%uncharacterized%', b_emb))
                            
                        orphans = cur.fetchall()
                        if not orphans:
                            print(f"No close orphan found for {family}.")
                            continue
                            
                        # Calculate baseline from top 5 (if we have at least 1)
                        baseline_dist = sum(o[3] for o in orphans) / len(orphans)
                        o_id, o_name, o_seq, o_dist = orphans[0]
                        
                        print(f"[{family}] Top dist: {o_dist:.4f} | Top {len(orphans)} Avg Baseline: {baseline_dist:.4f}")
                        
                        # Only consider it relevant if distance is very near to zero (<= 0.06)
                        if o_dist > 0.06:
                            print(f"Distance {o_dist:.4f} is not near enough to zero (must be <= 0.06). Skipping.")
                            continue
                            
                        # 3. Enrich Bait via SRF
                        cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{b_id}",))
                        b_org_res = cur.fetchone()
                        b_org = b_org_res[0] if b_org_res else "Unknown"
                        
                        # 4. Enrich Orphan via SRF
                        cur.execute("SELECT organism FROM bio_search_uniprot(%s) LIMIT 1", (f"accession:{o_id}",))
                        o_org_res = cur.fetchone()
                        o_org = o_org_res[0] if o_org_res else "Unknown"
                        
                        # 5. Fetch PDBs
                        b_ok = fetch_pdb_or_fold(cur, b_id, b_seq, f"/Users/jonatas/code/ideia.me/assets/models/AF-{b_id}-F1-model_v4.pdb")
                        o_ok = fetch_pdb_or_fold(cur, o_id, o_seq, f"/Users/jonatas/code/ideia.me/assets/models/AF-{o_id}-F1-model_v4.pdb")
                        
                        if b_ok and o_ok:
                            # Limit 3 blog posts per day
                            import glob
                            today_str = datetime.now().strftime("%Y-%m-%d")
                            generated_today = len(glob.glob(f"/Users/jonatas/code/ideia.me/_posts/{today_str}-*.md"))
                            
                            if generated_today < 3:
                                generate_blog_post(family, b_id, b_org, b_seq, o_id, o_org, o_seq, o_dist)
                                print(f"Generated blog post for {family}.")
                            else:
                                print(f"Max 3 blog posts per day reached ({generated_today} so far). Skipping post generation for {family}.")
                            
                            # Calculate confidence score (cosine distance is typically 0 to 2, 0 is identical)
                            confidence_score = max(0.0, 1.0 - (o_dist / 2.0))
                            
                            # Insert into the database catalog (we always catalog discoveries < 0.06)
                            cur.execute('''
                                INSERT INTO orphan_discoveries 
                                (family_name, bait_id, bait_organism, orphan_id, orphan_organism, vector_distance, status, confidence_score)
                                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                            ''', (family, b_id, b_org, o_id, o_org, o_dist, 'Pending', confidence_score))
                            conn.commit()
                            
                            print(f"Successfully cataloged {family}. Sleeping for 15 minutes...")
                            time.sleep(900) # Sleep 15 minutes between successful hits
    
                        else:
                            print(f"Failed to fetch PDBs for {family}. Skipping.")
                            
            except Exception as e:
                print(f"Error on {family}: {e}")
                time.sleep(60)

if __name__ == "__main__":
    run_loop()
