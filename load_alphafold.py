# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "psycopg>=3.1",
#     "requests",
#     "urllib3"
# ]
# ///

import psycopg
import requests
import argparse
import sys
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DB_URI = "postgresql://localhost:28818/bio_demo"

# Configure robust requests session
session = requests.Session()
retry = Retry(connect=5, read=5, backoff_factor=1.0, status_forcelist=[429, 500, 502, 503, 504])
adapter = HTTPAdapter(max_retries=retry)
session.mount('http://', adapter)
session.mount('https://', adapter)

def get_uniprot_ids(tax_id, limit=100, reviewed_only=True):
    """Fetch UniProt IDs for a given organism taxonomy ID."""
    review_flag = "reviewed:true" if reviewed_only else "reviewed:false"
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=list&query=({review_flag})+AND+(taxonomy_id:{tax_id})"
    
    print(f"Fetching UniProt IDs from: {url}")
    resp = session.get(url, stream=True)
    if resp.status_code != 200:
        print(f"Failed to fetch from UniProt. Status: {resp.status_code}")
        return []
        
    ids = []
    for line in resp.iter_lines(decode_unicode=True):
        if line:
            ids.append(line.strip())
            if limit > 0 and len(ids) >= limit:
                break
    return ids

def fetch_alphafold_data(uniprot_id):
    """Fetch 3D structure from AlphaFold DB API."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    api_resp = session.get(api_url)
    
    if api_resp.status_code != 200:
        return None, []
        
    data = api_resp.json()
    if not data:
        return None, []
        
    prediction = data[0]
    version = prediction.get("latestVersion", 4)
    sequence = prediction.get("sequence", "")
    
    pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.pdb"
    pdb_resp = session.get(pdb_url)
    
    if pdb_resp.status_code != 200:
        return sequence, []
        
    atoms = []
    for line in pdb_resp.text.splitlines():
        if line.startswith("ATOM  "):
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                # AlphaFold stores pLDDT confidence in the B-factor column
                plddt = float(line[60:66].strip())
                
                atoms.append({
                    "name": atom_name, 
                    "x": x, "y": y, "z": z, 
                    "plddt": plddt
                })
            except ValueError:
                continue
                
    return sequence, atoms

def load_alphafold_to_db(tax_id, limit):
    print(f"--- AlphaFold Loader: Taxonomy ID {tax_id} ---")
    
    # 1. Get IDs
    uniprot_ids = get_uniprot_ids(tax_id, limit, reviewed_only=True)
    if not uniprot_ids:
        print("No UniProt IDs found or API failed.")
        return
        
    print(f"Found {len(uniprot_ids)} UniProt entries. Connecting to DB...")
    
    # 2. Connect and Load
    inserted_proteins = 0
    inserted_atoms = 0
    
    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            for uid in uniprot_ids:
                # Check if exists to allow resuming
                cur.execute("SELECT 1 FROM proteins WHERE uniprot_id = %s", (uid,))
                if cur.fetchone():
                    print(f"Skipping {uid} (Already in database)")
                    continue
                
                print(f"Fetching AlphaFold structure for {uid}...", end=" ", flush=True)
                seq, atoms = fetch_alphafold_data(uid)
                
                if not seq:
                    print("Failed (Not found in AFDB)")
                    continue
                    
                # Insert protein and generate embedding
                cur.execute(
                    """
                    INSERT INTO proteins (uniprot_id, name, sequence, embedding) 
                    VALUES (%s, %s, %s, get_esm_embedding(%s))
                    ON CONFLICT (uniprot_id) DO NOTHING
                    RETURNING uniprot_id;
                    """,
                    (uid, f"AlphaFold {uid}", seq, seq)
                )
                
                if not cur.fetchone():
                    print("Skipped (Conflict)")
                    continue
                    
                inserted_proteins += 1
                
                # Insert atoms
                if atoms:
                    print(f"Success ({len(atoms)} atoms). Inserting...", end="")
                    records = []
                    for a in atoms:
                        # Serialize safely
                        coord_str = f'{{"x":{a["x"]},"y":{a["y"]},"z":{a["z"]},"name":"{a["name"]}"}}'
                        records.append((uid, coord_str))
                        
                    with cur.copy("COPY protein_atoms (uniprot_id, coord) FROM STDIN") as copy:
                        for r in records:
                            copy.write_row(r)
                    inserted_atoms += len(atoms)
                    print(" Done.")
                else:
                    print(" Done (No atoms found).")
                    
    print("\n--- Summary ---")
    print(f"Inserted Proteins: {inserted_proteins}")
    print(f"Inserted Atoms: {inserted_atoms}")
    print("Database is up to date!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load AlphaFold structures directly into pg_bio")
    parser.add_argument("--tax-id", type=int, default=9606, help="Organism Taxonomy ID (e.g. 9606 for Human)")
    parser.add_argument("--limit", type=int, default=10, help="Max proteins to fetch (0 for all)")
    args = parser.parse_args()
    
    load_alphafold_to_db(args.tax_id, args.limit)
