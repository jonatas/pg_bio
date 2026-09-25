# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "rich",
#     "requests",
# ]
# ///

import psycopg
import sys
import argparse
import zlib
import urllib.request
import json
import os
from rich.console import Console

console = Console()
DB_URI = "postgresql://localhost:28818/bio_demo"
CHECKPOINT_FILE = ".ingest_checkpoint.json"

# We use a mocked embedding function since we can't run a 10GB ESM model in realtime
def get_esm_embedding(sequence: str) -> list[float]:
    import hashlib
    # Deterministic pseudo-embedding (1280 dimensions)
    h = hashlib.sha256(sequence.encode()).digest()
    base = [float(b) / 255.0 for b in h]
    return (base * 40)[:1280]

def init_db():
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute("CREATE EXTENSION IF NOT EXISTS vector;")
            cur.execute("""
                CREATE TABLE IF NOT EXISTS proteins (
                    uniprot_id VARCHAR(15) PRIMARY KEY,
                    name TEXT NOT NULL,
                    sequence TEXT NOT NULL,
                    embedding halfvec(1280) NOT NULL
                );
            """)

def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, 'r') as f:
            return json.load(f)
    return {}

def save_checkpoint(dataset, count):
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump({dataset: count}, f)

def stream_fasta(dataset: str):
    if dataset == "sprot":
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    elif dataset.startswith("trembl_"):
        division = dataset.split("_")[1]
        taxonomies = {
            "archaea": "2157",
            "bacteria": "2",
            "fungi": "4751",
            "viruses": "10239",
            "plants": "3193"
        }
        if division not in taxonomies:
            console.print(f"[red]Unknown TrEMBL division: {division}. Try archaea, bacteria, fungi, viruses, or plants.[/red]")
            return
        tax_id = taxonomies[division]
        url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=true&query=reviewed:false+AND+taxonomy_id:{tax_id}"
    else:
        console.print(f"[red]Unknown dataset: {dataset}. Use 'sprot' or 'trembl_archaea', 'trembl_viruses', etc.[/red]")
        return
        
    console.print(f"[cyan]Streaming dataset '{dataset}' from {url}...[/cyan]")
    
    init_db()
    chk = load_checkpoint()
    skip_count = chk.get(dataset, 0)
    
    if skip_count > 0:
        console.print(f"[yellow]Resuming {dataset} from sequence {skip_count}...[/yellow]")

    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    
    batch = []
    count = 0
    BATCH_SIZE = 10000
    
    with urllib.request.urlopen(req) as response:
        dctx = zlib.decompressobj(16 + zlib.MAX_WBITS)
        
        current_id = ""
        current_name = ""
        current_seq = []
        
        # Read stream in chunks
        for chunk in iter(lambda: response.read(8192 * 4), b""):
            text = dctx.decompress(chunk).decode('utf-8')
            lines = text.split('\n')
            
            for line in lines:
                if not line: continue
                
                if line.startswith('>'):
                    if current_id:
                        count += 1
                        if count > skip_count:
                            seq_str = "".join(current_seq)
                            batch.append((current_id, current_name, seq_str, f"[{','.join(map(str, get_esm_embedding(seq_str)))}]"))
                            
                            if len(batch) >= BATCH_SIZE:
                                insert_batch(batch)
                                save_checkpoint(dataset, count)
                                console.print(f"[+] Inserted {count} proteins from {dataset}...")
                                batch = []
                                
                    parts = line.split('|')
                    if len(parts) >= 3:
                        current_id = parts[1]
                        current_name = parts[2].strip()
                    else:
                        current_id = "UNKNOWN"
                        current_name = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
                    
        # Final flush
        if current_id:
            count += 1
            if count > skip_count:
                seq_str = "".join(current_seq)
                batch.append((current_id, current_name, seq_str, f"[{','.join(map(str, get_esm_embedding(seq_str)))}]"))
                
        if batch:
            insert_batch(batch)
            save_checkpoint(dataset, count)
            console.print(f"[+] Final flush. Inserted {count} total proteins.")
            
    console.print(f"🎉 INGESTION COMPLETE for {dataset}!")
    
    # Rebuild index
    console.print("[cyan]Rebuilding HNSW halfvec index (this may take a few minutes)...[/cyan]")
    with psycopg.connect(DB_URI) as conn:
        conn.autocommit = True
        with conn.cursor() as cur:
            cur.execute("DROP INDEX IF EXISTS idx_protein_embedding;")
            cur.execute("CREATE INDEX idx_protein_embedding ON proteins USING hnsw (embedding halfvec_cosine_ops);")
    console.print("✅ Index built successfully.")

def insert_batch(batch):
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.executemany(
                "INSERT INTO proteins (uniprot_id, name, sequence, embedding) VALUES (%s, %s, %s, %s) ON CONFLICT (uniprot_id) DO NOTHING",
                batch
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stream UniProt/TrEMBL slices into pg_bio")
    parser.add_argument("--dataset", type=str, default="sprot", help="sprot, trembl_archaea, trembl_viruses, trembl_bacteria, trembl_fungi, etc.")
    args = parser.parse_args()
    stream_fasta(args.dataset)
