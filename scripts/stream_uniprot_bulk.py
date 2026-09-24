# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "requests",
# ]
# ///

import psycopg
import requests
import zlib
import json
import os
import sys
import time

DB_URI = "postgresql://localhost:28818/bio_demo"

# Swiss-Prot (570k proteins, ~90MB compressed) -> Safe for laptop testing
SWISSPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
# TrEMBL (250M proteins, ~40GB compressed) -> WILL FILL 747GB HARD DRIVE
TREMBL_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"

CHECKPOINT_FILE = ".ingest_checkpoint.json"
BATCH_SIZE = 10000

def get_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r") as f:
            return json.load(f).get("last_uniprot_id")
    return None

def save_checkpoint(uniprot_id):
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump({"last_uniprot_id": uniprot_id}, f)

def flush_batch(conn, batch):
    """Bulk insert a batch of proteins using COPY for maximum performance."""
    if not batch: return
    
    with conn.cursor() as cur:
        with cur.copy("COPY proteins (uniprot_id, name, sequence, embedding) FROM STDIN") as copy:
            for pid, name, seq in batch:
                # We offload the embedding generation to the database via a trigger or default, 
                # but COPY bypasses defaults. Instead, we can write a raw string that PostgreSQL parses.
                # However, COPYing function calls is tricky. 
                pass
        
        # ACTUALLY, for embedding calculation inside Postgres, execute_batch is safer than COPY 
        # when we need to evaluate `get_esm_embedding(seq)`.
        cur.executemany("""
            INSERT INTO proteins (uniprot_id, name, sequence, embedding) 
            VALUES (%s, %s, %s, get_esm_embedding(%s)::vector(1280))
            ON CONFLICT (uniprot_id) DO NOTHING;
        """, [(p[0], p[1], p[2], p[2]) for p in batch])
    conn.commit()

def stream_fasta(url):
    print(f"📡 Initiating HTTP Stream from: {url}")
    print("   (Data is decompressed entirely in-memory. Nothing is saved to disk!)")
    
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    
    buffer = ""
    for chunk in response.iter_content(chunk_size=1024 * 1024): # 1MB chunks
        if chunk:
            text = d.decompress(chunk).decode('utf-8')
            buffer += text
            
            # Yield full lines
            while "\n" in buffer:
                line, buffer = buffer.split("\n", 1)
                yield line
                
    # Yield remaining
    if buffer:
        yield buffer

def run_ingestion(url):
    checkpoint_id = get_checkpoint()
    skipping = checkpoint_id is not None
    
    if skipping:
        print(f"🔄 Resuming from checkpoint: {checkpoint_id}. Scanning stream (this may take a while)...")
    else:
        print("🚀 Starting fresh ingestion pipeline...")

    conn = psycopg.connect(DB_URI, autocommit=False)
    
    batch = []
    current_id = None
    current_name = None
    current_seq = []
    
    total_inserted = 0
    start_time = time.time()
    
    try:
        for line in stream_fasta(url):
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                # Save previous protein
                if current_id:
                    seq = "".join(current_seq)
                    
                    if skipping:
                        if current_id == checkpoint_id:
                            print(f"✅ Found checkpoint '{checkpoint_id}'. Resuming database inserts!")
                            skipping = False
                    else:
                        batch.append((current_id, current_name, seq))
                        
                        if len(batch) >= BATCH_SIZE:
                            flush_batch(conn, batch)
                            total_inserted += len(batch)
                            save_checkpoint(batch[-1][0])
                            
                            elapsed = time.time() - start_time
                            rate = total_inserted / elapsed
                            print(f"   [+] Inserted {total_inserted} proteins... ({rate:.0f} prot/sec)")
                            batch = []

                # Parse new header
                # e.g., >tr|A0A024R1R8|A0A024R1R8_HUMAN HCG1982144
                parts = line.split("|")
                current_id = parts[1] if len(parts) > 1 else "UNKNOWN"
                current_name = parts[2].split(" OS=")[0] if len(parts) > 2 else "Unknown"
                current_seq = []
            else:
                current_seq.append(line)
                
        # Flush final batch
        if batch and not skipping:
            flush_batch(conn, batch)
            total_inserted += len(batch)
            save_checkpoint(batch[-1][0])
            
        print(f"\n🎉 INGESTION COMPLETE! Successfully loaded {total_inserted} new proteins.")
        
    except KeyboardInterrupt:
        print("\n🛑 Pipeline paused by user. Checkpoint saved. Run again to resume.")
    except Exception as e:
        print(f"\n❌ Error during stream: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    # Choose your destiny:
    target_url = SWISSPROT_URL # Change to TREMBL_URL if you dare
    run_ingestion(target_url)
