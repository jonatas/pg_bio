# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "psycopg>=3.1",
#     "requests",
#     "torch",
#     "transformers",
#     "tqdm",
# ]
# ///

import psycopg
import requests
import torch
from transformers import AutoTokenizer, EsmModel
import sys

DB_URI = "postgresql://localhost:28818/bio_demo"

# We use the extremely fast 8-Million parameter ESM-2 model. 
# It runs lightning fast locally but still captures deep evolutionary biology!
MODEL_NAME = "facebook/esm2_t6_8M_UR50D" 

def map_all_locally():
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"Loading AI Model ({MODEL_NAME}) on {device}...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmModel.from_pretrained(MODEL_NAME).to(device)
    model.eval()

    print("\nConnecting to bio_demo database...")
    try:
        conn = psycopg.connect(DB_URI, autocommit=True)
    except Exception as e:
        print(f"Failed to connect to Postgres. Is it running on 28818? Error: {e}")
        sys.exit(1)

    with conn.cursor() as cur:
        # Clear previous mocked data
        # We comment out TRUNCATE so the script can safely resume if it was interrupted!
        # cur.execute("TRUNCATE TABLE proteins CASCADE;")
        
        print("\nStreaming FULL Human Proteome from UniProt...")
        url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:9606)"
        response = requests.get(url, stream=True)
        
        proteins = []
        current_id, current_name, current_seq = None, None, []
        
        for line in response.iter_lines(decode_unicode=True):
            if not line: continue
            if line.startswith(">"):
                if current_id:
                    proteins.append((current_id, current_name, "".join(current_seq)))
                parts = line.split("|")
                current_id = parts[1] if len(parts) > 1 else "UNKNOWN"
                current_name = parts[2].split(" OS=")[0] if len(parts) > 2 else "Unknown"
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_id:
            proteins.append((current_id, current_name, "".join(current_seq)))
            
        print(f"✅ Downloaded {len(proteins)} real human proteins.")
        
        # Fetch already inserted proteins to avoid wasting GPU time on resumption
        cur.execute("SELECT uniprot_id FROM proteins")
        existing_ids = {row[0] for row in cur.fetchall()}
        
        proteins = [p for p in proteins if p[0] not in existing_ids]
        print(f"⏩ Resuming... {len(proteins)} proteins remaining to map.")
        
        print("\nStarting Local AI Mapping (Generating biological vector embeddings)...")
        
        # Process in batches to leverage GPU parallelism
        batch_size = 50
        total_batches = (len(proteins) + batch_size - 1) // batch_size
        
        for i in range(0, len(proteins), batch_size):
            batch = proteins[i:i+batch_size]
            seqs = [p[2] for p in batch]
            
            # Truncate sequences to 1022 amino acids to prevent Out-Of-Memory (OOM) on massive proteins
            seqs = [s[:1022] for s in seqs] 
            
            inputs = tokenizer(seqs, return_tensors="pt", padding=True, truncation=True)
            inputs = {k: v.to(device) for k, v in inputs.items()}
            
            with torch.no_grad():
                outputs = model(**inputs)
                # Mean pooling: Average the sequence tokens to get a single 320-dim vector per protein
                embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()
            
            records = []
            for j, (pid, name, seq) in enumerate(batch):
                emb_list = embeddings[j].tolist()
                records.append((pid, name, seq, emb_list))
            
            cur.executemany(
                "INSERT INTO proteins (uniprot_id, name, sequence, embedding) VALUES (%s, %s, %s, %s) ON CONFLICT (uniprot_id) DO NOTHING",
                records
            )
            
            if (i // batch_size) % 10 == 0:
                print(f"Mapped and inserted batch {(i // batch_size) + 1}/{total_batches} ({(i + len(batch))} proteins)...")
                
            # Clear the PyTorch MPS cache so it doesn't hoard memory and crash on massive proteomes!
            if torch.backends.mps.is_available():
                torch.mps.empty_cache()
                
    print("\n🎉 Mapped ALL proteins locally and saved AI embeddings to the database!")
    conn.close()

if __name__ == "__main__":
    map_all_locally()
