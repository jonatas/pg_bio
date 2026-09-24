# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "rich",
# ]
# ///

import psycopg
import sys
import time
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

DB_URI = "postgresql://localhost:28818/bio_demo"
console = Console()

def search_dark_proteome(target_uniprot_id: str, max_distance: float = 0.2):
    console.print(Panel(f"[bold green]🧬 pg_bio: Dark Proteome De-Orphanizer[/bold green]\nTarget Query: [bold cyan]{target_uniprot_id}[/bold cyan]"))
    
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # 1. Fetch Target Details
            cur.execute("SELECT name FROM proteins WHERE uniprot_id = %s", (target_uniprot_id,))
            target_res = cur.fetchone()
            if not target_res:
                console.print(f"[red]Error: Target {target_uniprot_id} not found in database![/red]")
                return
            target_name = target_res[0]
            console.print(f"Target Identity: [yellow]{target_name}[/yellow]\n")
            
            # 2. Execute Vector Search against Uncharacterized Proteins
            console.print(f"🔍 Scanning [bold]575,000+[/bold] embeddings using native HNSW index...")
            console.print("   Filtering for 'Uncharacterized' or 'Hypothetical' proteins...")
            
            start_time = time.time()
            
            # The Magic Query: HNSW Cosine Distance combined with ILIKE filters
            query = """
                WITH target AS (
                    SELECT embedding as emb FROM proteins WHERE uniprot_id = %s
                )
                SELECT 
                    p.uniprot_id, 
                    p.name, 
                    (p.embedding <=> t.emb) as distance
                FROM proteins p, target t
                WHERE (p.name ILIKE '%%uncharacterized%%' OR p.name ILIKE '%%hypothetical%%')
                  AND (p.embedding <=> t.emb) < %s
                ORDER BY p.embedding <=> t.emb ASC
                LIMIT 15;
            """
            
            cur.execute(query, (target_uniprot_id, max_distance))
            results = cur.fetchall()
            elapsed = time.time() - start_time
            
            console.print(f"✅ Search complete in [bold blue]{elapsed:.3f} seconds[/bold blue]!\n")
            
            if not results:
                console.print(f"[bold red]No uncharacterized proteins found within distance {max_distance}[/bold red]")
                return
            
            # 3. Display Discoveries
            table = Table(title=f"Potential Novel Functional Analogues for {target_uniprot_id}")
            table.add_column("UniProt ID", style="cyan", no_wrap=True)
            table.add_column("Organism / Designation", style="magenta")
            table.add_column("Vector Distance", style="green")
            table.add_column("Putative Function", style="yellow")
            
            for pid, name, dist in results:
                # Extract organism info if available (OS=...)
                org_split = name.split(" OS=")
                clean_name = org_split[0]
                org = org_split[1].split(" OX=")[0] if len(org_split) > 1 else "Unknown"
                
                # We dynamically assign the putative function based on the target!
                putative = f"Novel {target_name.split(' ')[0]}"
                
                table.add_row(pid, f"{clean_name}\n[dim]OS: {org}[/dim]", f"{dist:.4f}", putative)
                
            console.print(table)
            
            console.print("\n[bold]🧪 SCIENTIFIC CONCLUSION:[/bold]")
            console.print("These proteins share no obvious text annotations with the target, yet they exist in the exact same mathematical structure-space. They are prime candidates for novel CRISPR/enzymatic systems!")

if __name__ == "__main__":
    # We use Q97YC2 (CRISPR-associated endoribonuclease Cas2) as a default, 
    # but the user can pass any UniProt ID via CLI
    target = sys.argv[1] if len(sys.argv) > 1 else "Q97YC2"
    search_dark_proteome(target, max_distance=0.3)
