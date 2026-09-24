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
from rich.progress import track

DB_URI = "postgresql://localhost:28818/bio_demo"
console = Console()

def batch_discover(keyword: str, max_distance: float = 0.35):
    console.print(f"[bold green]🧬 pg_bio: Batch Dark Proteome De-Orphanizer[/bold green]")
    console.print(f"Targeting Family: [bold cyan]{keyword}[/bold cyan] (Distance threshold < {max_distance})\n")
    
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # 1. Fetch all known targets in this family
            cur.execute("SELECT uniprot_id, name FROM proteins WHERE name ILIKE %s", (f"%{keyword}%",))
            targets = cur.fetchall()
            
            if not targets:
                console.print(f"[red]No characterized proteins found for keyword '{keyword}'.[/red]")
                return
                
            console.print(f"Found [yellow]{len(targets)}[/yellow] known '{keyword}' proteins to use as bait...")
            
            discoveries = []
            
            # 2. Loop over targets and do KNN search against uncharacterized proteins
            # We use LIMIT 1 for maximum HNSW speed
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
                ORDER BY p.embedding <=> t.emb ASC
                LIMIT 1;
            """
            
            start_time = time.time()
            
            for t_id, t_name in track(targets, description="Mining Dark Proteome..."):
                cur.execute(query, (t_id,))
                res = cur.fetchone()
                
                if res:
                    u_id, u_name, dist = res
                    if dist <= max_distance:
                        discoveries.append({
                            "target_id": t_id,
                            "target_name": t_name.split(" OS=")[0],
                            "orphan_id": u_id,
                            "orphan_name": u_name.split(" OS=")[0],
                            "distance": dist
                        })
            
            elapsed = time.time() - start_time
            
            # 3. Report Results
            console.print(f"\n✅ Batch Mining Complete in [bold blue]{elapsed:.2f} seconds[/bold blue]!\n")
            
            if not discoveries:
                console.print(f"[bold yellow]No structural homologues found in the Dark Proteome for '{keyword}'.[/bold yellow]")
                return
                
            table = Table(title=f"Novel Orphan Discoveries for '{keyword}'")
            table.add_column("Known Target (Bait)", style="cyan")
            table.add_column("Orphan Discovery", style="magenta")
            table.add_column("Vector Distance", style="green", justify="right")
            
            # Sort by closest distance
            discoveries.sort(key=lambda x: x["distance"])
            
            for d in discoveries[:20]: # Show top 20
                table.add_row(
                    f"{d['target_id']}\n[dim]{d['target_name']}[/dim]",
                    f"{d['orphan_id']}\n[dim]{d['orphan_name']}[/dim]",
                    f"{d['distance']:.4f}"
                )
                
            console.print(table)

if __name__ == "__main__":
    kw = sys.argv[1] if len(sys.argv) > 1 else "CRISPR"
    # The default distance is 0.35 to cast a slightly wider net
    batch_discover(kw, max_distance=0.35)
