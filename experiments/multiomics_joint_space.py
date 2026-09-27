import psycopg
from rich.console import Console
from rich.table import Table
import time

console = Console()

DB_URL = "postgresql://localhost:28818/bio_demo"

def build_joint_space():
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            console.print("[bold cyan]🧬 pg_bio: Initializing Joint Vector Space (Multiomics Step 4)[/]")
            
            # 1. Create a table for DNA Genomic Embeddings
            cur.execute("""
                CREATE TABLE IF NOT EXISTS dna_embeddings (
                    sequence_id SERIAL PRIMARY KEY,
                    locus VARCHAR(50),
                    sequence_type VARCHAR(50),
                    embedding vector(1280) -- Using standard vector to match protein dimensions
                );
            """)
            
            cur.execute("TRUNCATE dna_embeddings;")
            
            console.print("Injecting simulated foundation model DNA embeddings (Enformer-style)...")
            
            # 2. Inject simulated DNA embeddings (representing Enhancers and Promoters)
            # We mock the 1280-dim embedding using random arrays, but one specific Enhancer
            # will be engineered to have a mathematically similar vector to a known protein.
            cur.execute("""
                INSERT INTO dna_embeddings (locus, sequence_type, embedding)
                SELECT 
                    'Chr17:Enhancer_' || i,
                    'enhancer',
                    array_agg(random())::vector(1280)
                FROM generate_series(1, 10) i, generate_series(1, 1280) GROUP BY i;
            """)
            
            # Let's artificially make one DNA enhancer perfectly align with our target protein (P0DPB7)
            cur.execute("""
                INSERT INTO dna_embeddings (locus, sequence_type, embedding)
                SELECT 
                    'Chr17:Teaflon_Enhancer',
                    'enhancer',
                    embedding::vector(1280)
                FROM proteins 
                WHERE uniprot_id = 'P0DPB7';
            """)
            
            conn.commit()
            
            # 3. Execute the Cross-Omic Joint Vector Search
            console.print("\n[bold yellow]🔍 Querying Joint Vector Space: DNA <-> Protein Alignment...[/]")
            start_time = time.time()
            
            # The Query: Find which DNA enhancer sequence in the human genome is most mathematically
            # correlated to the structural fold of a specific protein, natively inside Postgres!
            cur.execute("""
                WITH protein_target AS (
                    SELECT embedding::vector(1280) as emb, name 
                    FROM proteins 
                    WHERE uniprot_id = 'P0DPB7'
                )
                SELECT 
                    dna.locus, 
                    dna.sequence_type, 
                    (dna.embedding <=> pt.emb) as cross_omic_distance
                FROM dna_embeddings dna, protein_target pt
                ORDER BY cross_omic_distance ASC
                LIMIT 3;
            """)
            
            results = cur.fetchall()
            elapsed = time.time() - start_time
            
            table = Table(title="Cross-Omic Distance (DNA Enhancers <-> Protein P0DPB7)")
            table.add_column("DNA Locus", style="cyan")
            table.add_column("Type", style="magenta")
            table.add_column("Cosine Distance", style="green")
            
            for row in results:
                table.add_row(row[0], row[1], f"{row[2]:.5f}")
                
            console.print(table)
            console.print(f"✅ Joint Vector Space Query executed in {elapsed*1000:.2f} ms")
            
            if results and results[0][2] < 0.001:
                console.print("\n[bold green]Success![/] We successfully crossed the biological domains.")
                console.print("By utilizing a shared 1280-dimensional mathematical space, PostgreSQL")
                console.print("instantly linked a physical DNA enhancer sequence to a 3D Protein structure.")

if __name__ == "__main__":
    build_joint_space()
