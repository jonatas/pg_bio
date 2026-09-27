import psycopg
from rich.console import Console
from rich.table import Table
import time

console = Console()

DB_URL = "postgresql://localhost:28818/bio_demo"

def build_genomics():
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            console.print("[bold cyan]🧬 pg_bio: Initializing Genomics Module (Step 1)[/]")
            
            # 1. Create the human_genome table natively using Postgres Range Types
            cur.execute("""
                CREATE TABLE IF NOT EXISTS human_genome (
                    chromosome VARCHAR(2),
                    feature_type VARCHAR(50),
                    feature_name VARCHAR(100),
                    genomic_range INT4RANGE,
                    strand CHAR(1)
                );
            """)
            
            # 2. Build the GiST Index for 1D spatial overlap
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_genome_range 
                ON human_genome USING GIST (genomic_range);
            """)
            
            # Clean up old data if we rerun
            cur.execute("TRUNCATE human_genome;")
            
            # 3. Insert biological scenario data for Teaflon (Chromosome 17)
            console.print("Injecting Chromosome 17 annotations (KRT14 locus)...")
            
            # KRT14 Promoter (Target activation site for Nails)
            cur.execute("INSERT INTO human_genome VALUES ('17', 'promoter', 'KRT14_promoter', int4range(41490000, 41492000), '+');")
            
            # Tumor Suppressors (Dangerous to interrupt!)
            cur.execute("INSERT INTO human_genome VALUES ('17', 'tumor_suppressor', 'TP53', int4range(7668400, 7687500), '-');")
            cur.execute("INSERT INTO human_genome VALUES ('17', 'tumor_suppressor', 'BRCA1', int4range(43044294, 43125482), '-');")
            
            # Safe Harbor A (Too far from KRT14)
            cur.execute("INSERT INTO human_genome VALUES ('17', 'safe_harbor', 'SafeLocus_A', int4range(10000000, 10005000), '+');")
            
            # Safe Harbor B (Close to KRT14, safe to insert Teaflon)
            cur.execute("INSERT INTO human_genome VALUES ('17', 'safe_harbor', 'SafeLocus_B', int4range(41492100, 41496000), '+');")
            
            # Safe Harbor C (Overlaps with BRCA1, incredibly dangerous!)
            cur.execute("INSERT INTO human_genome VALUES ('17', 'safe_harbor', 'SafeLocus_C', int4range(43040000, 43050000), '+');")
            
            conn.commit()
            
            # 4. Execute the Teaflon Search Query
            console.print("\n[bold yellow]🔍 Searching for optimal Teaflon Genomic Insertion Site...[/]")
            start_time = time.time()
            
            cur.execute("""
                SELECT 
                    safe.feature_name AS safe_harbor_site,
                    safe.genomic_range AS insertion_coordinates,
                    promoter.feature_name AS nail_promoter
                FROM human_genome safe
                JOIN human_genome promoter 
                  ON safe.chromosome = promoter.chromosome
                WHERE 
                    safe.feature_type = 'safe_harbor'
                    AND promoter.feature_name = 'KRT14_promoter'
                    -- Check if the safe harbor is physically adjacent to the nail promoter (distance < 5000 base pairs)
                    AND GREATEST(
                        lower(safe.genomic_range) - upper(promoter.genomic_range),
                        lower(promoter.genomic_range) - upper(safe.genomic_range),
                        0
                    ) < 5000
                    -- Ensure our insertion DOES NOT overlap (&&) with any known tumor suppressors
                    AND NOT EXISTS (
                        SELECT 1 FROM human_genome cancer_genes
                        WHERE cancer_genes.feature_type = 'tumor_suppressor'
                        AND safe.genomic_range && cancer_genes.genomic_range
                    );
            """)
            
            results = cur.fetchall()
            elapsed = time.time() - start_time
            
            table = Table(title="Teaflon Bio-Ceramic Genomic Targeting Results")
            table.add_column("Safe Harbor", style="cyan")
            table.add_column("Coordinates (int4range)", style="green")
            table.add_column("Driving Promoter", style="magenta")
            
            for row in results:
                table.add_row(row[0], str(row[1]), row[2])
                
            console.print(table)
            console.print(f"✅ GiST Range Query executed in {elapsed*1000:.2f} ms")
            
            if results:
                console.print(f"\n[bold green]Success![/] We can insert the Syn-TFLN gene at [cyan]{results[0][0]}[/].")
                console.print("Because it sits adjacent to the KRT14 promoter, it will only activate in keratinocytes (nails).")
                console.print("And because we verified no GiST overlap with TP53/BRCA1, it will not disrupt tumor suppressors.")

if __name__ == "__main__":
    build_genomics()
