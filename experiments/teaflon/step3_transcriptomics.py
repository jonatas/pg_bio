import psycopg
import random
from rich.console import Console
from rich.table import Table
import time

console = Console()

DB_URL = "postgresql://localhost:28818/bio_demo"

# Gene ID for Syn-TFLN (Teaflon)
TEAFLON_GENE_ID = 9999

def build_transcriptomics():
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            console.print("[bold cyan]🧬 pg_bio: Initializing Transcriptomics Module (Step 3)[/]")
            
            # 1. Create the Single-Cell RNA CSR Table
            cur.execute("""
                CREATE TABLE IF NOT EXISTS single_cell_expression (
                    cell_id SERIAL PRIMARY KEY,
                    cell_type VARCHAR(50),
                    expressed_gene_ids INT[],
                    expression_counts REAL[]
                );
            """)
            
            cur.execute("TRUNCATE single_cell_expression;")
            
            console.print("Simulating Single-Cell RNA sequencing for 10,000 human cells...")
            
            # 2. Inject simulated data
            # We simulate Keratinocytes (Nail cells) which should express Teaflon at the target dosage
            # and Osteoblasts (Bone cells) which should have near-zero expression (since we targeted the nail promoter).
            
            insert_query = """
                INSERT INTO single_cell_expression (cell_type, expressed_gene_ids, expression_counts) 
                VALUES (%s, %s, %s)
            """
            
            for i in range(1, 10001):
                is_keratinocyte = (i <= 5000)
                cell_type = "Keratinocyte (Nail)" if is_keratinocyte else "Osteoblast (Bone)"
                
                # Standard housekeeping genes that every cell expresses
                gene_ids = [101, 102, 103]
                expr_counts = [random.uniform(50, 100), random.uniform(20, 40), random.uniform(200, 300)]
                
                # Teaflon Gene Dosage Simulation (Target: 15 to 50 copies)
                teaflon_expr = 0.0
                if is_keratinocyte:
                    # Normally distributed dosage around 35 copies
                    teaflon_expr = max(0.0, random.gauss(35.0, 5.0))
                else:
                    # Accidental leaky expression in bone cells (should be near zero)
                    teaflon_expr = max(0.0, random.gauss(0.5, 0.2)) if random.random() > 0.9 else 0.0
                
                if teaflon_expr > 0:
                    gene_ids.append(TEAFLON_GENE_ID)
                    expr_counts.append(teaflon_expr)
                
                cur.execute(insert_query, (cell_type, gene_ids, expr_counts))
            
            conn.commit()
            
            # 3. Execute the Dosage Analysis Query using PostgreSQL Array Unnesting
            console.print("\n[bold yellow]🔍 Analyzing Syn-TFLN (Teaflon) Expression Dosage via CSR Arrays...[/]")
            start_time = time.time()
            
            # We UNNEST the sparse arrays in parallel with WITH ORDINALITY to map genes to their exact expression counts.
            # Then we aggregate to find the average dosage per cell type.
            cur.execute(f"""
                WITH unnested_expression AS (
                    SELECT 
                        cell_id, 
                        cell_type,
                        gene_ids.gene_id,
                        expr_counts.expr
                    FROM single_cell_expression e
                    CROSS JOIN LATERAL unnest(e.expressed_gene_ids) WITH ORDINALITY AS gene_ids(gene_id, idx1)
                    JOIN LATERAL unnest(e.expression_counts) WITH ORDINALITY AS expr_counts(expr, idx2)
                      ON gene_ids.idx1 = expr_counts.idx2
                    WHERE gene_ids.gene_id = {TEAFLON_GENE_ID}
                )
                SELECT 
                    cell_type,
                    COUNT(cell_id) as cells_expressing,
                    AVG(expr) as average_dosage_rna_copies
                FROM unnested_expression
                GROUP BY cell_type
                ORDER BY average_dosage_rna_copies DESC;
            """)
            
            results = cur.fetchall()
            elapsed = time.time() - start_time
            
            table = Table(title="Teaflon (Syn-TFLN) Single-Cell RNA Dosage Analysis")
            table.add_column("Tissue Cell Type", style="cyan")
            table.add_column("Cells Expressing", style="magenta")
            table.add_column("Average RNA Dosage (Target: 15-50)", style="green")
            
            for row in results:
                table.add_row(row[0], str(row[1]), f"{row[2]:.2f}")
                
            console.print(table)
            console.print(f"✅ Sparse Array (CSR) Query executed in {elapsed*1000:.2f} ms")
            
            # Validation Output
            keratinocyte_dosage = next((row[2] for row in results if "Keratinocyte" in row[0]), 0.0)
            if 15.0 <= keratinocyte_dosage <= 50.0:
                console.print(f"\n[bold green]Success![/] The average dosage in Nails is {keratinocyte_dosage:.2f}.")
                console.print("This is perfectly inside the Goldilocks Zone (15-50). The bio-ceramic will reinforce the nail without shattering it!")
            else:
                console.print("\n[bold red]Warning:[/] Dosage is outside the Goldilocks Zone!")

if __name__ == "__main__":
    build_transcriptomics()
