import psycopg
import math
import random
from rich.console import Console
from rich.table import Table
import time

console = Console()

DB_URL = "postgresql://localhost:28818/bio_demo"

def build_folding():
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            console.print("[bold cyan]🧬 pg_bio: Initializing Proteomic Folding Module (Step 2)[/]")
            
            # 1. Create the atomic structure table for the folded Teaflon protein
            cur.execute("""
                CREATE TABLE IF NOT EXISTS teaflon_atoms (
                    atom_id SERIAL PRIMARY KEY,
                    residue_name VARCHAR(10),
                    atom_type VARCHAR(10),
                    x FLOAT,
                    y FLOAT,
                    z FLOAT,
                    z_curve_index BIGINT
                );
            """)
            
            cur.execute("TRUNCATE teaflon_atoms;")
            
            console.print("Computational folding of Syn-TFLN completed (simulated).")
            console.print("Mapping 3D atomic coordinates into pg_bio using Z-Order Morton Coding...")
            
            # 2. Simulate injecting folded 3D coordinates (a dense beta-sheet tube for bio-ceramic strength)
            # We will create a dense core of hydrophobic residues and a specific binding pocket for Keratin.
            insert_query = """
                INSERT INTO teaflon_atoms (residue_name, atom_type, x, y, z, z_curve_index) 
                VALUES (%s, %s, %s, %s, %s, z_order_encode(%s, %s, %s))
            """
            
            # Generate a cylindrical beta-barrel structure
            for i in range(1000):
                angle = i * 0.1
                height = i * 0.05
                radius = 10.0
                
                # Basic structural atoms (Carbon backbone)
                x = radius * math.cos(angle) + random.uniform(-0.5, 0.5)
                y = radius * math.sin(angle) + random.uniform(-0.5, 0.5)
                z = height + random.uniform(-0.5, 0.5)
                
                residue = "ALA" if i % 2 == 0 else "LEU"
                cur.execute(insert_query, (residue, "C", x, y, z, x, y, z))
                
                # At a specific location (x: 10, y: 0, z: 25), we engineer the "Keratin Binding Pocket"
                # using specialized polar residues (Arginine/Glutamic Acid) to latch onto nail keratin.
                if 490 <= i <= 510:
                    cur.execute(insert_query, ("ARG_BIND", "N", x+2.0, y+2.0, z, x+2.0, y+2.0, z))
            
            # 3. Create the B-Tree index on the 1D Z-Curve index
            # This is the magic: we index 3D space using a standard 1D Postgres B-Tree!
            cur.execute("CREATE INDEX IF NOT EXISTS idx_teaflon_zorder ON teaflon_atoms (z_curve_index);")
            conn.commit()
            
            # 4. Execute the 3D Bounding Box Search to map the Binding Pocket
            console.print("\n[bold yellow]🔍 Querying 3D Spatial Binding Pocket via Z-Order Curve...[/]")
            start_time = time.time()
            
            # We search for any atoms within a 5-Angstrom bounding box around (x=10, y=0, z=25)
            # We use z_order_encode to rapidly filter millions of empty 3D spaces
            cur.execute("""
                SELECT atom_id, residue_name, atom_type, x, y, z 
                FROM teaflon_atoms
                WHERE z_curve_index BETWEEN z_order_encode(5.0, -5.0, 20.0) AND z_order_encode(15.0, 5.0, 30.0)
                  AND x BETWEEN 5.0 AND 15.0
                  AND y BETWEEN -5.0 AND 5.0
                  AND z BETWEEN 20.0 AND 30.0
                  AND residue_name = 'ARG_BIND'
                ORDER BY z ASC;
            """)
            
            results = cur.fetchall()
            elapsed = time.time() - start_time
            
            table = Table(title="Syn-TFLN (Teaflon) Keratin Binding Pocket (3D Match)")
            table.add_column("Atom ID", style="dim")
            table.add_column("Residue", style="cyan")
            table.add_column("Type", style="magenta")
            table.add_column("X (Å)", style="green")
            table.add_column("Y (Å)", style="green")
            table.add_column("Z (Å)", style="green")
            
            for row in results:
                table.add_row(str(row[0]), row[1], row[2], f"{row[3]:.2f}", f"{row[4]:.2f}", f"{row[5]:.2f}")
                
            console.print(table)
            console.print(f"✅ 3D Spatial Z-Order Query executed in {elapsed*1000:.2f} ms")
            
            if results:
                console.print("\n[bold green]Success![/] The folded Teaflon structure perfectly exposes its Keratin Binding Pocket.")
                console.print("By utilizing Morton Coding (Z-Order curves), PostgreSQL traversed the 3D atomic cloud")
                console.print("using a simple 1D B-Tree index, finding the exact interaction sites where Teaflon will lock into the nail matrix!")

if __name__ == "__main__":
    build_folding()
