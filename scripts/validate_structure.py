# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
#     "requests",
#     "rich",
# ]
# ///

import os
import sys
import argparse
import requests
from rich.console import Console

console = Console()

def fetch_alphafold_cif(uniprot_id: str, out_dir: str = "data"):
    """Fetch AlphaFold model for a UniProt ID."""
    os.makedirs(out_dir, exist_ok=True)
    cif_path = os.path.join(out_dir, f"AF-{uniprot_id}-F1-model_v4.cif")
    
    if os.path.exists(cif_path):
        console.print(f"[green]Using cached AlphaFold model for {uniprot_id}[/green]")
        return cif_path
        
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
    console.print(f"Downloading AlphaFold model for {uniprot_id}...")
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(cif_path, "w") as f:
            f.write(response.text)
        return cif_path
    else:
        console.print(f"[red]Failed to fetch AlphaFold model for {uniprot_id} (Status: {response.status_code}). Note: Very large proteins (>2700aa) or Archaea orphans might not have a single AFDB file.[/red]")
        return None

def align_structures(uniprot_a: str, uniprot_b: str):
    cif_a = fetch_alphafold_cif(uniprot_a)
    cif_b = fetch_alphafold_cif(uniprot_b)
    
    if not cif_a or not cif_b:
        return
    
    # Set environment variable for headless rendering
    os.environ["PYOPENGL_PLATFORM"] = "osmesa"
    os.makedirs("output", exist_ok=True)
    
    console.print("[cyan]Initializing PyMOL for structural superposition...[/cyan]")
    import pymol
    pymol.pymol_argv = ["pymol", "-cq"]
    pymol.finish_launching()
    
    from pymol import cmd
    
    cmd.reinitialize()
    
    # Load structures
    cmd.load(cif_a, "target")
    cmd.load(cif_b, "orphan")
    
    console.print(f"\n[bold]Structural Validation: {uniprot_a} vs {uniprot_b}[/bold]")
    
    # Try cealign first since sequence identity is low (~20%)
    try:
        result = cmd.cealign("target", "orphan")
        rmsd = result['RMSD']
        aligned_atoms = result['alignment_length']
        console.print(f"Method: cealign (low sequence identity)")
        console.print(f"RMSD: [bold green]{rmsd:.3f} Å[/bold green] over {aligned_atoms} atoms")
    except Exception as e:
        console.print(f"[red]cealign failed: {e}[/red]")
        cmd.quit()
        return

    # Render configuration
    cmd.show("cartoon", "all")
    cmd.color("cyan", "target")
    cmd.color("salmon", "orphan")
    
    # Focus on the aligned region
    cmd.orient()
    
    png_path = f"output/superposition_{uniprot_a}_{uniprot_b}.png"
    pse_path = f"output/session_{uniprot_a}_{uniprot_b}.pse"
    
    console.print("Rendering high-quality image...")
    cmd.set("ray_opaque_background", 1)
    cmd.png(png_path, width=1200, height=900, dpi=150)
    cmd.save(pse_path)
    
    console.print(f"[green]Saved render to {png_path}[/green]")
    console.print(f"[green]Saved PyMOL session to {pse_path}[/green]")
    
    cmd.quit()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", nargs='+', default=["P97499:A0A8G1A137", "A0A0R4IZ84:A0A9X9T8I6"])
    args = parser.parse_args()
    
    for pair in args.pairs:
        a, b = pair.split(":")
        align_structures(a, b)
        print("-" * 50)
