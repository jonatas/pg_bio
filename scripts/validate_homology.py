# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
#     "biopython",
#     "rich",
# ]
# ///

import sys
import argparse
import requests
from Bio import Align
from Bio.Align import substitution_matrices
from rich.console import Console
from rich.table import Table

console = Console()

def fetch_sequence(uniprot_id: str) -> str:
    """Fetch protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.splitlines()
        seq = "".join(lines[1:])
        return seq
    
    # Try UniParc if not in UniProtKB
    url_parc = f"https://rest.uniprot.org/uniparc/{uniprot_id}.fasta"
    response_parc = requests.get(url_parc)
    if response_parc.status_code == 200:
        lines = response_parc.text.splitlines()
        seq = "".join(lines[1:])
        return seq
        
    console.print(f"[red]Failed to fetch sequence for {uniprot_id}[/red]")
    return ""

def validate_pair(bait_id: str, orphan_id: str):
    console.print(f"\n[bold]Validating {bait_id} vs {orphan_id}...[/bold]")
    seq1 = fetch_sequence(bait_id)
    seq2 = fetch_sequence(orphan_id)
    
    if not seq1 or not seq2:
        return
        
    console.print(f"Bait Length: {len(seq1)} aa | Orphan Length: {len(seq2)} aa")
    
    aligner = Align.PairwiseAligner()
    # Use BLOSUM62 for realistic protein alignment scoring
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    
    # Calculate identity
    matches = best_alignment.counts().identities
    # Use the shorter sequence length for coverage-based identity, or alignment length
    length = max(len(seq1), len(seq2)) 
    identity = (matches / length) * 100
    
    score = best_alignment.score
    console.print(f"Global Alignment Score: [bold cyan]{score}[/bold cyan]")
    console.print(f"Sequence Identity: [bold green]{identity:.2f}%[/bold green]")
    
    # Heuristics for Homology
    if identity < 20:
        console.print("[yellow]Warning: Sequence identity is < 20%. This is in the 'twilight zone' of sequence homology. The structural similarity could still be valid (distant homology), but pure sequence alignment doesn't strongly confirm it.[/yellow]")
    elif identity < 30:
        console.print("[blue]Notice: Sequence identity is 20-30%. Remote homology is possible. You might want to validate with HMMER or structural alignment (Foldseek).[/blue]")
    else:
        console.print("[green]Strong sequence homology confirmed! (>30% identity)[/green]")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate homology between two UniProt IDs")
    parser.add_argument("--pairs", nargs='+', help="Pairs of UniProt IDs in format BAIT:ORPHAN", default=[
        "Q99973:A0A060HRR0",
        "P97499:A0A060HRR0",
        "A0A0R4IZ84:A0A166BY34",
        "Q5RAK6:A0A202E3R8"
    ])
    args = parser.parse_args()
    
    console.print("[bold green]🧬 pg_bio Validation Script[/bold green]")
    console.print("Using Needleman-Wunsch global alignment with BLOSUM62 to cross-validate findings.\n")
    
    for pair in args.pairs:
        if ":" in pair:
            bait, orphan = pair.split(":", 1)
            validate_pair(bait, orphan)
        else:
            console.print(f"[red]Invalid pair format: {pair}. Use BAIT:ORPHAN[/red]")
