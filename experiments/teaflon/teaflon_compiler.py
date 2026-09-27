import time
from rich.console import Console
from rich.layout import Layout
from rich.panel import Panel
from rich.align import Align
from rich.table import Table
from rich import box

console = Console()

def generate_dashboard():
    layout = Layout()
    layout.split_column(
        Layout(name="header", size=3),
        Layout(name="main")
    )
    layout["main"].split_row(
        Layout(name="left"),
        Layout(name="right")
    )
    layout["left"].split_column(
        Layout(name="genomics"),
        Layout(name="transcriptomics")
    )
    layout["right"].split_column(
        Layout(name="proteomics"),
        Layout(name="joint_space")
    )

    # Header
    layout["header"].update(Panel(Align.center("[bold cyan]🧬 pg_bio Multiomics Compiler : Teaflon (Syn-TFLN) Project[/]"), style="on black"))

    # Genomics
    g_table = Table(box=box.SIMPLE)
    g_table.add_column("Target", style="cyan")
    g_table.add_column("Status", style="green")
    g_table.add_row("Safe Harbor", "[17: 41492100-41496000]")
    g_table.add_row("Promoter", "KRT14 (Nails)")
    g_table.add_row("Tumor Suppressors", "0 Overlap (SAFE)")
    layout["genomics"].update(Panel(g_table, title="[bold magenta]Step 1: Genomics (GiST)[/]", border_style="magenta"))

    # Proteomics
    p_table = Table(box=box.SIMPLE)
    p_table.add_column("Atomic Site", style="cyan")
    p_table.add_column("Z-Order Bounding Box", style="green")
    p_table.add_row("ARG_BIND_508", "Matched (11.19, -2.97, 24.40)")
    p_table.add_row("ARG_BIND_506", "Matched (9.96, -3.68, 24.49)")
    layout["proteomics"].update(Panel(p_table, title="[bold blue]Step 2: Proteomics (Z-Order)[/]", border_style="blue"))

    # Transcriptomics
    t_table = Table(box=box.SIMPLE)
    t_table.add_column("Tissue", style="cyan")
    t_table.add_column("Dosage", style="green")
    t_table.add_row("Keratinocyte", "34.94 copies (Optimal)")
    t_table.add_row("Osteoblast", "0.51 copies (Safe)")
    layout["transcriptomics"].update(Panel(t_table, title="[bold yellow]Step 3: Transcriptomics (CSR)[/]", border_style="yellow"))

    # Joint Space
    j_table = Table(box=box.SIMPLE)
    j_table.add_column("Locus", style="cyan")
    j_table.add_column("Alignment (Cosine)", style="green")
    j_table.add_row("Chr17:Teaflon_Enhancer", "0.00000")
    layout["joint_space"].update(Panel(j_table, title="[bold red]Step 4: Joint Embeddings (vector)[/]", border_style="red"))

    return layout

if __name__ == "__main__":
    console.clear()
    with console.status("[bold green]Compiling synthetic multiomics genome..."):
        time.sleep(1)
    console.print(generate_dashboard())
