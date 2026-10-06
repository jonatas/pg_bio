# Continuous Night Discovery 🌙

This example script acts as an autonomous background daemon that searches for uncharacterized (orphan) proteins in the UniProt database, folds them, and automatically publishes a blog post about its findings.

It is designed to run unattended overnight, constantly discovering new "dark proteome" proteins by looping through a predefined list of biological bait keywords (e.g., "Amidase", "Luciferase").

## How it Works
1. **Search:** It uses `pg_bio` to search UniProt for uncharacterized proteins matching a bait keyword.
2. **Embed:** It generates a high-dimensional vector embedding for the sequence.
3. **Compare:** It queries your local PostgreSQL database to find the closest known structurally similar proteins using pgvector cosine distance.
4. **Publish:** If the similarity passes a strict threshold (distance <= 0.06), it writes a Jekyll blog post summarizing the structural homology and fetches/folds the 3D models using AlphaFold or ESMFold.

## Usage

Make sure you have your virtual environment activated and `uv` installed.

```bash
uv run python continuous_night_discovery.py
```

The script will run infinitely, pausing periodically between discoveries to avoid rate limits.
