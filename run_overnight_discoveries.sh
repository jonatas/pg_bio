#!/bin/bash
echo "Starting overnight Dark Proteome discovery pipeline..."

families=("Kinase" "Polymerase" "Cutinase" "CRISPR" "Cas9" "Argonaute")

for family in "${families[@]}"; do
    echo "========================================"
    echo "Starting batch discovery for: $family"
    echo "========================================"
    uv run scripts/batch_deorphanizer.py "$family" >> "overnight_discovery_${family}.log" 2>&1
done

echo "Overnight pipeline complete!"
