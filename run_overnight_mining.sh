#!/bin/bash
echo "Starting Overnight Dark Proteome Mining..." > overnight_discoveries.log
date >> overnight_discoveries.log

families=("Argonaute" "Telomerase" "Polymerase" "PETase" "Luciferase" "Helicase" "Integrase" "Isomerase")

for family in "${families[@]}"; do
    echo "========================================" >> overnight_discoveries.log
    echo "Mining Family: $family" >> overnight_discoveries.log
    echo "========================================" >> overnight_discoveries.log
    uv run scripts/batch_deorphanizer.py "$family" >> overnight_discoveries.log 2>&1
done

echo "Overnight Mining Complete!" >> overnight_discoveries.log
date >> overnight_discoveries.log
