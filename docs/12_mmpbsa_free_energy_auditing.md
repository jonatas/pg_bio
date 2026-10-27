# 12. MM/GBSA & MM/PBSA Free-Energy Auditing & Allosteric Hotspot Discovery

> 🌐 **Interactive 3Dmol.js Blog Post:** [Auditing 500 ns Molecular Dynamics in SQL: Triterpenoids, 5-Lipoxygenase, and Cryptic Allosteric Pockets with pg_bio](https://ideia.me/5lox-triterpenoids-mmpbsa-free-energy-audit-pg-bio)  
> 📓 **Executed Jupyter Tutorial:** [`experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb`](../experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb) (validated on UFAM / COARACI [Zenodo 23565306](https://doi.org/10.5281/zenodo.23565306))

## The Scientist's Advantage

When medicinal chemists and computational biologists evaluate ligand binding across multi-hundred-nanosecond Molecular Dynamics (MD) ensembles using **MM/GBSA** (Generalized Born) and **MM/PBSA** (Poisson-Boltzmann) free-energy solvers, three recurring pitfalls obscure mechanism-of-action analysis:

1. **Unverified Endpoint Decomposition Sums (*The Double-Entry Accounting Trap*)**:
   Per-residue energy decomposition (`idecomp=2`) splits binding into `Internal`, `van der Waals`, `Electrostatic`, `Polar Solvation`, and `Non-Polar Solvation`. However, in standard `gmx_MMPBSA` / `MMPBSA.py` workflows:
   - In **MM/GBSA**, per-residue surface-area non-polar sums do not close exactly against the global LCPO `ESURF` term (`~0.8–1.1 kcal/mol` discrepancy).
   - In **MM/PBSA**, the per-residue `Non-Polar Solv.` column is identically `0.000 kcal/mol`, while the global PB endpoint includes a massive **`+20 to +26 kcal/mol`** cavity/dispersion shift (`ENPOLAR + EDISPER`). Summing raw per-residue `TOTAL` columns leads to erroneous claims of exact global-score decomposition.
2. **Pre-Defined Pocket Tunnel Vision (*Missing Cryptic Sub-Pockets*)**:
   Studies typically monitor a fixed list of pocket residues (e.g., 9 operational residues in human 5-lipoxygenase, `ALOX5`), missing strong peripheral stabilizing anchors (`Ile126`, `Val29`, `Leu66`) that engage when triterpenoids breathe or reorient over 500 ns.
3. **Mean-Value Blindness Over Long Trajectories**:
   Reporting only a 500 ns average hides whether a ligand stays deeply buried (`SASA > 4.2 nm²`) for 89% of the trajectory (`β-Amyrin`) or drifts out of the pocket (`34.9%` in `CID 101663382`).

`pg_bio` solves all three problems directly inside PostgreSQL 18 in compiled Rust.

---

## Native Rust Functions & Types

### 1. `audit_mmpbsa_decomposition(...) -> EnergyDecompositionAudit`
Performs a single-pass $O(N)$ double-entry thermodynamic audit across all $N$ entities (`676` protein residues + `1` ligand):
- Verifies per-entity identity: `MM_plus_polar = Internal + vdW + Electrostatic + Polar Solvation`.
- Verifies gas-phase closure: $|\sum (\text{Int} + \text{vdW} + \text{Elec}) - \text{GGAS}| \le \text{tolerance}$.
- Verifies polar solvation closure: $|\sum \text{Polar} - \text{E}_{\text{polar}}^{\text{global}}| \le \text{tolerance}$.
- Quantifies the exact `nonpolar_discrepancy` (`global_nonpolar_net - SUM(Non-Polar Solv.)`).

### 2. `solvent_shift_penalty(gb_polar, gb_nonpolar, pb_polar, pb_enpolar, pb_edisper) -> f64`
Computes the exact net solvation shift $\Delta \Delta G_{\text{solv}}^{\text{PB}-\text{GB}}$ reconciling `MM/GBSA` and `MM/PBSA` endpoints.

### 3. `build_energy_hotspot_map(...)` & `discover_cryptic_hotspots(...)`
Converts per-residue `MM_plus_polar` arrays into a Compressed Sparse Row (CSR) `SparseAttentionMap` and traverses it to return stabilizing residues outside the known operational pocket.

### 4. `allosteric_mechanism_index(...) -> f64`
Computes the normalized Polar/Ionic-to-Hydrophobic Routing Index in $[-1.0, +1.0]$, separating hydrophobic-clamped ligands (`< 0`, e.g., `Lupeol`, `β-Amyrin`) from H-bond/salt-bridge-switched ligands (`> 0`, e.g., `AKBA`).

---

## SQL Example (Zenodo 23565306: Human 5-LOX + 9 Triterpenoids)

```sql
SELECT
    t.ligand_label,
    UPPER(d.model) AS model,
    (audit_mmpbsa_decomposition(
        d.model,
        array_agg(d.internal_kcal),
        array_agg(d.vdw_kcal),
        array_agg(d.electrostatic_kcal),
        array_agg(d.polar_solv_kcal),
        array_agg(d.nonpolar_solv_kcal),
        array_agg(d.mm_plus_polar_kcal),
        g.ggas, g.global_polar, g.global_nonpolar, 0.005
    ))::text AS thermodynamic_audit,
    discover_cryptic_hotspots(
        build_energy_hotspot_map(
            677,
            array_agg(d.simulation_resid) FILTER (WHERE d.side = 'R'),
            array_agg(d.mm_plus_polar_kcal) FILTER (WHERE d.side = 'R'),
            999, -0.05
        ),
        999,
        ARRAY[68, 108, 109, 110, 129, 130, 133, 134, 137],
        0.45, 4
    ) AS cryptic_hotspots,
    round(allosteric_mechanism_index(
        build_energy_hotspot_map(
            677,
            array_agg(d.simulation_resid) FILTER (WHERE d.side = 'R'),
            array_agg(d.mm_plus_polar_kcal) FILTER (WHERE d.side = 'R'),
            999, -0.05
        ),
        999,
        ARRAY[109, 110, 130],
        ARRAY[68, 133, 134]
    )::numeric, 4) AS mechanism_index
FROM lox5_residue_decomposition d
JOIN lox5_triterpenoids t ON t.system_id = d.system_id
JOIN lox5_global_summary_view g ON g.system_id = d.system_id AND g.model = d.model
GROUP BY t.system_id, t.ligand_label, d.model, g.ggas, g.global_polar, g.global_nonpolar;
```
