# Swarm Discovery & Unified Pipeline System UI (Migrated to [`jonatas/bio-pipelines`](https://github.com/jonatas/bio-pipelines))

To keep `pg_bio` focused strictly on the high-performance Rust `pgrx` PostgreSQL 18 extension (3D Morton Z-Order spatial indexing, ESM vector homology, Kabsch superposition, RangeFinder 10 Å pocket `numrange` GiST operators, and in-database RK4 bioprocess solvers), the **Autonomous Swarm Discovery Engine**, **Dark-Proteome (`DUF`) & Hypothesis Agent**, and **Unified Pipeline System UI & `/telemetry`** now live in [**`jonatas/bio-pipelines`** (`base/`)](https://github.com/jonatas/bio-pipelines/tree/main/base).

## Canonical Entrypoints in [`jonatas/bio-pipelines`](https://github.com/jonatas/bio-pipelines)

- **Autonomous Multi-Worker Swarm & In-Process Hypothesis Evolution (`every 25 waves`)**:
  ```bash
  git clone https://github.com/jonatas/bio-pipelines.git
  cd bio-pipelines
  uv run python base/swarm_discovery.py
  ```
- **Unified Pipeline System UI & `/telemetry` Server**:
  ```bash
  cd bio-pipelines
  uv run python base/swarm_web.py
  ```

### Unified Routes
- `/` & `/pipelines` — 8-Feedstock 4-Module ($M_1 \rightarrow M_2 \rightarrow M_3 \rightarrow M_4$) 48h RK4 Bioprocess Digital Twin
- `/discovery` — Multi-Worker Protein Correlation Swarm Heatmap
- `/telemetry` — Unified Pipeline System Telemetry, 5-Stage Closed-Loop Architecture, Live 8-Feedstock RangeFinder GiST Reverse Target Fishing Matrix, Autonomous Hypothesis Ledger (`bio_hypothesis_ledger`), and Evolving Dark-Proteome (`DUF`) Annotations (`bio_annotation_evolution`)
- `/cards` — 1-to-N Pokémon-Card Hub & Biophysical Similarity Arena
- `/compare` — Side-by-Side 3D Kabsch & RangeFinder Pocket Compare Lab
- `/galaxy` — 3D Organism Correlation Galaxy
- `/tpp` — CETSA / TPP Microprotein Knockout & Ideal Operating Temp ($T_{\text{opt}}$) Lab

*(Note: `swarm_discovery.py` and `swarm_web.py` in this directory are lightweight forwarding shims for backward compatibility when `BIO_PIPELINES_ROOT` or a sibling `../bio-pipelines` checkout is present).*
