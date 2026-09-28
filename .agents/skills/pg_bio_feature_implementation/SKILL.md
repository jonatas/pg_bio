---
name: pg_bio_feature_implementation
description: Guides the implementation of new features or GitHub issues for the pg_bio extension, ensuring rigorous testing, benchmarking, documentation, and end-to-end tutorials.
---

# Feature Implementation Guidelines for pg_bio

You have been invoked to implement a new feature (often from a GitHub issue) for the `pg_bio` PostgreSQL extension. Because `pg_bio` operates at the bleeding edge of database performance, systems engineering, and bioinformatics, you must strictly follow this 5-phase lifecycle for every feature you build.

## Phase 1: Rust Implementation (`pgrx`)
1. Implement the feature in `src/lib.rs` (or the appropriate Rust module).
2. Ensure the code aligns with the **"Compute to the Data"** philosophy: keep the heavy mathematical lifting (e.g., Morton coding, CSR matrix unrolling, Smith-Waterman) entirely in Rust to prevent Python I/O bottlenecks.

## Phase 2: Unit & Regression Testing
1. For every new Postgres type or `#[pg_extern]` function, you MUST write a corresponding unit test inside the `#[cfg(any(test, feature = "pg_test"))]` schema block.
2. Ensure you run the tests via `cargo pgrx test pg18` to verify SQL mappings and Rust logic.
3. Check for regressions. If modifying existing embeddings or spatial logic, verify it does not break the downstream pipelines (e.g., Teaflon).

## Phase 3: Performance Benchmarking
1. You must write a benchmark script in the `scripts/` directory to measure the new feature.
2. Compare the native `pg_bio` execution time against the "Default Python" approach (e.g., using `BioPython`, `NumPy`, or `SciPy`).
3. **Crucial Memory Lesson:** Always test limits. If the feature involves massive graph indexing or heavy memory structures (like HNSW), explicitly consider `maintenance_work_mem` ceilings. Ensure the feature handles SSD disk-spilling gracefully without causing PostgreSQL parallel worker deadlocks.

## Phase 4: Python SDK & Documentation
1. If the SQL interface changed, update the Python client (`pgbio-py/pgbio/client.py`) so scientists have a frictionless wrapper.
2. Update the README or the `docs/` folder to explain the biological advantage of the feature.
3. **Docs Style:** Write docs focused on *The Scientist's Advantage*. Explain the biological/chemical problem, how the native Rust/SQL math solves it, and provide a clear SQL code example.

## Phase 5: End-to-End Tutorial Testing
1. A feature is not complete until it is proven in an end-to-end workflow.
2. Create or update an interactive `.ipynb` Jupyter notebook in the `experiments/` directory showcasing a real-world scientific scenario for the feature.
3. **Execution Rule:** Do not just write the notebook and leave it blank. You MUST dynamically install dependencies (via `uv`) and execute the notebook in the background using `uv run jupyter nbconvert --to notebook --execute --inplace <notebook_name>.ipynb` to bake the live DataFrames and Matplotlib charts directly into the file, exactly as we did for the multiomics Teaflon project.

---
### Strict Rules to Remember
* **NEVER** use bash redirection (`cat << EOF > file`) to create or modify files. Use your native IDE tools (`write_to_file`, `replace_file_content`).
* Always verify that the database builds cleanly (`cargo pgrx run pg18`) before declaring the implementation complete.
