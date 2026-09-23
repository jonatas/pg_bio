---
name: pg_bio_coding_standards
description: Enforces strict coding standards, IO usage, and testing methodologies for the pg_bio extension.
trigger: always_on
---

# pg_bio Extension Guidelines

## 1. File I/O Operations
- **CRITICAL:** NEVER use bash redirection (`cat << EOF > file`) to create or modify files. This can cause terminal hangs and is strictly forbidden by the user.
- **ALWAYS** use the native agent IDE tools (`write_to_file` and `replace_file_content`) for all file operations.

## 2. pgrx Extension Development
- All PostgreSQL extension code is written in Rust using the `pgrx` framework.
- For every new Postgres type or `#[pg_extern]` function created in `src/lib.rs`, you MUST write a corresponding unit test inside the `#[cfg(any(test, feature = "pg_test"))]` schema block.
- Tests must be executable via `cargo pgrx test pg18`.

## 3. Architecture Context
- `pg_bio` is designed for high-performance bioinformatics inside the database.
- It leverages:
  - **Z-Order (Morton Coding)** for 3D atomic spatial indexing.
  - **Vector Embeddings** for functional and structural homology (ESM models).
  - **Sparse Attention Maps** for Compressed Sparse Row (CSR) neural network weight storage.
