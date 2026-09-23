# pg_bio Backend Engine

A high-performance FastAPI backend that serves as the API layer over the `pg_bio` PostgreSQL database. It seamlessly integrates the `pgbio-py` Python SDK to expose advanced computational biology workflows.

## Setup

1. Install dependencies (we recommend using `uv`):
```bash
uv venv
source .venv/bin/activate
uv pip install -e ../pgbio-py
uv pip install -r requirements.txt
```

2. Run the server:
```bash
uvicorn main:app --reload --port 8000
```

## Features

* **`POST /api/homologues`**: Find structural homologues using K-mer vector embeddings.
* **`GET /api/proteins/{id}/attention/{idx}`**: Traverse Sparse Attention Maps natively in Postgres.
* **`GET /api/spatial/radius`**: Ultra-fast Z-Order indexed 3D spatial queries.
* **`POST /api/research/workflow`**: A massive aggregated endpoint that combines all three features to demonstrate a full research pipeline: predicting homologues, finding the active site via AI attention, and fetching the 3D binding pocket coordinates!
