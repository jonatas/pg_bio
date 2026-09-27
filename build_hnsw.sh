#!/bin/bash
psql "postgresql://localhost:28818/bio_demo" -c "SET maintenance_work_mem = '4GB'; CREATE INDEX idx_protein_embedding ON proteins USING hnsw (embedding halfvec_cosine_ops);" > hnsw_build.log 2>&1
echo "HNSW Build complete!" >> hnsw_build.log
