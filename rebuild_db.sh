#!/bin/bash
echo "Starting database optimization..." > db_optimize.log
date >> db_optimize.log

echo "Converting to halfvec..." >> db_optimize.log
psql "postgresql://localhost:28818/bio_demo" -c "ALTER TABLE proteins ALTER COLUMN embedding TYPE halfvec(1280) USING embedding::halfvec;" >> db_optimize.log 2>&1

echo "Building HNSW index with 4GB work mem..." >> db_optimize.log
psql "postgresql://localhost:28818/bio_demo" -c "CREATE INDEX idx_protein_embedding ON proteins USING hnsw (embedding halfvec_cosine_ops);" >> db_optimize.log 2>&1

echo "Restarting overnight mining..." >> db_optimize.log
nohup ./run_overnight_mining.sh &

echo "Done!" >> db_optimize.log
date >> db_optimize.log
