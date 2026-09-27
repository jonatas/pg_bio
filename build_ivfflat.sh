#!/bin/bash
echo "Building IVFFLAT index..." > db_optimize.log
psql "postgresql://localhost:28818/bio_demo" -c "CREATE INDEX idx_protein_embedding ON proteins USING ivfflat (embedding halfvec_cosine_ops) WITH (lists = 1500);" >> db_optimize.log 2>&1
echo "Index built! Restarting overnight mining..." >> db_optimize.log
nohup ./run_overnight_mining.sh &
