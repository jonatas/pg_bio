# shellcheck shell=bash
# CREATE EXTENSION pg_bio CASCADE pulls in vector and pg_vectorsearch.
set -euo pipefail
: "${out:?set by nix}"
export PGDATA="$TMPDIR/pgdata" PGHOST="$TMPDIR/run"
mkdir -p "$PGHOST"
initdb -D "$PGDATA" --auth-local=trust --no-instructions >/dev/null
pg_ctl -D "$PGDATA" -w -o "-c listen_addresses='' -k $PGHOST" start
# Stop the server however the script ends, so a failure leaks nothing.
trap 'pg_ctl -D "$PGDATA" -m immediate stop >/dev/null 2>&1 || true' EXIT
# CASCADE fails, under ON_ERROR_STOP, if a required extension is missing.
# The two SELECTs run code from the libraries, not just their catalog rows.
psql -d postgres -v ON_ERROR_STOP=1 <<'SQL'
CREATE EXTENSION pg_bio CASCADE;
SELECT '[1,2]'::vector <-> '[1,3]'::vector;
SELECT count(*) FROM chaperone_core_examples();
SQL
touch "$out"
