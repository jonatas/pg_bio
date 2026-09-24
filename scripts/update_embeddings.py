# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
# ]
# ///
import psycopg
import time
from pgbio import PgBioClient

print("Waiting for cargo pgrx install to finish...")
# We will wait until the task completes (just a blind wait or polling)
