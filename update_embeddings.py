import psycopg
import time
from pgbio import PgBioClient

print("Waiting for cargo pgrx install to finish...")
# We will wait until the task completes (just a blind wait or polling)
