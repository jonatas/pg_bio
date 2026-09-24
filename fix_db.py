import psycopg
DB_URI = "postgresql://localhost:28818/bio_demo"
with psycopg.connect(DB_URI, autocommit=True) as conn:
    with conn.cursor() as cur:
        cur.execute("DROP TABLE protein_atoms;")
        cur.execute("""
            CREATE TABLE protein_atoms (
                atom_id SERIAL PRIMARY KEY,
                uniprot_id VARCHAR(20) REFERENCES proteins(uniprot_id),
                coord ResidueCoord,
                z_index BIGINT GENERATED ALWAYS AS (residue_z_index(coord)) STORED
            );
            CREATE INDEX idx_spatial_z_order ON protein_atoms (z_index);
        """)
