-- ============================================================================
-- pg_bio Layer 15: RangeFinder (arXiv:2611.03379v1) Pocket-Interval GiST &
-- Composite B-Tree Virtual Screening Example
-- ============================================================================

-- 1. Create the Proteome Binding Pockets table with 8-column GiST index over numrange
DROP TABLE IF EXISTS rf_example_pockets CASCADE;
CREATE TABLE rf_example_pockets (
    pocket_id            TEXT PRIMARY KEY,
    uniprot_id           TEXT NOT NULL,
    target_name          TEXT NOT NULL,
    rf_profile           RangeFinderPocket NOT NULL,
    mw_range             NUMRANGE NOT NULL,
    logp_range           NUMRANGE NOT NULL,
    tpsa_range           NUMRANGE NOT NULL,
    hba_range            NUMRANGE NOT NULL,
    hbd_range            NUMRANGE NOT NULL,
    rot_range            NUMRANGE NOT NULL,
    arom_range           NUMRANGE NOT NULL,
    fsp3_range           NUMRANGE NOT NULL,
    charge_range         NUMRANGE NOT NULL,
    aromatic_probability FLOAT8 NOT NULL
);

CREATE INDEX idx_rf_example_pockets_gist
ON rf_example_pockets
USING gist (
    mw_range, logp_range, tpsa_range, hba_range,
    hbd_range, rot_range, arom_range, fsp3_range
);

-- 2. Create the Ligand Library table with Covering Composite B-Tree index
DROP TABLE IF EXISTS rf_example_ligands CASCADE;
CREATE TABLE rf_example_ligands (
    mol_id          BIGSERIAL PRIMARY KEY,
    compound_name   TEXT NOT NULL,
    smiles          TEXT NOT NULL,
    mw              FLOAT8 NOT NULL,
    logp            FLOAT8 NOT NULL,
    hbd             FLOAT8 NOT NULL,
    hba            FLOAT8 NOT NULL,
    tpsa            FLOAT8 NOT NULL,
    rotatable_bonds FLOAT8 NOT NULL,
    formal_charge   FLOAT8 NOT NULL,
    aromatic_rings  FLOAT8 NOT NULL,
    fsp3            FLOAT8 NOT NULL
);

CREATE INDEX idx_rf_example_ligands_btree
ON rf_example_ligands (mw, logp, tpsa, hba)
INCLUDE (mol_id, compound_name, smiles, hbd, rotatable_bonds, aromatic_rings, fsp3);

-- 3. Populate reference drug/bioprocess ligands using pure-Rust compute_ligand_properties_row()
INSERT INTO rf_example_ligands (
    compound_name, smiles, mw, logp, hbd, hba, tpsa,
    rotatable_bonds, formal_charge, aromatic_rings, fsp3
)
SELECT
    v.name, v.smiles,
    p.mw, p.logp, p.hbd, p.hba, p.tpsa,
    p.rotatable_bonds, p.formal_charge, p.aromatic_rings, p.fsp3
FROM (VALUES
    ('Aspirin',        'CC(=O)Oc1ccccc1C(=O)O'),
    ('Ibuprofen',      'CC(C)Cc1ccc(C(C)C(=O)O)cc1'),
    ('Naproxen',       'COc1ccc2cc(C(C)C(=O)O)ccc2c1'),
    ('Diclofenac',     'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl'),
    ('Vanillin',       'COc1cc(C=O)ccc1O'),
    ('Ferulic Acid',   'COc1cc(C=CC(=O)O)ccc1O'),
    ('2,5-FDCA',       'O=C(O)c1ccc(C(=O)O)o1'),
    ('Bisabolene',     'CC1=CCC(C(=CCC=C(C)C)C)CC1'),
    ('Cyclohexane',    'C1CCCCC1'),
    ('Stearic Acid',   'CCCCCCCCCCCCCCCCCC(=O)O')
) AS v(name, smiles)
CROSS JOIN LATERAL compute_ligand_properties_row(v.smiles) p;

-- 4. Populate empirical & predicted RangeFinder pockets
INSERT INTO rf_example_pockets (
    pocket_id, uniprot_id, target_name, rf_profile,
    mw_range, logp_range, tpsa_range, hba_range, hbd_range,
    rot_range, arom_range, fsp3_range, charge_range, aromatic_probability
)
SELECT
    nr.pocket_id, src.uniprot_id, src.target_name, src.profile,
    nr.mw_range, nr.logp_range, nr.tpsa_range, nr.hba_range,
    nr.hbd_range, nr.rot_range, nr.arom_range, nr.fsp3_range,
    nr.charge_range, nr.aromatic_probability
FROM (
    SELECT
        'COX2_NSAID_POCKET' AS pocket_id,
        'P35354' AS uniprot_id,
        'Prostaglandin G/H synthase 2 (COX-2)' AS target_name,
        empirical_pocket_ranges(
            'COX2_NSAID_POCKET',
            ARRAY[
                'CC(=O)Oc1ccccc1C(=O)O',
                'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
                'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl'
            ]
        ) AS profile
    UNION ALL
    SELECT
        'PHENOLIC_UPCYCLING_POCKET',
        'Q0CVS2',
        'Feruloyl Esterase / Vanillin Synthase Pocket',
        empirical_pocket_ranges(
            'PHENOLIC_UPCYCLING_POCKET',
            ARRAY[
                'COc1cc(C=O)ccc1O',
                'COc1cc(C=CC(=O)O)ccc1O',
                'O=C(O)c1ccc(C(=O)O)o1'
            ]
        )
) src
CROSS JOIN LATERAL pocket_numranges(src.profile) nr;

-- 5. Forward Screen: Find all ligands compatible with COX2_NSAID_POCKET via Composite B-Tree
SELECT l.compound_name, ROUND(l.mw::numeric, 2) AS mw, ROUND(l.logp::numeric, 2) AS logp, ROUND(l.tpsa::numeric, 1) AS tpsa
FROM pocket_interval_bounds(
    (SELECT rf_profile FROM rf_example_pockets WHERE pocket_id = 'COX2_NSAID_POCKET')
) b
JOIN rf_example_ligands l
  ON l.mw              BETWEEN b.mw_min   AND b.mw_max
 AND l.logp            BETWEEN b.logp_min AND b.logp_max
 AND l.tpsa            BETWEEN b.tpsa_min AND b.tpsa_max
 AND l.hba             BETWEEN b.hba_min  AND b.hba_max
 AND l.hbd             BETWEEN b.hbd_min  AND b.hbd_max
 AND l.rotatable_bonds BETWEEN b.rot_min  AND b.rot_max
 AND l.aromatic_rings  BETWEEN b.arom_min AND b.arom_max
 AND l.fsp3            BETWEEN b.fsp3_min AND b.fsp3_max
ORDER BY l.mw;

-- 6. Reverse Target Fishing: Given Ibuprofen SMILES, find compatible pockets via Multi-Column GiST `@>`
WITH lig AS (
    SELECT * FROM compute_ligand_properties_row('CC(C)Cc1ccc(C(C)C(=O)O)cc1')
)
SELECT p.pocket_id, p.uniprot_id, p.target_name, ROUND(p.aromatic_probability::numeric, 3) AS arom_prob
FROM rf_example_pockets p, lig
WHERE p.mw_range   @> lig.mw::numeric
  AND p.logp_range @> lig.logp::numeric
  AND p.tpsa_range @> lig.tpsa::numeric
  AND p.hba_range  @> lig.hba::numeric
  AND p.hbd_range  @> lig.hbd::numeric
  AND p.rot_range  @> lig.rotatable_bonds::numeric
  AND p.arom_range @> lig.aromatic_rings::numeric
  AND p.fsp3_range @> lig.fsp3::numeric;
