import urllib.parse
import urllib.request
import json
import psycopg

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

# 1. Canonical benchmark enzymes from AmazonZyme's published papers & theses
BENCHMARK_ACCESSIONS = [
    ("B1B533", "Penicillium citrinum (LMI01 / UFAM2 35kDa Xylanase)"),
    ("P69327", "Aspergillus awamori (SSF Glucoamylase GA1/GA2)"),
    ("P56271", "Aspergillus niger (SSF Acid Alpha-Amylase)"),
    ("P55329", "Aspergillus niger (UFAM1 Endo-1,4-beta-xylanase A)"),
    ("P48825", "Aspergillus aculeatus (UFAM1/2 Beta-glucosidase 1)"),
    ("P07987", "Trichoderma reesei (Industrial CBH II Exoglucanase 2)"),
    ("P07981", "Trichoderma reesei (Industrial EG-1 Endoglucanase)"),
    ("P36217", "Trichoderma reesei (Industrial Endo-1,4-beta-xylanase 2)"),
    ("P06278", "Bacillus licheniformis (Thermostable Liquefaction Alpha-Amylase)"),
    ("P00692", "Bacillus amyloliquefaciens (Liquefaction Alpha-Amylase)"),
    ("P18429", "Bacillus subtilis (Endo-1,4-beta-xylanase A)"),
    ("Q0KD98", "Cupriavidus necator (TPI FosLip2 93%-id Lipase)"),
    ("P23608", "Cupriavidus necator (TPI Bioplastic PHA Synthase PhaC)"),
    ("A0A8E4UR43", "Uncultured Bacteroidetes (TPI Clone Fos4 / Ade1 Carboxylesterase, PDB 8V16)"),
    ("A0ABW4BNN5", "Uncultured Bacterium (TPI Clone FosAmy5 / Max104 GH57 Alpha-Amylase)"),
]

# 2. Targeted UniProt queries across AmazonZyme's most-used organisms (uncharacterized + family homologs)
ORGANISM_QUERIES = [
    # Penicillium citrinum & Penicillium chrysogenum/rubens (UFAM2 LMI01 fungal workhorses)
    ("Penicillium citrinum", '(organism_name:"Penicillium citrinum" OR organism_name:"Penicillium") AND (protein_name:uncharacterized OR protein_name:hydrolase OR protein_name:xylanase OR protein_name:cellulase OR protein_name:glucosidase OR protein_name:esterase) AND length:[200 TO 650]'),
    # Aspergillus niger / awamori / oryzae / fumigatus (UFAM1 LMI03 & SSF Corn Ethanol workhorses)
    ("Aspergillus spp.", '(organism_name:"Aspergillus niger" OR organism_name:"Aspergillus awamori" OR organism_name:"Aspergillus oryzae" OR organism_name:"Aspergillus terreus") AND (protein_name:uncharacterized OR protein_name:glycoside OR protein_name:amylase OR protein_name:glucoamylase OR protein_name:xylanase OR protein_name:cellulase) AND length:[200 TO 680]'),
    # Trichoderma reesei & Trichoderma harzianum (Celluclast / CTec2 benchmark & Amazonian soil isolates)
    ("Trichoderma spp.", '(organism_name:"Trichoderma reesei" OR organism_name:"Trichoderma harzianum" OR organism_name:"Trichoderma virens") AND (protein_name:uncharacterized OR protein_name:glycoside OR protein_name:cellulase OR protein_name:xylanase OR protein_name:glucanase) AND length:[200 TO 600]'),
    # Bacillus licheniformis / subtilis / amyloliquefaciens (Thermostable liquefaction & soil isolates)
    ("Bacillus spp.", '(organism_name:"Bacillus licheniformis" OR organism_name:"Bacillus subtilis" OR organism_name:"Bacillus amyloliquefaciens" OR organism_name:"Bacillus") AND (protein_name:uncharacterized OR protein_name:amylase OR protein_name:xylanase OR protein_name:esterase OR protein_name:lipase) AND length:[200 TO 600]'),
    # Cupriavidus necator / Ralstonia (FosLip2 lipase & PHA synthase bioplastics)
    ("Cupriavidus necator", '(organism_name:"Cupriavidus necator" OR organism_name:"Cupriavidus") AND (protein_name:lipase OR protein_name:esterase OR protein_name:"alpha/beta hydrolase" OR protein_name:"poly(3-hydroxyalkanoate)") AND length:[220 TO 620]'),
    # Soil Bacteroidota & Acidobacteriota (Terra Preta de Indio metagenomic phyla for Ade1 / Fos4)
    ("Amazonian Soil Bacteroidota/Acidobacteriota", '(taxonomy_name:"Bacteroidota" OR taxonomy_name:"Acidobacteriota") AND (protein_name:carboxylesterase OR protein_name:esterase OR protein_name:"alpha/beta hydrolase") AND length:[240 TO 320]'),
]


def fetch_uniprot_batch(query: str, size: int = 80):
    params = urllib.parse.urlencode({
        "query": query,
        "format": "json",
        "size": str(size),
    })
    url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "pg_bio/1.0"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    out = []
    for entry in data.get("results", []):
        acc = entry.get("primaryAccession", "")
        org = entry.get("organism", {}).get("scientificName", "Unknown")
        tax_id = entry.get("organism", {}).get("taxonId", 0)
        seq = entry.get("sequence", {}).get("value", "")
        pname = "Uncharacterized protein"
        desc = entry.get("proteinDescription", {})
        if "recommendedName" in desc:
            pname = desc["recommendedName"].get("fullName", {}).get("value", pname)
        elif "submissionNames" in desc and desc["submissionNames"]:
            pname = desc["submissionNames"][0].get("fullName", {}).get("value", pname)
        if acc and seq:
            full_header = f"{pname} OS={org} OX={tax_id}"
            out.append((acc, full_header, org, seq))
    return out


def main():
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            print("=== 1. Ingesting & Indexing AmazonZyme's 15 Canonical Benchmark Enzymes ===")
            for acc, label in BENCHMARK_ACCESSIONS:
                try:
                    rows = fetch_uniprot_batch(f"accession:{acc}", size=1)
                    if rows:
                        uid, header, org, seq = rows[0]
                        cur.execute("""
                            INSERT INTO proteins (uniprot_id, name, sequence, embedding)
                            VALUES (%s, %s, %s, get_esm_embedding(%s)::halfvec(1280))
                            ON CONFLICT (uniprot_id) DO UPDATE
                            SET name = EXCLUDED.name,
                                sequence = EXCLUDED.sequence,
                                embedding = get_esm_embedding(EXCLUDED.sequence)::halfvec(1280);
                        """, (uid, f"{header} [{label}] [AmazonZyme-Organism]", seq, seq))
                        print(f"  + Indexed {uid:10s} | {label} ({len(seq)} aa)")
                except Exception as e:
                    print(f"  ! Failed {acc}: {e}")
            conn.commit()

            print("\n=== 2. Streaming Candidate Enzymes from AmazonZyme's Core Organisms ===")
            total_inserted = 0
            for group_name, q in ORGANISM_QUERIES:
                try:
                    entries = fetch_uniprot_batch(q, size=100)
                    group_count = 0
                    for uid, header, org, seq in entries:
                        tag_header = f"{header} [AmazonZyme-Organism: {group_name}]"
                        cur.execute("""
                            INSERT INTO proteins (uniprot_id, name, sequence, embedding)
                            VALUES (%s, %s, %s, get_esm_embedding(%s)::halfvec(1280))
                            ON CONFLICT (uniprot_id) DO UPDATE
                            SET name = EXCLUDED.name,
                                embedding = get_esm_embedding(EXCLUDED.sequence)::halfvec(1280);
                        """, (uid, tag_header, seq, seq))
                        group_count += 1
                        total_inserted += 1
                    conn.commit()
                    print(f"  + {group_name:42s}: {group_count:3d} proteins indexed with get_esm_embedding()")
                except Exception as e:
                    print(f"  ! Error querying {group_name}: {e}")

            print(f"\nTotal AmazonZyme organism proteins indexed/updated: {total_inserted}")

if __name__ == "__main__":
    main()
