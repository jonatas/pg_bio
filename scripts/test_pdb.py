import requests

def parse_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url)
    if resp.status_code != 200:
        print(f"Failed to fetch {pdb_id}")
        return
    
    seqres = []
    atoms = []
    
    for line in resp.text.splitlines():
        if line.startswith("SEQRES"):
            # SEQRES   1 A  310  PRO ASN VAL VAL GLU LYS
            parts = line[19:].split()
            seqres.extend(parts)
        elif line.startswith("ATOM  "):
            # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
            atom_name = line[12:16].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            atoms.append({"name": atom_name, "x": x, "y": y, "z": z})
            
    print(f"Seq length (tokens): {len(seqres)}")
    print(f"Atoms: {len(atoms)}")
    if atoms:
        print(f"First atom: {atoms[0]}")

parse_pdb("1CRN")
