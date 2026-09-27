import requests
import json
res = requests.get("https://alphafold.ebi.ac.uk/api/prediction/P97499")
print(res.status_code, res.text)
