import requests

query = {
  "query": {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
      "operator": "exact_match",
      "value": "Escherichia coli"
    }
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {
      "start": 0,
      "rows": 10
    }
  }
}

resp = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
print(resp.status_code)
print(resp.text[:200])
