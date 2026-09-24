# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "mcp[cli]",
#     "torch",
#     "transformers",
# ]
# ///

from mcp.server.fastmcp import FastMCP
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
import os

# Initialize FastMCP Server
mcp = FastMCP("ESMFold Local Server")

# Global variables to cache the massive 11GB model in memory
_model = None
_tokenizer = None

def load_model():
    """Lazily loads the model only when a fold is requested."""
    global _model, _tokenizer
    if _model is None:
        device = "mps" if torch.backends.mps.is_available() else "cpu"
        _tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
        _model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).to(device)
    return _model, _tokenizer

@mcp.tool()
def fold_sequence(sequence: str, output_pdb_path: str) -> str:
    """
    Predicts the 3D atomic structure of a protein sequence using the local ESMFold Computational model.
    Saves the result to the specified output_pdb_path and returns a success message.
    """
    model, tokenizer = load_model()
    device = next(model.parameters()).device
    
    # Tokenize and run inference
    inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    
    with torch.no_grad():
        outputs = model(**inputs)
        
    # Convert tensor outputs to PDB format
    pdb_string = model.output_to_pdb(outputs)[0]
    
    # Save the file
    with open(output_pdb_path, "w") as f:
        f.write(pdb_string)
        
    return f"Success! Folded {len(sequence)} amino acids natively on {device}. Saved to {os.path.abspath(output_pdb_path)}."

if __name__ == "__main__":
    # Start the stdio transport for the agent to communicate with
    mcp.run(transport='stdio')
