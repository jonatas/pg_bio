from dataclasses import dataclass
from typing import List, Optional, Tuple

@dataclass
class Protein:
    uniprot_id: str
    name: str
    sequence: str
    embedding_distance: Optional[float] = None

@dataclass
class Atom:
    atom_id: int
    uniprot_id: str
    x: float
    y: float
    z: float
    name: str

@dataclass
class SpatialBounds:
    min_x: float
    max_x: float
    min_y: float
    max_y: float
    min_z: float
    max_z: float
