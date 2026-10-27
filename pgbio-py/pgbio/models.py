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

@dataclass
class PreorganizationProfile:
    metric_name: str
    bin_edges: List[float]
    bin_counts: List[float]
    mean_val: float
    std_dev: float

@dataclass
class EnergyDecompositionAudit:
    model: str
    entity_count: int
    max_mm_polar_identity_error: float
    gas_closure_error: float
    polar_closure_error: float
    nonpolar_discrepancy: float
    total_mm_plus_polar: float
    gas_and_polar_closed: bool
    nonpolar_closed: bool


