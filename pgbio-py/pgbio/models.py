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

@dataclass
class ExpressionCassette:
    host_organism: str
    promoter_name: str
    coding_dna: str
    full_plasmid_dna: str
    length_bp: int
    gc_content: float
    cai_score: float
    restriction_sites_free: bool

@dataclass
class GlobularFoldingProfile:
    ca_count: int
    mean_plddt: float
    radius_of_gyration: float
    compactness_ratio: float
    elongation: float
    is_globular: bool



