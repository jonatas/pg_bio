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

@dataclass
class RasaProfile:
    residue_count: int
    total_sasa: float
    polar_sasa: float
    apolar_sasa: float
    mean_rasa: float
    core_buried_fraction: float
    surface_hydrophobic_fraction: float
    per_residue_rasa: List[float]

@dataclass
class HBondAudit:
    total_hbonds: int
    backbone_hbonds: int
    sidechain_hbonds: int
    cross_interface_hbonds: int
    buried_unsatisfied_polars: int
    hbond_satisfaction_rate: float

@dataclass
class ChaperoneEvaluation:
    ca_count: int
    sc_rmsd: float
    mean_plddt: float
    mean_rasa: float
    core_buried_fraction: float
    surface_hydrophobic_fraction: float
    interface_delta_sasa: float
    total_hbonds: int
    cross_interface_hbonds: int
    buried_unsatisfied_polars: int
    radius_of_gyration: float
    is_globular: bool
    design_composite_score: float
    passes_chaperone_gate: bool

@dataclass
class ChaperoneCoreExampleAudit:
    example_name: str
    modality: str
    sc_rmsd: float
    mean_plddt: float
    interface_delta_sasa: float
    cross_interface_hbonds: int
    buried_unsatisfied_polars: int
    active_site_or_motif_rmsd: float
    ligand_or_cofactor_verified: bool
    composite_score: float
    passes_example_gate: bool


@dataclass
class LigandProperties:
    mw: float
    logp: float
    hbd: float
    hba: float
    tpsa: float
    rotatable_bonds: float
    formal_charge: float
    aromatic_rings: float
    fsp3: float


@dataclass
class RangeFinderPocket:
    pocket_id: str
    atom_count: int
    edge_count_5a: int
    mw_mu: float
    mw_sigma: float
    mw_min: float
    mw_max: float
    logp_mu: float
    logp_sigma: float
    logp_min: float
    logp_max: float
    hbd_mu: float
    hbd_sigma: float
    hbd_min: float
    hbd_max: float
    hba_mu: float
    hba_sigma: float
    hba_min: float
    hba_max: float
    tpsa_mu: float
    tpsa_sigma: float
    tpsa_min: float
    tpsa_max: float
    rot_mu: float
    rot_sigma: float
    rot_min: float
    rot_max: float
    charge_mu: float
    charge_sigma: float
    charge_min: float
    charge_max: float
    arom_mu: float
    arom_sigma: float
    arom_min: float
    arom_max: float
    fsp3_mu: float
    fsp3_sigma: float
    fsp3_min: float
    fsp3_max: float
    aromatic_logit: float
    aromatic_probability: float


@dataclass
class RangeFinderEvaluation:
    pocket_id: str
    include_formal_charge: bool
    range_iou: float
    giou_loss: float
    library_captured_fraction: float
    reduction_factor: float
    discard_fraction: float
    recall: float
    rrf_score: float





