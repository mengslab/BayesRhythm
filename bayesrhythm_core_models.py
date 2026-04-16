from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd
class AnalysisMode(str, Enum): SINGLE_CONDITION="single_condition"; COMPARATIVE="comparative"
class RhythmMode(str, Enum): CIRCADIAN="circadian"; ULTRADIAN="ultradian"; BOTH="both"
class ReferenceStrategy(str, Enum): CONDITION_MATCHED="condition_matched"; CHOW_ONLY="chow_only"; TRF_ONLY="trf_only"; CROSS_CONDITION_SENSITIVITY="cross_condition_sensitivity"
@dataclass
class DatasetConfig:
    study_name:str; species:str; analysis_mode:AnalysisMode; conditions:List[str]; expected_rhythm_modes:RhythmMode; period_search_range:Tuple[float,float]; reference_strategy:ReferenceStrategy; notes:str=""
@dataclass
class ImportAuditTrail:
    source_files:List[str]=field(default_factory=list); notes:List[str]=field(default_factory=list)
@dataclass
class SingleCellReferenceInput:
    expression_matrix:Any
    cell_metadata:pd.DataFrame
    gene_metadata:pd.DataFrame
    dataset_metadata:Dict[str,Any]=field(default_factory=dict)
    layers:Dict[str,Any]=field(default_factory=dict)
    audit_trail:Optional[ImportAuditTrail]=None
@dataclass
class AnchorDatasetRecord:
    dataset_id:str; condition:str; zeitgeber_time:float; source_file:str; n_cells:int; n_genes:int; annotation_status:str; validation_status:str; readiness_status:str; import_mode:str; dataset_label:str; low_confidence_fraction:float=0.0; unresolved_fraction:float=0.0; disagreement_fraction:float=0.0
@dataclass
class AnchorPanel:
    anchor_records:List[AnchorDatasetRecord]; datasets:Dict[str,SingleCellReferenceInput]; panel_metadata:Dict[str,Any]=field(default_factory=dict)
    def to_summary_table(self)->pd.DataFrame:
        return pd.DataFrame([{"dataset_id":r.dataset_id,"dataset_label":r.dataset_label,"condition":r.condition,"zeitgeber_time":r.zeitgeber_time,"n_cells":r.n_cells,"n_genes":r.n_genes,"annotation_status":r.annotation_status,"validation_status":r.validation_status,"readiness_status":r.readiness_status,"low_confidence_fraction":r.low_confidence_fraction,"unresolved_fraction":r.unresolved_fraction,"disagreement_fraction":r.disagreement_fraction,"import_mode":r.import_mode,"source_file":r.source_file} for r in self.anchor_records])
@dataclass
class AnchorPanelValidationResult:
    is_complete:bool; completeness_table:pd.DataFrame; missing_anchor_keys:List[str]; warnings:List[str]
