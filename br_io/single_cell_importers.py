from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Any
import pandas as pd
from bayesrhythm_core_models import ImportAuditTrail, SingleCellReferenceInput
from br_io.mex_loader import load_10x_mex_directory, load_10x_mex_archive
@dataclass
class SingleCellImportResult:
    single_cell: SingleCellReferenceInput
    preview: Dict[str,pd.DataFrame]
    annotation_required: bool=False
class SingleCellImporter:
    def _read_table(self, path):
        path=Path(path)
        if path.suffix.lower()==".xlsx": return pd.read_excel(path)
        return pd.read_csv(path, sep="," if path.suffix.lower()==".csv" else "\t")
    def import_bundle(self, expression_file, cell_metadata_file, gene_metadata_file=None, dataset_metadata: Optional[Dict[str,Any]]=None):
        expr=self._read_table(expression_file); meta=self._read_table(cell_metadata_file); genes=self._read_table(gene_metadata_file) if gene_metadata_file is not None else pd.DataFrame()
        if "cell_id" not in expr.columns: expr=expr.rename(columns={expr.columns[0]:"cell_id"})
        expr["cell_id"]=expr["cell_id"].astype(str); expr=expr.set_index("cell_id")
        if "cell_id" not in meta.columns: meta=meta.rename(columns={meta.columns[0]:"cell_id"})
        meta["cell_id"]=meta["cell_id"].astype(str)
        if "cell_type" not in meta.columns: meta["cell_type"]="unknown"
        sc=SingleCellReferenceInput(expr, meta, genes, dataset_metadata or {}, audit_trail=ImportAuditTrail(source_files=[str(expression_file),str(cell_metadata_file)]))
        preview={"cell_metadata_head":meta.head(20),"cell_type_counts":meta["cell_type"].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n_cells")}
        return SingleCellImportResult(sc, preview, annotation_required=meta["cell_type"].astype(str).str.lower().eq("unknown").all())
    def import_mex_directory(self, directory_path, dataset_metadata: Optional[Dict[str,Any]]=None):
        sc=load_10x_mex_directory(directory_path, dataset_metadata=dataset_metadata); meta=sc.cell_metadata
        preview={"cell_metadata_head":meta.head(20),"cell_type_counts":meta["cell_type"].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n_cells")}
        return SingleCellImportResult(sc, preview, annotation_required=True)
    def import_mex_archive(self, archive_path, dataset_metadata: Optional[Dict[str,Any]]=None):
        sc=load_10x_mex_archive(archive_path, dataset_metadata=dataset_metadata); meta=sc.cell_metadata
        preview={"cell_metadata_head":meta.head(20),"cell_type_counts":meta["cell_type"].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n_cells")}
        return SingleCellImportResult(sc, preview, annotation_required=True)
