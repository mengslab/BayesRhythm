from dataclasses import dataclass
import pandas as pd

MARKER_PRIORITY = ["Alb","Ttr","Apoa1","Adgre1","Csf1r","Lyz2","Kdr","Pecam1","Klf2","Krt19","Krt7","Sox9","Col1a1","Col1a2","Rgs5","Cd3d","Cd3e","Trbc2"]

@dataclass
class MultiAnchorReferenceBuildConfig:
    min_cells_per_cell_type: int = 20
    max_shared_genes: int = 3000

@dataclass
class ReferenceAtlas:
    shared_genes: list
    reference_table: pd.DataFrame
    support_table: pd.DataFrame

@dataclass
class MultiAnchorReferenceBuildResult:
    atlas: ReferenceAtlas
    previews: dict
    warnings: list

class MultiAnchorReferenceBuilder:
    def _select_shared_genes(self, anchor_manager, panel, max_shared_genes):
        gene_sets = []
        for rec in panel.anchor_records:
            expr = anchor_manager.get_dataset(rec.dataset_id).expression_matrix
            gene_sets.append(set(map(str, expr.columns)))
        shared = sorted(set.intersection(*gene_sets)) if gene_sets else []
        if len(shared) <= max_shared_genes:
            return shared
        priority = [g for g in MARKER_PRIORITY if g in shared]
        remaining = [g for g in shared if g not in priority]
        return priority + remaining[:max(0, max_shared_genes - len(priority))]

    def run(self, anchor_manager, config=None):
        cfg = config or MultiAnchorReferenceBuildConfig()
        panel = anchor_manager.get_anchor_panel()
        if not panel.anchor_records:
            return MultiAnchorReferenceBuildResult(ReferenceAtlas([], pd.DataFrame(), pd.DataFrame()), {}, ["No anchors registered."])
        shared_genes = self._select_shared_genes(anchor_manager, panel, cfg.max_shared_genes)
        ref_rows=[]; support_rows=[]; warnings=[]
        if len(shared_genes) < 100: warnings.append(f"Reference atlas built with only {len(shared_genes)} shared genes.")
        elif len(shared_genes) == cfg.max_shared_genes: warnings.append(f"Shared genes capped at {cfg.max_shared_genes} for memory safety.")
        for rec in panel.anchor_records:
            ds = anchor_manager.get_dataset(rec.dataset_id)
            expr = ds.expression_matrix.loc[:, shared_genes]
            dense = expr.sparse.to_dense() if hasattr(expr, "sparse") else expr.copy()
            meta = ds.cell_metadata.copy()
            for cell_type, sub in meta.groupby("cell_type"):
                ids = [c for c in sub["cell_id"].astype(str).tolist() if c in dense.index]
                if len(ids) < cfg.min_cells_per_cell_type:
                    continue
                centroid = dense.loc[ids].mean(axis=0)
                row = centroid.to_dict()
                row.update({"anchor_id": rec.dataset_id,"condition": rec.condition,"zeitgeber_time": rec.zeitgeber_time,"cell_type": str(cell_type),"n_cells": len(ids)})
                ref_rows.append(row)
                support_rows.append({"anchor_id": rec.dataset_id,"condition": rec.condition,"zeitgeber_time": rec.zeitgeber_time,"cell_type": str(cell_type),"n_cells": len(ids)})
            del dense
        ref_df=pd.DataFrame(ref_rows); support_df=pd.DataFrame(support_rows)
        atlas=ReferenceAtlas(shared_genes, ref_df, support_df)
        return MultiAnchorReferenceBuildResult(atlas, {"support_table": support_df, "reference_preview": ref_df.head(50)}, warnings)
