from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy.optimize import nnls

@dataclass
class DeconvolutionConfig:
    min_shared_genes: int = 50

@dataclass
class DeconvolutionResult:
    cell_fractions: pd.DataFrame
    fit_metrics: pd.DataFrame

@dataclass
class DeconvolutionRunResult:
    result: DeconvolutionResult
    previews: dict
    warnings: list

class TrajectoryAwareDeconvolutionEngine:
    def run(self, bulk_matrix, bulk_metadata, atlas, config=None):
        cfg = config or DeconvolutionConfig()
        ref = atlas.reference_table.copy()
        if bulk_matrix is None or bulk_metadata is None or ref.empty:
            return DeconvolutionRunResult(DeconvolutionResult(pd.DataFrame(), pd.DataFrame()), {}, ["Bulk data or reference atlas missing."])
        gene_cols = [g for g in atlas.shared_genes if g in bulk_matrix.index and g in ref.columns]
        if len(gene_cols) < cfg.min_shared_genes:
            return DeconvolutionRunResult(DeconvolutionResult(pd.DataFrame(), pd.DataFrame()), {}, [f"Too few shared genes for deconvolution: {len(gene_cols)}"])
        X = ref[gene_cols].astype(float).T.values
        labels = (ref["anchor_id"].astype(str) + "|" + ref["cell_type"].astype(str)).tolist()
        frac_rows = []
        fit_rows = []
        for sample_id in bulk_matrix.columns:
            y = bulk_matrix.loc[gene_cols, sample_id].astype(float).values
            coef, _ = nnls(X, y)
            fracs = coef / coef.sum() if coef.sum() > 0 else coef
            recon = X @ fracs
            ss_res = float(np.sum((y - recon) ** 2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2)) if len(y) else 0.0
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
            row = {"sample_id": sample_id}
            for lab, val in zip(labels, fracs):
                row[lab] = float(val)
            frac_rows.append(row)
            fit_rows.append({"sample_id": sample_id, "fit_r2": r2, "n_shared_genes": len(gene_cols)})
        frac_df = pd.DataFrame(frac_rows)
        fit_df = pd.DataFrame(fit_rows)
        return DeconvolutionRunResult(DeconvolutionResult(frac_df, fit_df), {"fit_preview": fit_df.head(50), "fractions_preview": frac_df.head(20)}, [])
