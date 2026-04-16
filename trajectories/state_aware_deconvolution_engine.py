from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy.optimize import nnls

@dataclass
class StateAwareDeconvolutionResult:
    cell_fractions: pd.DataFrame
    fit_metrics: pd.DataFrame
    state_program_scores: pd.DataFrame

@dataclass
class StateAwareDeconvolutionRunResult:
    result: StateAwareDeconvolutionResult
    previews: dict
    warnings: list

class StateAwareDeconvolutionEngine:
    def run(self, bulk_matrix, bulk_metadata, atlas, min_shared_genes=100):
        ref = atlas.reference_table.copy()
        if bulk_matrix is None or bulk_metadata is None or ref.empty:
            return StateAwareDeconvolutionRunResult(StateAwareDeconvolutionResult(pd.DataFrame(), pd.DataFrame(), pd.DataFrame()), {}, ["Bulk data or reference atlas missing."])
        gene_cols = [g for g in atlas.shared_genes if g in bulk_matrix.index and g in ref.columns]
        if len(gene_cols) < min_shared_genes:
            return StateAwareDeconvolutionRunResult(StateAwareDeconvolutionResult(pd.DataFrame(), pd.DataFrame(), pd.DataFrame()), {}, [f"Too few shared genes: {len(gene_cols)}"])
        X = ref[gene_cols].astype(float).T.values
        labels = (ref["anchor_id"].astype(str) + "|" + ref["cell_type"].astype(str)).tolist()
        frac_rows=[]; fit_rows=[]; state_rows=[]
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
            fit_rows.append({"sample_id": sample_id, "fit_r2": r2, "n_shared_genes": len(gene_cols), "dominant_component": labels[int(np.argmax(fracs))] if len(fracs) else ""})
            # state-program scores collapsed by cell_type across anchors
            state_score = {"sample_id": sample_id}
            for cell_type in sorted(ref["cell_type"].astype(str).unique().tolist()):
                idx = [i for i, lab in enumerate(labels) if lab.endswith("|" + cell_type)]
                state_score[cell_type] = float(np.sum(fracs[idx])) if idx else 0.0
            state_rows.append(state_score)
        frac_df = pd.DataFrame(frac_rows)
        fit_df = pd.DataFrame(fit_rows)
        state_df = pd.DataFrame(state_rows)
        return StateAwareDeconvolutionRunResult(StateAwareDeconvolutionResult(frac_df, fit_df, state_df), {"fit_preview": fit_df.head(50), "fractions_preview": frac_df.head(20), "state_program_preview": state_df.head(20)}, [])
