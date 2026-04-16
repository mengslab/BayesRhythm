from dataclasses import dataclass
import pandas as pd
import numpy as np

@dataclass
class StateProgramDecompositionResult:
    gene_results: pd.DataFrame

@dataclass
class StateProgramDecompositionRunResult:
    result: StateProgramDecompositionResult
    previews: dict
    warnings: list

class StateProgramDecompositionEngine:
    def run(self, bulk_matrix, state_aware_deconv_result, rhythm_results=None):
        if bulk_matrix is None or state_aware_deconv_result is None or state_aware_deconv_result.state_program_scores.empty:
            return StateProgramDecompositionRunResult(StateProgramDecompositionResult(pd.DataFrame()), {}, ["Bulk matrix or state-program scores missing."])
        state_df = state_aware_deconv_result.state_program_scores.copy()
        if "sample_id" not in state_df.columns:
            return StateProgramDecompositionRunResult(StateProgramDecompositionResult(pd.DataFrame()), {}, ["State-program table missing sample_id."])
        state_df = state_df.set_index("sample_id")
        common_samples = [s for s in bulk_matrix.columns if s in state_df.index]
        if len(common_samples) < 3:
            return StateProgramDecompositionRunResult(StateProgramDecompositionResult(pd.DataFrame()), {}, ["Need at least 3 matched samples."])
        X = state_df.loc[common_samples].astype(float)
        results = []
        for gene_id in bulk_matrix.index.astype(str):
            y = bulk_matrix.loc[gene_id, common_samples].astype(float)
            try:
                beta, *_ = np.linalg.lstsq(X.values, y.values, rcond=None)
                yhat = X.values @ beta
                ss_res = float(np.sum((y.values - yhat)**2))
                ss_tot = float(np.sum((y.values - np.mean(y.values))**2))
                r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 0.0
                rhythm_support = None
                if rhythm_results is not None and not rhythm_results.empty and "gene_id" in rhythm_results.columns:
                    sub = rhythm_results[rhythm_results["gene_id"].astype(str) == str(gene_id)]
                    if not sub.empty and "fit_r2" in sub.columns:
                        rhythm_support = float(sub.iloc[0]["fit_r2"])
                label = "state-program-mediated" if r2 >= 0.45 else "residual/intrinsic"
                if rhythm_support is not None and rhythm_support >= 0.35 and r2 < 0.45:
                    label = "intrinsic-rhythmic"
                results.append({"gene_id": gene_id, "state_program_r2": r2, "rhythm_support_r2": rhythm_support, "decomposition_label": label})
            except Exception:
                continue
        df = pd.DataFrame(results)
        return StateProgramDecompositionRunResult(StateProgramDecompositionResult(df), {"state_program_decomposition_preview": df.head(50)}, [])
