from dataclasses import dataclass
import pandas as pd
import numpy as np

@dataclass
class BiologicalDecompositionResult:
    gene_results: pd.DataFrame

@dataclass
class BiologicalDecompositionRunResult:
    result: BiologicalDecompositionResult
    previews: dict
    warnings: list

class BiologicalDecompositionEngine:
    def run(self, bulk_matrix, deconv_result, rhythm_results=None):
        if bulk_matrix is None or deconv_result is None or deconv_result.cell_fractions.empty:
            return BiologicalDecompositionRunResult(BiologicalDecompositionResult(pd.DataFrame()), {}, ["Bulk matrix or deconvolution result missing."])
        frac = deconv_result.cell_fractions.copy()
        if "sample_id" not in frac.columns:
            return BiologicalDecompositionRunResult(BiologicalDecompositionResult(pd.DataFrame()), {}, ["Deconvolution fractions missing sample_id."])
        frac = frac.set_index("sample_id")
        common_samples = [s for s in bulk_matrix.columns if s in frac.index]
        if len(common_samples) < 3:
            return BiologicalDecompositionRunResult(BiologicalDecompositionResult(pd.DataFrame()), {}, ["Need at least 3 matched samples for decomposition."])
        frac = frac.loc[common_samples]
        results = []
        for gene_id in bulk_matrix.index.astype(str):
            y = bulk_matrix.loc[gene_id, common_samples].astype(float).values
            xvar = frac.var(axis=1).mean()
            yvar = float(np.var(y))
            comp_score = float(min(1.0, xvar / (yvar + 1e-8))) if yvar > 0 else 0.0
            residual_score = float(max(0.0, 1.0 - comp_score))
            label = "composition-mediated" if comp_score >= 0.5 else "residual/intrinsic"
            rhythm_support = None
            if rhythm_results is not None and not rhythm_results.empty and "gene_id" in rhythm_results.columns:
                sub = rhythm_results[rhythm_results["gene_id"].astype(str) == str(gene_id)]
                if not sub.empty and "fit_r2" in sub.columns:
                    rhythm_support = float(sub.iloc[0]["fit_r2"])
            if rhythm_support is not None and rhythm_support >= 0.3 and comp_score < 0.5:
                label = "intrinsic-rhythmic"
            results.append({"gene_id": gene_id, "composition_score": comp_score, "residual_score": residual_score, "rhythm_support_r2": rhythm_support, "decomposition_label": label})
        df = pd.DataFrame(results)
        return BiologicalDecompositionRunResult(BiologicalDecompositionResult(df), {"decomposition_preview": df.head(50)}, [])
