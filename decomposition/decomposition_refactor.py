from dataclasses import dataclass
import pandas as pd

@dataclass
class DecompositionRefactorResult:
    unified_gene_table: pd.DataFrame

def build_unified_decomposition_table(basic_df=None, state_df=None):
    if (basic_df is None or basic_df.empty) and (state_df is None or state_df.empty):
        return DecompositionRefactorResult(pd.DataFrame())
    base = None
    if basic_df is not None and not basic_df.empty:
        base = basic_df.copy()
        base["source"] = "basic"
    if state_df is not None and not state_df.empty:
        s = state_df.copy()
        s["source"] = "state_program"
        if base is None:
            base = s
        else:
            common = [c for c in ["gene_id"] if c in base.columns and c in s.columns]
            if common:
                base = base.merge(s, on=common, how="outer", suffixes=("_basic", "_state"))
            else:
                base = pd.concat([base, s], ignore_index=True, sort=False)
    return DecompositionRefactorResult(base)
