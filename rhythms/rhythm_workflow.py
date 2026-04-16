from dataclasses import dataclass
import pandas as pd

@dataclass
class RhythmWorkflowSummary:
    summary_table: pd.DataFrame

def summarize_rhythm_results(simple_df=None, advanced_df=None, comparison_df=None):
    rows = []
    if simple_df is not None and not simple_df.empty:
        rows.append({
            "engine": "rhythm_discovery",
            "n_genes": int(simple_df.shape[0]),
            "n_rhythmic": int(simple_df["rhythmic_call"].astype(str).eq("rhythmic").sum()) if "rhythmic_call" in simple_df.columns else 0,
        })
    if advanced_df is not None and not advanced_df.empty:
        rows.append({
            "engine": "advanced_rhythm",
            "n_genes": int(advanced_df.shape[0]),
            "n_rhythmic": int(advanced_df["rhythmic_call"].astype(str).eq("rhythmic").sum()) if "rhythmic_call" in advanced_df.columns else 0,
        })
    if comparison_df is not None and not comparison_df.empty:
        rows.append({
            "engine": "rhythm_comparison",
            "n_genes": int(comparison_df.shape[0]),
            "n_rhythmic": int(comparison_df["comparison_call"].astype(str).eq("changed").sum()) if "comparison_call" in comparison_df.columns else 0,
        })
    return RhythmWorkflowSummary(pd.DataFrame(rows))
