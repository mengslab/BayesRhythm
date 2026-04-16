import pandas as pd, numpy as np
def score_markers_per_cell(sc_input, marker_library, use_log1p=True):
    expr=sc_input.expression_matrix
    rows=[]
    for label,markers in marker_library.items():
        present=[m for m in markers if m in expr.columns]
        if present:
            sub=expr[present]
            dense=sub.sparse.to_dense().astype(float)
            if use_log1p: dense=np.log1p(dense)
            score=dense.mean(axis=1)
            detected=(dense>0).sum(axis=1)
            frac=detected/max(1,len(present))
        else:
            score=pd.Series(0.0,index=expr.index); frac=pd.Series(0.0,index=expr.index)
        rows.append(pd.DataFrame({"cell_id":expr.index.astype(str),"label":label,"marker_support_score":score.values,"marker_fraction_detected":frac.values}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
def assign_marker_labels(marker_scores_df, min_support=0.1, ambiguity_margin=0.05):
    if marker_scores_df.empty: return pd.DataFrame()
    piv=marker_scores_df.pivot_table(index="cell_id", columns="label", values="marker_support_score", aggfunc="mean", fill_value=0.0); frac_piv=marker_scores_df.pivot_table(index="cell_id", columns="label", values="marker_fraction_detected", aggfunc="mean", fill_value=0.0); out=[]
    for cell_id,row in piv.iterrows():
        s=row.sort_values(ascending=False); top_label=s.index[0]; top=float(s.iloc[0]); second=float(s.iloc[1]) if len(s)>1 else 0.0; margin=top-second; frac=float(frac_piv.loc[cell_id, top_label]) if top_label in frac_piv.columns else 0.0; combined=0.65*top+0.35*frac; conf="high" if combined>=0.5 and margin>=ambiguity_margin else ("moderate" if combined>=0.2 else "low"); out.append({"cell_id":cell_id,"marker_predicted_label":top_label if combined>=0.2 else "unknown","marker_support_score":top,"marker_fraction_detected":frac,"marker_second_best_score":second,"marker_ambiguity_score":margin,"marker_combined_confidence_score":combined,"marker_confidence_label":conf})
    return pd.DataFrame(out)
