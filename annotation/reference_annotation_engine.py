import pandas as pd, numpy as np
def build_reference_centroids(reference_expr, reference_labels):
    ref=reference_expr.copy(); labels=pd.Series(reference_labels).astype(str); cols={}
    for label in labels.unique():
        ids=labels[labels==label].index.astype(str).tolist(); overlap=[i for i in ids if i in ref.index.astype(str)]
        if overlap: cols[label]=ref.loc[overlap].mean(axis=0)
    return pd.DataFrame(cols)
def score_cells_against_reference(sc_input, reference_centroids, max_shared_genes=5000):
    if reference_centroids.empty: return pd.DataFrame()
    expr=sc_input.expression_matrix; shared=[g for g in reference_centroids.index if g in expr.columns]
    if not shared: return pd.DataFrame()
    shared=shared[:max_shared_genes]
    dense=expr[shared].sparse.to_dense().astype(float); refs=reference_centroids.loc[shared].astype(float); rows=[]
    for cell_id,row in dense.iterrows():
        for label in refs.columns:
            if row.std()==0 or refs[label].std()==0: score=0.0
            else:
                score=float(np.corrcoef(row.values, refs[label].values)[0,1]); score=0.0 if np.isnan(score) else score
            rows.append({"cell_id":str(cell_id),"label":str(label),"reference_score":score})
    return pd.DataFrame(rows)
def assign_reference_labels(score_df, ambiguity_margin=0.05):
    if score_df.empty: return pd.DataFrame()
    piv=score_df.pivot_table(index="cell_id", columns="label", values="reference_score", aggfunc="mean", fill_value=0.0); out=[]
    for cell_id,row in piv.iterrows():
        s=row.sort_values(ascending=False); top_label=s.index[0]; top=float(s.iloc[0]); second=float(s.iloc[1]) if len(s)>1 else 0.0; margin=top-second; combined=max(0.0,min(1.0,0.5*(top+1.0)))*0.7+max(0.0,min(1.0,margin+0.5))*0.3; conf="high" if combined>=0.7 and margin>=ambiguity_margin else ("moderate" if combined>=0.45 else "low"); out.append({"cell_id":cell_id,"reference_predicted_label":top_label if top>0 else "unknown","reference_top_score":top,"reference_second_score":second,"reference_margin":margin,"reference_combined_confidence_score":combined,"reference_confidence_label":conf})
    return pd.DataFrame(out)
