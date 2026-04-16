import pandas as pd
def combine_annotation_confidence(marker_df=None, reference_df=None):
    base=None
    for df in [marker_df, reference_df]:
        if df is not None and not df.empty: base=df[["cell_id"]].copy(); break
    if base is None: return pd.DataFrame()
    if marker_df is not None and not marker_df.empty: base=base.merge(marker_df[["cell_id","marker_predicted_label","marker_combined_confidence_score","marker_confidence_label"]], on="cell_id", how="left")
    if reference_df is not None and not reference_df.empty: base=base.merge(reference_df[["cell_id","reference_predicted_label","reference_combined_confidence_score","reference_confidence_label"]], on="cell_id", how="left")
    base["methods_agree"]=base.apply(lambda r: str(r.get("marker_predicted_label",""))==str(r.get("reference_predicted_label","")) and str(r.get("marker_predicted_label","")) not in ["","unknown","nan"], axis=1)
    m=base.get("marker_combined_confidence_score", pd.Series([0.0]*len(base))); r=base.get("reference_combined_confidence_score", pd.Series([0.0]*len(base))); base["combined_annotation_confidence_score"]=(0.5*m.fillna(0.0)+0.5*r.fillna(0.0)+(0.1*base["methods_agree"].astype(float))).clip(0,1); base["combined_annotation_confidence_label"]=base["combined_annotation_confidence_score"].apply(lambda x:"high" if x>=0.75 else ("moderate" if x>=0.45 else "low")); return base
