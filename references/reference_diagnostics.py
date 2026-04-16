import pandas as pd
def compute_reference_diagnostics(sc_input, marker_df=None, reference_df=None):
    outputs={}; meta=sc_input.cell_metadata.copy()
    if "cell_type" in meta.columns: outputs["label_support"]=meta["cell_type"].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n_cells")
    if marker_df is not None and not marker_df.empty and reference_df is not None and not reference_df.empty:
        merged=marker_df[["cell_id","marker_predicted_label"]].merge(reference_df[["cell_id","reference_predicted_label"]], on="cell_id", how="inner"); merged["agree"]=merged["marker_predicted_label"].astype(str)==merged["reference_predicted_label"].astype(str); outputs["annotation_agreement_summary"]=merged["agree"].value_counts().rename_axis("agree").reset_index(name="n_cells")
    return outputs
