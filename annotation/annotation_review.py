from dataclasses import dataclass
import streamlit as st
import pandas as pd
@dataclass
class AnnotationReviewState:
    review_table: pd.DataFrame
    summary: dict
def merge_annotation_sources(cell_metadata, marker_df=None, reference_df=None, combined_conf_df=None):
    base=cell_metadata[["cell_id"]].copy()
    for df in [marker_df, reference_df, combined_conf_df]:
        if df is not None and not df.empty: base=base.merge(df, on="cell_id", how="left")
    if "reference_predicted_label" in base.columns: base["final_cell_type"]=base["reference_predicted_label"].fillna(base.get("marker_predicted_label","unknown"))
    elif "marker_predicted_label" in base.columns: base["final_cell_type"]=base["marker_predicted_label"].fillna("unknown")
    else: base["final_cell_type"]="unknown"
    base["review_status"]="pending"; base["review_comment"]=""; return base
def initialize_annotation_review(review_table): return AnnotationReviewState(review_table.copy(), {"total_cells":len(review_table)})
def render_annotation_review_panel(review_state): review_state.review_table=st.data_editor(review_state.review_table, use_container_width=True, num_rows="fixed"); return review_state
def apply_reviewed_labels(cell_metadata, review_state):
    final=review_state.review_table[["cell_id","final_cell_type"]].rename(columns={"final_cell_type":"cell_type"}); merged=cell_metadata.drop(columns=["cell_type"], errors="ignore").merge(final,on="cell_id",how="left"); merged["cell_type"]=merged["cell_type"].fillna("unknown"); return merged
