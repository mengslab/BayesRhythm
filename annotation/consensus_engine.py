import pandas as pd

def build_consensus_labels(marker_df=None, reference_df=None):
    if (marker_df is None or marker_df.empty) and (reference_df is None or reference_df.empty):
        return pd.DataFrame()

    base = None
    for df in [marker_df, reference_df]:
        if df is not None and not df.empty:
            base = df[["cell_id"]].copy()
            break

    if marker_df is not None and not marker_df.empty:
        keep = [c for c in [
            "cell_id", "marker_predicted_label", "marker_confidence_label",
            "marker_combined_confidence_score"
        ] if c in marker_df.columns]
        base = base.merge(marker_df[keep], on="cell_id", how="left")

    if reference_df is not None and not reference_df.empty:
        keep = [c for c in [
            "cell_id", "reference_predicted_label", "reference_confidence_label",
            "reference_combined_confidence_score"
        ] if c in reference_df.columns]
        base = base.merge(reference_df[keep], on="cell_id", how="left")

    def score_label(conf):
        conf = str(conf).lower()
        return {"high": 3, "moderate": 2, "low": 1}.get(conf, 0)

    rows = []
    for _, r in base.iterrows():
        m_label = str(r.get("marker_predicted_label", "unknown"))
        r_label = str(r.get("reference_predicted_label", "unknown"))
        m_conf = str(r.get("marker_confidence_label", "low"))
        r_conf = str(r.get("reference_confidence_label", "low"))
        m_score = score_label(m_conf)
        r_score = score_label(r_conf)

        if m_label == r_label and m_label not in ["unknown", "", "nan"]:
            final = m_label
            rationale = "agree"
            conf = "high" if min(m_score, r_score) >= 2 else "moderate"
        elif r_label not in ["unknown", "", "nan"] and r_score > m_score:
            final = r_label
            rationale = "reference_preferred"
            conf = r_conf
        elif m_label not in ["unknown", "", "nan"] and m_score > r_score:
            final = m_label
            rationale = "marker_preferred"
            conf = m_conf
        elif r_label not in ["unknown", "", "nan"] and m_label in ["unknown", "", "nan"]:
            final = r_label
            rationale = "reference_only"
            conf = r_conf
        elif m_label not in ["unknown", "", "nan"] and r_label in ["unknown", "", "nan"]:
            final = m_label
            rationale = "marker_only"
            conf = m_conf
        else:
            final = "unknown"
            rationale = "unresolved"
            conf = "low"

        rows.append({
            "cell_id": r["cell_id"],
            "consensus_cell_type": final,
            "consensus_rationale": rationale,
            "consensus_confidence_label": conf,
            "marker_predicted_label": m_label,
            "reference_predicted_label": r_label,
        })

    return pd.DataFrame(rows)

def summarize_consensus(consensus_df):
    if consensus_df is None or consensus_df.empty:
        return pd.DataFrame()
    return consensus_df["consensus_cell_type"].astype(str).value_counts().rename_axis("consensus_cell_type").reset_index(name="n_cells")
