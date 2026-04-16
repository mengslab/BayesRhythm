import pandas as pd

def compute_anchor_preflight(sc_input, review_table=None, validation_report=None, min_cells=500, min_detected_genes_median=200):
    expr = sc_input.expression_matrix
    n_cells = int(expr.shape[0])
    n_genes = int(expr.shape[1])

    severity = "ready"
    reasons = []

    if validation_report is not None and not validation_report.empty:
        if (validation_report["severity"] == "error").any():
            severity = "blocked"
            reasons.append("validation_error")
        elif (validation_report["severity"] == "warning").any() and severity != "blocked":
            severity = "warning"
            reasons.append("validation_warning")

    if n_cells < min_cells:
        severity = "blocked"
        reasons.append("low_cell_count")

    # sparse-safe detected genes per cell
    try:
        sparse_arr = expr.sparse.to_coo().tocsr()
        detected = pd.Series(sparse_arr.indptr[1:] - sparse_arr.indptr[:-1])
        median_detected = float(detected.median()) if len(detected) else 0.0
    except Exception:
        median_detected = 0.0
        reasons.append("qc_compute_failed")
        if severity == "ready":
            severity = "warning"

    if median_detected < min_detected_genes_median:
        if severity != "blocked":
            severity = "warning"
        reasons.append("low_median_detected_genes")

    unresolved_fraction = 1.0
    low_conf_fraction = 1.0
    disagreement_fraction = 1.0

    if review_table is not None and not review_table.empty:
        if "final_cell_type" in review_table.columns:
            unresolved_fraction = float(review_table["final_cell_type"].astype(str).str.lower().eq("unknown").mean())
        if "combined_annotation_confidence_label" in review_table.columns:
            low_conf_fraction = float(review_table["combined_annotation_confidence_label"].astype(str).str.lower().eq("low").mean())
        if "methods_agree" in review_table.columns:
            disagreement_fraction = float((~review_table["methods_agree"].fillna(False)).mean())

        if unresolved_fraction > 0.3:
            severity = "blocked"
            reasons.append("too_many_unresolved_labels")
        elif unresolved_fraction > 0.15 and severity != "blocked":
            severity = "warning"
            reasons.append("moderate_unresolved_labels")

        if low_conf_fraction > 0.4 and severity != "blocked":
            severity = "warning"
            reasons.append("many_low_confidence_labels")

        if disagreement_fraction > 0.5 and severity != "blocked":
            severity = "warning"
            reasons.append("high_method_disagreement")

    return pd.DataFrame([{
        "n_cells": n_cells,
        "n_genes": n_genes,
        "median_detected_genes_per_cell": median_detected,
        "unresolved_fraction": unresolved_fraction,
        "low_confidence_fraction": low_conf_fraction,
        "disagreement_fraction": disagreement_fraction,
        "preflight_status": severity,
        "preflight_reasons": "; ".join(reasons) if reasons else "passed",
    }])


def compute_panel_preflight_summary(anchor_manager):
    rows = []
    panel = anchor_manager.get_anchor_panel()
    for record in panel.anchor_records:
        rows.append({
            "dataset_id": record.dataset_id,
            "dataset_label": record.dataset_label,
            "condition": record.condition,
            "zeitgeber_time": record.zeitgeber_time,
            "n_cells": record.n_cells,
            "n_genes": record.n_genes,
            "annotation_status": record.annotation_status,
            "validation_status": record.validation_status,
            "readiness_status": record.readiness_status,
            "low_confidence_fraction": record.low_confidence_fraction,
            "unresolved_fraction": record.unresolved_fraction,
            "disagreement_fraction": record.disagreement_fraction,
        })
    return pd.DataFrame(rows)

def compute_panel_preflight_status_counts(panel_df):
    if panel_df is None or panel_df.empty or "readiness_status" not in panel_df.columns:
        return pd.DataFrame()
    return panel_df["readiness_status"].astype(str).value_counts().rename_axis("readiness_status").reset_index(name="n_anchors")
