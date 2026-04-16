import pandas as pd

def build_validation_report(issues):
    if not issues:
        return pd.DataFrame([{"severity":"ok","component":"validation","message":"No validation issues detected.","recommended_fix":""}])
    return pd.DataFrame([{
        "severity": sev,
        "component": comp,
        "message": msg,
        "recommended_fix": "Review uploaded files and metadata." if sev == "error" else "Check carefully before continuing."
    } for sev, comp, msg in issues])

def validate_single_cell_input(sc_input):
    issues=[]; expr=sc_input.expression_matrix; meta=sc_input.cell_metadata
    if expr is None or getattr(expr, "empty", False):
        issues.append(("error","expression_matrix","Expression matrix is empty."))
        return build_validation_report(issues)
    if expr.index.duplicated().any():
        issues.append(("error","expression_matrix","Duplicate cell IDs found in expression matrix index."))
    if expr.columns.duplicated().any():
        issues.append(("warning","expression_matrix","Duplicate gene names found in expression matrix columns."))
    if expr.shape[0] < 50:
        issues.append(("warning","expression_matrix",f"Very low cell count detected: {expr.shape[0]}"))
    if expr.shape[1] < 100:
        issues.append(("warning","expression_matrix",f"Very low gene count detected: {expr.shape[1]}"))
    if expr.shape[1] > 200000:
        issues.append(("warning","expression_matrix",f"Extremely large gene dimension detected: {expr.shape[1]} genes."))
    if meta is None or meta.empty:
        issues.append(("error","cell_metadata","Cell metadata is empty."))
    else:
        for col in ["cell_id","condition","zeitgeber_time","cell_type"]:
            if col not in meta.columns:
                sev = "error" if col == "cell_id" else "warning"
                issues.append((sev,"cell_metadata",f"Missing expected column: {col}"))
        if "cell_id" in meta.columns:
            if meta["cell_id"].astype(str).duplicated().any():
                issues.append(("error","cell_metadata","Duplicate cell_id values found in metadata."))
            expr_ids = set(map(str, expr.index))
            meta_ids = set(meta["cell_id"].astype(str))
            if expr_ids != meta_ids:
                issues.append(("warning","alignment",f"Expression and metadata cell IDs differ: {len(expr_ids - meta_ids)} only in matrix, {len(meta_ids - expr_ids)} only in metadata."))
        if "condition" in meta.columns:
            missing_condition = meta["condition"].astype(str).replace("", pd.NA).isna().mean()
            if missing_condition > 0:
                issues.append(("warning","cell_metadata",f"{missing_condition:.1%} cells have missing condition labels."))
        if "zeitgeber_time" in meta.columns:
            zt = pd.to_numeric(meta["zeitgeber_time"], errors="coerce")
            missing_zt = zt.isna().mean()
            if missing_zt > 0:
                issues.append(("warning","cell_metadata",f"{missing_zt:.1%} cells have missing zeitgeber_time labels."))
            bad_zt = (~zt.isna()) & ((zt < 0) | (zt > 24))
            if bad_zt.any():
                issues.append(("warning","cell_metadata","Found zeitgeber_time values outside [0, 24]."))
        if "cell_type" in meta.columns:
            unknown_fraction = meta["cell_type"].astype(str).str.lower().isin(["unknown","","nan"]).mean()
            if unknown_fraction > 0.8:
                issues.append(("warning","cell_metadata",f"High unknown cell_type fraction: {unknown_fraction:.1%}"))
    return build_validation_report(issues)

def validate_anchor_metadata(cell_metadata):
    issues=[]
    for col in ["cell_id","condition","zeitgeber_time"]:
        if col not in cell_metadata.columns:
            issues.append(("error","anchor_metadata",f"Missing required column: {col}"))
    if "condition" in cell_metadata.columns:
        vals = cell_metadata["condition"].astype(str).replace("", pd.NA).dropna().unique().tolist()
        if len(vals) != 1:
            issues.append(("error","anchor_metadata",f"Anchor upload must contain exactly one condition. Found {vals}"))
    if "zeitgeber_time" in cell_metadata.columns:
        vals = pd.to_numeric(cell_metadata["zeitgeber_time"], errors="coerce").dropna().unique().tolist()
        if len(vals) != 1:
            issues.append(("error","anchor_metadata",f"Anchor upload must contain exactly one zeitgeber time. Found {vals}"))
    return build_validation_report(issues)

def validate_bulk_input(expression_df, metadata_df):
    issues=[]
    if expression_df is None or expression_df.empty:
        issues.append(("error","bulk_expression","Bulk expression matrix is empty."))
        return build_validation_report(issues)
    if metadata_df is None or metadata_df.empty:
        issues.append(("error","bulk_metadata","Bulk metadata is empty."))
        return build_validation_report(issues)
    for col in ["sample_id","condition","zeitgeber_time"]:
        if col not in metadata_df.columns:
            issues.append(("error","bulk_metadata",f"Missing required column: {col}"))
    if "sample_id" in metadata_df.columns:
        if metadata_df["sample_id"].astype(str).duplicated().any():
            issues.append(("error","bulk_metadata","Duplicate sample_id values found."))
        expr_ids = set(map(str, expression_df.columns))
        meta_ids = set(metadata_df["sample_id"].astype(str))
        if expr_ids != meta_ids:
            issues.append(("warning","alignment",f"Bulk metadata and expression columns differ: {len(expr_ids-meta_ids)} only in matrix, {len(meta_ids-expr_ids)} only in metadata."))
    return build_validation_report(issues)
