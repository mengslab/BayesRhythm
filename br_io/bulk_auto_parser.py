import re
import pandas as pd

SUMMARY_COLUMNS = ["GeneID", "Base mean", "log2(FC)", "StdErr", "Wald-Stats", "P-value", "P-adj", "Chow", "TRF"]

def parse_joined_bulk_excel(df: pd.DataFrame):
    cols = [str(c) for c in df.columns]
    gene_col = "GeneID" if "GeneID" in cols else cols[0]
    sample_cols = [c for c in cols if c not in SUMMARY_COLUMNS and c != gene_col]
    expr = df[[gene_col] + sample_cols].copy()
    expr = expr.rename(columns={gene_col: "gene_id"}).set_index("gene_id")
    meta_rows = []
    for c in sample_cols:
        cond = "Chow" if "Chow" in c else ("TRF" if "TRF" in c else None)
        m = re.search(r'ZT(\d+)', c)
        zt = int(m.group(1)) if m else None
        m2 = re.search(r'(SRR\d+)', c)
        accession = m2.group(1) if m2 else ""
        meta_rows.append({
            "sample_id": c,
            "condition": cond if cond is not None else "unknown",
            "zeitgeber_time": zt,
            "accession": accession,
            "source_label": c,
        })
    meta = pd.DataFrame(meta_rows)
    return expr, meta

def merge_bulk_metadata_overlay(metadata_df: pd.DataFrame, overlay_df: pd.DataFrame):
    if metadata_df is None or metadata_df.empty or overlay_df is None or overlay_df.empty:
        return metadata_df
    out = metadata_df.copy()
    ov = overlay_df.copy()
    # accepted keys: sample_id or accession
    if "sample_id" in ov.columns:
        out = out.drop(columns=["condition", "zeitgeber_time"], errors="ignore").merge(ov, on="sample_id", how="left")
    elif "accession" in ov.columns:
        out = out.drop(columns=["condition", "zeitgeber_time"], errors="ignore").merge(ov, on="accession", how="left")
    else:
        return metadata_df
    if "condition" not in out.columns:
        out["condition"] = "unknown"
    if "zeitgeber_time" not in out.columns:
        out["zeitgeber_time"] = pd.NA
    out["condition"] = out["condition"].fillna("unknown")
    return out

def summarize_bulk_metadata(metadata_df: pd.DataFrame):
    if metadata_df is None or metadata_df.empty:
        return pd.DataFrame()
    return metadata_df.groupby(["condition", "zeitgeber_time"]).size().reset_index(name="n_samples").sort_values(["condition", "zeitgeber_time"])
