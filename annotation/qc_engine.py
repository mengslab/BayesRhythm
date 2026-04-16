import pandas as pd, numpy as np
def compute_single_cell_qc(sc_input):
    expr=sc_input.expression_matrix
    # sparse-safe computations
    sparse_arr = expr.sparse.to_coo().tocsr()
    n_genes = np.diff(sparse_arr.indptr)
    lib_size = np.asarray(sparse_arr.sum(axis=1)).ravel()
    meta=sc_input.cell_metadata.copy().set_index("cell_id").reindex(expr.index.astype(str)).reset_index()
    summary=pd.DataFrame([{"n_cells":int(expr.shape[0]),"n_genes":int(expr.shape[1]),"median_detected_genes_per_cell":float(np.median(n_genes)) if len(n_genes) else 0.0,"median_library_size":float(np.median(lib_size)) if len(lib_size) else 0.0,"unknown_label_fraction":float(meta["cell_type"].astype(str).str.lower().eq("unknown").mean()) if "cell_type" in meta.columns else 1.0}])
    per_cell=pd.DataFrame({"cell_id":expr.index.astype(str),"detected_genes":n_genes,"library_size":lib_size})
    label_counts=meta["cell_type"].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n_cells") if "cell_type" in meta.columns else pd.DataFrame()
    condition_counts=meta.groupby(["condition","zeitgeber_time"]).size().reset_index(name="n_cells") if {"condition","zeitgeber_time"}.issubset(meta.columns) else pd.DataFrame()
    return {"dataset_summary":summary,"per_cell_qc":per_cell,"label_counts":label_counts,"condition_counts":condition_counts}
