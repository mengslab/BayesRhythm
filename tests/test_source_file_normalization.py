import pandas as pd
from bayesrhythm_core_models import SingleCellReferenceInput
from anchors.multi_anchor_manager import MultiAnchorManager

def test_source_file_normalization_prefers_archive():
    expr = pd.DataFrame({"g1":[1,2]}, index=["c1","c2"])
    meta = pd.DataFrame({"cell_id":["c1","c2"],"condition":["TRF","TRF"],"zeitgeber_time":[20,20],"cell_type":["A","A"]})
    sc = SingleCellReferenceInput(expr, meta, pd.DataFrame())
    mgr = MultiAnchorManager()
    review = pd.DataFrame({"final_cell_type":["A","A"],"combined_annotation_confidence_label":["high","high"],"methods_agree":[True,True]})
    src = "/Users/me/TRF_ZT20_filtered_feature_bc_matrix.tar.gz, /tmp/matrix.mtx.gz, /tmp/barcodes.tsv.gz, /tmp/features.tsv.gz"
    rec = mgr.register_anchor("TRF_ZT20", sc, src, "mex_archive_batch_local_path", dataset_label="TRF_ZT20", validation_report=pd.DataFrame([{"severity":"ok"}]), review_table=review, keep_dataset=False)
    assert rec.source_file.endswith(".tar.gz")
    assert "," not in rec.source_file
