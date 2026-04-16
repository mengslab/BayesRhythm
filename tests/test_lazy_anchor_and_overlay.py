import pandas as pd
from bayesrhythm_core_models import SingleCellReferenceInput
from anchors.multi_anchor_manager import MultiAnchorManager
from br_io.bulk_auto_parser import merge_bulk_metadata_overlay

def test_overlay_merge_smoke():
    md = pd.DataFrame({"sample_id":["a","b"],"condition":["unknown","unknown"],"zeitgeber_time":[None,None],"accession":["SRR1","SRR2"]})
    ov = pd.DataFrame({"accession":["SRR1","SRR2"],"condition":["Chow","TRF"],"zeitgeber_time":[0,8]})
    out = merge_bulk_metadata_overlay(md, ov)
    assert list(out["condition"]) == ["Chow","TRF"]

def test_manager_unique_zt_smoke():
    expr = pd.DataFrame({"g1":[1,2],"g2":[3,4]}, index=["c1","c2"])
    meta = pd.DataFrame({"cell_id":["c1","c2"],"condition":["Chow","Chow"],"zeitgeber_time":[8,8],"cell_type":["A","A"]})
    sc = SingleCellReferenceInput(expr, meta, pd.DataFrame())
    mgr = MultiAnchorManager()
    review = pd.DataFrame({"final_cell_type":["A","A"],"combined_annotation_confidence_label":["high","high"],"methods_agree":[True,True]})
    rec = mgr.register_anchor("Chow_ZT8", sc, "dummy.tar.gz", "mex_archive_batch_local_path", dataset_label="Chow_ZT8", validation_report=pd.DataFrame([{"severity":"ok"}]), review_table=review, keep_dataset=True)
    assert rec.zeitgeber_time == 8.0
