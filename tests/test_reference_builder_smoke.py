import pandas as pd
from bayesrhythm_core_models import SingleCellReferenceInput
from anchors.multi_anchor_manager import MultiAnchorManager
from references.multi_anchor_reference_builder import MultiAnchorReferenceBuilder

def test_reference_builder_smoke():
    expr = pd.DataFrame({"g1":[1,2,0,0],"g2":[0,0,1,2]}, index=["c1","c2","c3","c4"])
    meta = pd.DataFrame({"cell_id":["c1","c2","c3","c4"],"condition":["Chow"]*4,"zeitgeber_time":[8]*4,"cell_type":["A","A","B","B"]})
    sc = SingleCellReferenceInput(expr, meta, pd.DataFrame())
    mgr = MultiAnchorManager()
    review = pd.DataFrame({"final_cell_type":["A","A","B","B"],"combined_annotation_confidence_label":["high"]*4,"methods_agree":[True]*4})
    mgr.register_anchor("Chow_ZT8", sc, "x", "test", dataset_label="Chow_ZT8", validation_report=pd.DataFrame([{"severity":"ok"}]), review_table=review)
    res = MultiAnchorReferenceBuilder().run(mgr)
    assert res.atlas.reference_table.shape[0] >= 2
