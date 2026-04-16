import pandas as pd
from bayesrhythm_core_models import SingleCellReferenceInput
from annotation.marker_annotation_engine import score_markers_per_cell, assign_marker_labels
from annotation.marker_library import MARKER_LIBRARY

def test_marker_annotation_smoke():
    expr = pd.DataFrame({"Alb":[1,0],"Ttr":[1,0],"Adgre1":[0,1]}, index=["c1","c2"])
    meta = pd.DataFrame({"cell_id":["c1","c2"],"condition":["Chow","Chow"],"zeitgeber_time":[8,8],"cell_type":["unknown","unknown"]})
    sc = SingleCellReferenceInput(expr, meta, pd.DataFrame())
    scores = score_markers_per_cell(sc, MARKER_LIBRARY, use_log1p=False)
    assigns = assign_marker_labels(scores)
    assert not assigns.empty
