import pandas as pd
from bayesrhythm_core_models import SingleCellReferenceInput
from anchors.multi_anchor_manager import MultiAnchorManager
from references.multi_anchor_reference_builder import MultiAnchorReferenceBuilder
from trajectories.state_aware_deconvolution_engine import StateAwareDeconvolutionEngine
from rhythms.advanced_rhythm_engine import AdvancedRhythmEngine
from decomposition.state_program_decomposition_engine import StateProgramDecompositionEngine

def test_connected_workflow_smoke():
    expr = pd.DataFrame({"g1":[1,1,0,0],"g2":[0,0,1,1]}, index=["c1","c2","c3","c4"])
    meta = pd.DataFrame({"cell_id":["c1","c2","c3","c4"],"condition":["Chow"]*4,"zeitgeber_time":[8]*4,"cell_type":["A","A","B","B"]})
    sc = SingleCellReferenceInput(expr, meta, pd.DataFrame())
    mgr = MultiAnchorManager()
    review = pd.DataFrame({"final_cell_type":["A","A","B","B"],"combined_annotation_confidence_label":["high"]*4,"methods_agree":[True]*4})
    mgr.register_anchor("Chow_ZT8", sc, "x", "test", dataset_label="Chow_ZT8", validation_report=pd.DataFrame([{"severity":"ok"}]), review_table=review)
    ref = MultiAnchorReferenceBuilder().run(mgr).atlas
    bulk = pd.DataFrame({"s1":[1.0,0.0],"s2":[0.5,0.5],"s3":[0.0,1.0]}, index=["g1","g2"])
    bulk_md = pd.DataFrame({"sample_id":["s1","s2","s3"],"condition":["Chow","TRF","TRF"],"zeitgeber_time":[0,8,16]})
    deconv = StateAwareDeconvolutionEngine().run(bulk, bulk_md, ref, min_shared_genes=2).result
    assert not deconv.cell_fractions.empty
    adv = AdvancedRhythmEngine().run(bulk, bulk_md).result
    # allow empty if too few points, but engine should return object
    assert adv is not None
    sp = StateProgramDecompositionEngine().run(bulk, deconv, rhythm_results=adv.gene_results).result
    assert sp is not None
