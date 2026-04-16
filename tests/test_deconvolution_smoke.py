import pandas as pd
from trajectories.state_aware_deconvolution_engine import StateAwareDeconvolutionEngine
from types import SimpleNamespace

def test_state_aware_deconv_smoke():
    ref = pd.DataFrame({
        "anchor_id":["A","A"], "condition":["Chow","Chow"], "zeitgeber_time":[8,8], "cell_type":["Hep","Kupffer"], "n_cells":[10,10],
        "g1":[1.0,0.0], "g2":[0.0,1.0]
    })
    atlas = SimpleNamespace(shared_genes=["g1","g2"], reference_table=ref)
    bulk = pd.DataFrame({"s1":[0.8,0.2], "s2":[0.3,0.7]}, index=["g1","g2"])
    md = pd.DataFrame({"sample_id":["s1","s2"], "condition":["Chow","Chow"], "zeitgeber_time":[8,20]})
    res = StateAwareDeconvolutionEngine().run(bulk, md, atlas, min_shared_genes=2)
    assert not res.result.cell_fractions.empty
