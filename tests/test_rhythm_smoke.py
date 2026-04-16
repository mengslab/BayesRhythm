import pandas as pd
from rhythms.advanced_rhythm_engine import AdvancedRhythmEngine

def test_advanced_rhythm_smoke():
    bulk = pd.DataFrame({
        "s1":[1.0,2.0], "s2":[2.0,2.5], "s3":[3.0,2.0],
        "s4":[2.0,1.5], "s5":[1.0,1.0], "s6":[2.0,1.5]
    }, index=["g1","g2"])
    md = pd.DataFrame({
        "sample_id":["s1","s2","s3","s4","s5","s6"],
        "condition":["Chow"]*6,
        "zeitgeber_time":[0,4,8,12,16,20]
    })
    res = AdvancedRhythmEngine().run(bulk, md)
    assert not res.result.gene_results.empty
