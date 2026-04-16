import pandas as pd
from rhythms.rhythm_workflow import summarize_rhythm_results
from decomposition.decomposition_refactor import build_unified_decomposition_table

def test_rhythm_summary_smoke():
    s = pd.DataFrame({"gene_id":["g1","g2"],"rhythmic_call":["rhythmic","weak"]})
    c = pd.DataFrame({"gene_id":["g1","g2"],"comparison_call":["changed","stable"]})
    out = summarize_rhythm_results(simple_df=s, comparison_df=c)
    assert not out.summary_table.empty

def test_unified_decomposition_smoke():
    a = pd.DataFrame({"gene_id":["g1","g2"],"composition_score":[0.8,0.2]})
    b = pd.DataFrame({"gene_id":["g1","g2"],"state_program_r2":[0.7,0.1]})
    out = build_unified_decomposition_table(a, b)
    assert not out.unified_gene_table.empty
