import pandas as pd
from br_io.bulk_auto_parser import parse_joined_bulk_excel, summarize_bulk_metadata

def test_bulk_auto_parser_smoke():
    df = pd.DataFrame({
        "GeneID": ["g1", "g2"],
        "Base mean": [1, 2],
        "Chow SRR1 ZT0.fastqsanger.tabular": [10, 20],
        "TRF SRR2 ZT8.fastqsanger.tabular": [30, 40],
    })
    expr, meta = parse_joined_bulk_excel(df)
    assert expr.shape == (2, 2)
    assert meta.shape[0] == 2
    summary = summarize_bulk_metadata(meta)
    assert not summary.empty
