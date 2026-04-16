from dataclasses import dataclass
import pandas as pd

@dataclass
class RhythmComparisonResult:
    gene_results: pd.DataFrame

@dataclass
class RhythmComparisonRunResult:
    result: RhythmComparisonResult
    previews: dict
    warnings: list

class RhythmComparisonEngine:
    def run(self, bulk_matrix, bulk_metadata, period_hours=24.0):
        if bulk_matrix is None or bulk_metadata is None or bulk_matrix.empty or bulk_metadata.empty:
            return RhythmComparisonRunResult(RhythmComparisonResult(pd.DataFrame()), {}, ["Bulk matrix/metadata missing."])
        if "condition" not in bulk_metadata.columns or "sample_id" not in bulk_metadata.columns:
            return RhythmComparisonRunResult(RhythmComparisonResult(pd.DataFrame()), {}, ["Bulk metadata must contain condition and sample_id."])
        conditions = sorted(bulk_metadata["condition"].astype(str).dropna().unique().tolist())
        if len(conditions) < 2:
            return RhythmComparisonRunResult(RhythmComparisonResult(pd.DataFrame()), {}, ["Need at least two conditions for rhythm comparison."])
        c1, c2 = conditions[:2]
        md1 = bulk_metadata[bulk_metadata["condition"].astype(str) == c1].copy()
        md2 = bulk_metadata[bulk_metadata["condition"].astype(str) == c2].copy()
        common1 = [s for s in md1["sample_id"].astype(str).tolist() if s in bulk_matrix.columns]
        common2 = [s for s in md2["sample_id"].astype(str).tolist() if s in bulk_matrix.columns]
        if len(common1) < 3 or len(common2) < 3:
            return RhythmComparisonRunResult(RhythmComparisonResult(pd.DataFrame()), {}, ["Each condition needs at least 3 matched bulk samples."])
        rows = []
        for gene_id in bulk_matrix.index.astype(str):
            y1 = bulk_matrix.loc[gene_id, common1].astype(float)
            y2 = bulk_matrix.loc[gene_id, common2].astype(float)
            amp1 = float(y1.max() - y1.min()) / 2.0 if len(y1) else 0.0
            amp2 = float(y2.max() - y2.min()) / 2.0 if len(y2) else 0.0
            mesor1 = float(y1.mean()) if len(y1) else 0.0
            mesor2 = float(y2.mean()) if len(y2) else 0.0
            amp_delta = amp2 - amp1
            mesor_delta = mesor2 - mesor1
            call = "changed" if abs(amp_delta) > 0.1 * max(abs(amp1), 1e-8) or abs(mesor_delta) > 0.1 * max(abs(mesor1), 1e-8) else "stable"
            rows.append({"gene_id": gene_id, f"{c1}_amplitude_proxy": amp1, f"{c2}_amplitude_proxy": amp2, f"{c1}_mesor_proxy": mesor1, f"{c2}_mesor_proxy": mesor2, "amplitude_delta": amp_delta, "mesor_delta": mesor_delta, "comparison_call": call})
        df = pd.DataFrame(rows)
        return RhythmComparisonRunResult(RhythmComparisonResult(df), {"comparison_preview": df.head(50)}, [])
