from dataclasses import dataclass
import pandas as pd
import numpy as np

@dataclass
class RhythmDiscoveryResult:
    gene_results: pd.DataFrame

@dataclass
class RhythmDiscoveryRunResult:
    result: RhythmDiscoveryResult
    previews: dict
    warnings: list

class RhythmDiscoveryEngine:
    def run(self, bulk_matrix, bulk_metadata, period_hours=24.0):
        if bulk_matrix is None or bulk_metadata is None or bulk_matrix.empty or bulk_metadata.empty:
            return RhythmDiscoveryRunResult(RhythmDiscoveryResult(pd.DataFrame()), {}, ["Bulk matrix/metadata missing."])
        if "sample_id" not in bulk_metadata.columns or "zeitgeber_time" not in bulk_metadata.columns:
            return RhythmDiscoveryRunResult(RhythmDiscoveryResult(pd.DataFrame()), {}, ["Bulk metadata must contain sample_id and zeitgeber_time."])
        md = bulk_metadata.copy()
        md["sample_id"] = md["sample_id"].astype(str)
        md = md[md["sample_id"].isin(bulk_matrix.columns)].copy()
        if md.shape[0] < 4:
            return RhythmDiscoveryRunResult(RhythmDiscoveryResult(pd.DataFrame()), {}, ["Need at least 4 bulk samples for rhythm fitting."])
        t = md["zeitgeber_time"].astype(float).values
        omega = 2 * np.pi / period_hours
        X = np.column_stack([np.ones_like(t), np.cos(omega * t), np.sin(omega * t)])
        results = []
        sample_order = md["sample_id"].tolist()
        for gene_id in bulk_matrix.index.astype(str):
            y = bulk_matrix.loc[gene_id, sample_order].astype(float).values
            try:
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                mesor, bcos, bsin = beta
                amp = float(np.sqrt(bcos**2 + bsin**2))
                phase = float((np.arctan2(-bsin, bcos) % (2*np.pi)) / omega)
                yhat = X @ beta
                ss_res = float(np.sum((y - yhat)**2))
                ss_tot = float(np.sum((y - np.mean(y))**2))
                r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 0.0
                results.append({"gene_id": gene_id, "mesor": float(mesor), "amplitude": amp, "phase_hours": phase, "fit_r2": r2, "rhythmic_call": "rhythmic" if r2 >= 0.3 else "weak"})
            except Exception:
                continue
        df = pd.DataFrame(results)
        return RhythmDiscoveryRunResult(RhythmDiscoveryResult(df), {"rhythm_preview": df.head(50)}, [])
