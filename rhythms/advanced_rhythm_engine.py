from dataclasses import dataclass
import pandas as pd
import numpy as np

@dataclass
class AdvancedRhythmResult:
    gene_results: pd.DataFrame

@dataclass
class AdvancedRhythmRunResult:
    result: AdvancedRhythmResult
    previews: dict
    warnings: list

def _fit_model(t, y, period_hours=24.0, two_harmonic=False, linear_trend=False):
    omega = 2 * np.pi / period_hours
    cols = [np.ones_like(t), np.cos(omega*t), np.sin(omega*t)]
    names = ["intercept", "cos1", "sin1"]
    if two_harmonic:
        cols.extend([np.cos(2*omega*t), np.sin(2*omega*t)])
        names.extend(["cos2", "sin2"])
    if linear_trend:
        cols.append(t)
        names.append("trend")
    X = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat)**2))
    ss_tot = float(np.sum((y - np.mean(y))**2))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 0.0
    amp1 = float(np.sqrt(beta[names.index("cos1")]**2 + beta[names.index("sin1")]**2))
    phase1 = float((np.arctan2(-beta[names.index("sin1")], beta[names.index("cos1")]) % (2*np.pi)) / omega)
    return {
        "fit_r2": r2,
        "mesor": float(beta[names.index("intercept")]),
        "amplitude": amp1,
        "phase_hours": phase1,
        "model": "two_harmonic" if two_harmonic else "single_harmonic",
        "has_trend": linear_trend,
    }

class AdvancedRhythmEngine:
    def run(self, bulk_matrix, bulk_metadata, period_hours=24.0):
        if bulk_matrix is None or bulk_metadata is None or bulk_matrix.empty or bulk_metadata.empty:
            return AdvancedRhythmRunResult(AdvancedRhythmResult(pd.DataFrame()), {}, ["Bulk matrix/metadata missing."])
        if "sample_id" not in bulk_metadata.columns or "zeitgeber_time" not in bulk_metadata.columns:
            return AdvancedRhythmRunResult(AdvancedRhythmResult(pd.DataFrame()), {}, ["Bulk metadata must contain sample_id and zeitgeber_time."])
        md = bulk_metadata.copy()
        md["sample_id"] = md["sample_id"].astype(str)
        md = md[md["sample_id"].isin(bulk_matrix.columns)].copy()
        if md.shape[0] < 6:
            return AdvancedRhythmRunResult(AdvancedRhythmResult(pd.DataFrame()), {}, ["Need at least 6 bulk samples for advanced rhythm fitting."])
        t = md["zeitgeber_time"].astype(float).values
        sample_order = md["sample_id"].tolist()
        rows = []
        for gene_id in bulk_matrix.index.astype(str):
            y = bulk_matrix.loc[gene_id, sample_order].astype(float).values
            try:
                m1 = _fit_model(t, y, period_hours=period_hours, two_harmonic=False, linear_trend=False)
                m2 = _fit_model(t, y, period_hours=period_hours, two_harmonic=True, linear_trend=False)
                m3 = _fit_model(t, y, period_hours=period_hours, two_harmonic=True, linear_trend=True)
                best = max([m1, m2, m3], key=lambda x: x["fit_r2"])
                rows.append({
                    "gene_id": gene_id,
                    "best_model": best["model"] + ("+trend" if best["has_trend"] else ""),
                    "mesor": best["mesor"],
                    "amplitude": best["amplitude"],
                    "phase_hours": best["phase_hours"],
                    "fit_r2": best["fit_r2"],
                    "rhythmic_call": "rhythmic" if best["fit_r2"] >= 0.35 else ("weak" if best["fit_r2"] >= 0.2 else "nonrhythmic"),
                })
            except Exception:
                continue
        df = pd.DataFrame(rows)
        return AdvancedRhythmRunResult(AdvancedRhythmResult(df), {"advanced_rhythm_preview": df.head(50)}, [])
