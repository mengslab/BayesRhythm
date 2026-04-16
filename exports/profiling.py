from pathlib import Path
from time import perf_counter
import json

class StepProfiler:
    def __init__(self):
        self.records = []

    def start(self, step_name: str):
        return {"step_name": step_name, "t0": perf_counter()}

    def stop(self, token, extra=None):
        elapsed = perf_counter() - token["t0"]
        rec = {"step_name": token["step_name"], "elapsed_seconds": round(elapsed, 4)}
        if extra:
            rec.update(extra)
        self.records.append(rec)
        return rec

    def to_dataframe(self):
        try:
            import pandas as pd
            return pd.DataFrame(self.records)
        except Exception:
            return self.records

    def save_json(self, output_path):
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(self.records, indent=2), encoding="utf-8")
        return output_path
