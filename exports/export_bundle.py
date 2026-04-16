from pathlib import Path
from datetime import datetime
import zipfile, pandas as pd
def export_results_bundle(base_dir, tables):
    base=Path(base_dir); base.mkdir(parents=True, exist_ok=True); stamp=datetime.now().strftime("%Y%m%d_%H%M%S"); out_dir=base/f"results_{stamp}"; out_dir.mkdir(parents=True, exist_ok=True); written=[]
    for name,df in tables.items():
        if isinstance(df,pd.DataFrame) and not df.empty:
            p=out_dir/f"{name}.csv"; df.to_csv(p,index=False); written.append(p)
    (out_dir/"manifest.txt").write_text("\n".join([p.name for p in written]), encoding="utf-8")
    zip_path=base/f"results_{stamp}.zip"
    with zipfile.ZipFile(zip_path,"w",zipfile.ZIP_DEFLATED) as zf:
        for p in out_dir.rglob("*"):
            if p.is_file(): zf.write(p, p.relative_to(out_dir.parent))
    return zip_path
