from pathlib import Path
import gzip, tarfile, tempfile
import pandas as pd
from scipy.io import mmread
from bayesrhythm_core_models import SingleCellReferenceInput, ImportAuditTrail

def _open_text(path):
    path=Path(path)
    return gzip.open(path,"rt") if path.suffix==".gz" else open(path,"r",encoding="utf-8")

def _make_unique(names):
    counts={}; out=[]
    for n in names:
        if n not in counts: counts[n]=0; out.append(n)
        else: counts[n]+=1; out.append(f"{n}_{counts[n]}")
    return out

def detect_mex_triplet(input_dir):
    input_dir=Path(input_dir); files={p.name:p for p in input_dir.iterdir() if p.is_file()}
    def pick(cands):
        for c in cands:
            if c in files: return files[c]
        return None
    return {"matrix":pick(["matrix.mtx","matrix.mtx.gz"]),"barcodes":pick(["barcodes.tsv","barcodes.tsv.gz"]),"features":pick(["features.tsv","features.tsv.gz","genes.tsv","genes.tsv.gz"])}

def extract_mex_archive(archive_path):
    archive_path=Path(archive_path)
    tmpdir=Path(tempfile.mkdtemp())
    with tarfile.open(archive_path, "r:gz") as tf:
        tf.extractall(tmpdir)
    return tmpdir

def read_features_table(path):
    rows=[]
    with _open_text(path) as f:
        for line in f: rows.append(line.rstrip("\n").split("\t"))
    max_cols=max((len(r) for r in rows), default=0); cols=[f"col_{i+1}" for i in range(max_cols)]; df=pd.DataFrame(rows, columns=cols); out=pd.DataFrame()
    out["feature_id"]=df["col_1"] if "col_1" in df else pd.Series(dtype=str); out["gene_symbol"]=df["col_2"] if "col_2" in df else out["feature_id"]; out["feature_type"]=df["col_3"] if "col_3" in df else "Gene Expression"; out["feature_source"]=df["col_4"] if "col_4" in df else ""; out["chromosome"]=df["col_5"] if "col_5" in df else ""; out["extra_1"]=df["col_6"] if "col_6" in df else ""; return out

def read_barcodes_table(path):
    with _open_text(path) as f: vals=[line.strip() for line in f if line.strip()]
    return pd.DataFrame({"cell_id":vals})

def read_mex_matrix(matrix_path, features_df, barcodes_df):
    mat=mmread(str(matrix_path)).tocsr()
    if mat.shape[0]==len(features_df) and mat.shape[1]==len(barcodes_df): mat=mat.T.tocsr()
    genes=_make_unique(features_df["gene_symbol"].astype(str).tolist())
    expr=pd.DataFrame.sparse.from_spmatrix(mat, index=barcodes_df["cell_id"].astype(str).tolist(), columns=genes)
    return expr, features_df.copy(), barcodes_df.copy()

def load_10x_mex_directory(input_dir, dataset_metadata=None):
    triplet=detect_mex_triplet(input_dir)
    if not all(triplet.values()): raise FileNotFoundError("Missing MEX components.")
    features=read_features_table(triplet["features"]); barcodes=read_barcodes_table(triplet["barcodes"]); expr,gene_meta,cell_meta=read_mex_matrix(triplet["matrix"], features, barcodes)
    cell_meta["condition"]=""; cell_meta["zeitgeber_time"]=float("nan"); cell_meta["cell_type"]="unknown"
    return SingleCellReferenceInput(expr, cell_meta, gene_meta, dataset_metadata or {}, audit_trail=ImportAuditTrail(source_files=[str(triplet["matrix"]),str(triplet["barcodes"]),str(triplet["features"])]))

def load_10x_mex_archive(archive_path, dataset_metadata=None):
    tmpdir=extract_mex_archive(archive_path)
    sc=load_10x_mex_directory(tmpdir, dataset_metadata=dataset_metadata or {"mode":"mex_archive"})
    sc.audit_trail.source_files.insert(0, str(archive_path))
    return sc
