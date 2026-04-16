from pathlib import Path
def save_matplotlib_figure(fig, output_path_base: str):
    out = {}
    base = Path(output_path_base)
    base.parent.mkdir(parents=True, exist_ok=True)
    for ext in [".png", ".pdf", ".svg"]:
        path = base.with_suffix(ext)
        fig.savefig(path)
        out[ext.lstrip(".")] = str(path)
    return out
