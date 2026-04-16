from pathlib import Path
import matplotlib.pyplot as plt

def build_two_panel_bar_figure(df_left, x_left, y_left, title_left, df_right, x_right, y_right, title_right, output_base):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].bar(df_left[x_left].astype(str), df_left[y_left].astype(float))
    axes[0].set_title(title_left)
    axes[0].tick_params(axis='x', rotation=45)
    axes[1].bar(df_right[x_right].astype(str), df_right[y_right].astype(float))
    axes[1].set_title(title_right)
    axes[1].tick_params(axis='x', rotation=45)
    out = {}
    base = Path(output_base)
    base.parent.mkdir(parents=True, exist_ok=True)
    for ext in [".png", ".pdf", ".svg"]:
        path = base.with_suffix(ext)
        fig.savefig(path, bbox_inches="tight", dpi=300)
        out[ext.lstrip(".")] = str(path)
    plt.close(fig)
    return out
