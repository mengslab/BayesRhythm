from pathlib import Path
import matplotlib.pyplot as plt

def save_figure_bundle(fig, output_base):
    base = Path(output_base)
    base.parent.mkdir(parents=True, exist_ok=True)
    out = {}
    for ext in [".png", ".pdf", ".svg"]:
        p = base.with_suffix(ext)
        fig.savefig(p, bbox_inches="tight", dpi=300)
        out[ext.lstrip(".")] = str(p)
    plt.close(fig)
    return out

def make_bar_figure(df, x, y, title, rotation=45):
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(df[x].astype(str), df[y].astype(float))
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=rotation)
    return fig

def make_two_panel(df1, x1, y1, title1, df2, x2, y2, title2):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].bar(df1[x1].astype(str), df1[y1].astype(float))
    axes[0].set_title(title1)
    axes[0].tick_params(axis="x", rotation=45)
    axes[1].bar(df2[x2].astype(str), df2[y2].astype(float))
    axes[1].set_title(title2)
    axes[1].tick_params(axis="x", rotation=45)
    return fig
