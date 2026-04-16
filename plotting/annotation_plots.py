import matplotlib.pyplot as plt
from plotting.publication_theme import apply_publication_theme
def plot_label_counts(df, label_col, title):
    if df.empty or label_col not in df.columns: return None
    counts=df[label_col].astype(str).value_counts(); apply_publication_theme(); fig,ax=plt.subplots(figsize=(8,4.5)); ax.bar(counts.index.astype(str), counts.values.astype(float)); ax.set_title(title); ax.tick_params(axis="x", rotation=30); return fig
