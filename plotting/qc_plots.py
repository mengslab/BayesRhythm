import matplotlib.pyplot as plt
from plotting.publication_theme import apply_publication_theme
def plot_histogram(df, value_col, title):
    if df.empty or value_col not in df.columns: return None
    apply_publication_theme(); fig,ax=plt.subplots(figsize=(7,4.5)); ax.hist(df[value_col].astype(float), bins=25); ax.set_title(title); ax.set_xlabel(value_col); ax.set_ylabel("Count"); return fig
