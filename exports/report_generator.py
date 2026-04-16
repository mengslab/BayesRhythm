from pathlib import Path

def generate_markdown_report(output_path, sections: dict, figure_paths: dict | None = None):
    output_path = Path(output_path)
    lines = ["# BayesRhythm Analysis Report", ""]
    if figure_paths:
        lines.append("## Figures")
        lines.append("")
        for name, path in figure_paths.items():
            lines.append(f"- **{name}**: `{path}`")
        lines.append("")
    for title, body in sections.items():
        lines.append(f"## {title}")
        lines.append("")
        lines.append(body if isinstance(body, str) else str(body))
        lines.append("")
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path
