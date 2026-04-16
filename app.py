import tempfile
import re
from pathlib import Path
import streamlit as st
import pandas as pd
from bayesrhythm_core_models import AnalysisMode, DatasetConfig, ReferenceStrategy, RhythmMode
from bayesrhythm_dashboard_nature_methods_style import render_nature_methods_dashboard
from br_io.single_cell_importers import SingleCellImporter
from br_io.validators import validate_single_cell_input, validate_anchor_metadata, validate_bulk_input
from br_io.bulk_auto_parser import parse_joined_bulk_excel, summarize_bulk_metadata, merge_bulk_metadata_overlay
from exports.profiling import StepProfiler
from annotation.qc_engine import compute_single_cell_qc
from annotation.marker_library import MARKER_LIBRARY
from annotation.marker_annotation_engine import score_markers_per_cell, assign_marker_labels
from annotation.reference_annotation_engine import build_reference_centroids, score_cells_against_reference, assign_reference_labels
from annotation.confidence_engine import combine_annotation_confidence
from annotation.consensus_engine import build_consensus_labels, summarize_consensus
from annotation.annotation_review import merge_annotation_sources, initialize_annotation_review, render_annotation_review_panel, apply_reviewed_labels
from anchors.multi_anchor_manager import MultiAnchorManager
from anchors.preflight_qc import compute_anchor_preflight, compute_panel_preflight_summary, compute_panel_preflight_status_counts
from references.reference_diagnostics import compute_reference_diagnostics
from references.multi_anchor_reference_builder import MultiAnchorReferenceBuilder
from trajectories.state_aware_deconvolution_engine import StateAwareDeconvolutionEngine
from rhythms.rhythm_discovery_engine import RhythmDiscoveryEngine
from rhythms.rhythm_comparison_engine import RhythmComparisonEngine
from decomposition.biological_decomposition_engine import BiologicalDecompositionEngine
from decomposition.state_program_decomposition_engine import StateProgramDecompositionEngine
from decomposition.decomposition_refactor import build_unified_decomposition_table
from rhythms.rhythm_workflow import summarize_rhythm_results
from exports.report_generator import generate_markdown_report
from plotting.figure_studio import make_bar_figure, make_two_panel, save_figure_bundle
import matplotlib.pyplot as plt
from exports.export_bundle import export_results_bundle
from plotting.qc_plots import plot_histogram

APP_TITLE="BayesRhythm v3.2.0"
NAV=["Dashboard","Single-Cell Input","Single-Cell QC","Annotation Engine","Annotation Review","Reference Diagnostics","Anchor Preflight QC","Anchor Panel Summary","Batch Anchor Actions","Batch Reference & Consensus","Anchor Panel Manager","Bulk Input","Reference Builder","State-Aware Deconvolution","Rhythm Discovery","Rhythm Comparison","Biological Decomposition","Figure Studio","Report Generator","Validation & Profiling","Export Results"]

def save_uploaded_file(uploaded_file):
    suffix="".join(Path(uploaded_file.name).suffixes) or Path(uploaded_file.name).suffix
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(uploaded_file.getbuffer()); return Path(tmp.name)

def read_table(path):
    path=Path(path)
    if path.suffix.lower()==".xlsx": return pd.read_excel(path)
    return pd.read_csv(path, sep="," if path.suffix.lower()==".csv" else "\t")


def infer_condition_zt_from_name(name):
    name = Path(str(name)).name
    lower = name.lower()
    condition = None
    if "chow" in lower:
        condition = "Chow"
    elif "trf" in lower:
        condition = "TRF"
    zt = None
    m = re.search(r'zt[_-]?(\d{1,2})', lower)
    if m:
        try:
            zt = int(m.group(1))
        except Exception:
            zt = None
    return condition, zt

def infer_dataset_id_from_name(name):
    condition, zt = infer_condition_zt_from_name(name)
    if condition is not None and zt is not None:
        return f"{condition}_ZT{int(zt)}"
    stem = Path(str(name)).name
    stem = re.sub(r"\.tar\.gz$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"[^A-Za-z0-9_\-]+", "_", stem)
    return stem[:80]

def init_state():
    defaults={"dataset_config":DatasetConfig("BayesRhythm v3.2 project","mouse",AnalysisMode.COMPARATIVE,["Chow","TRF"],RhythmMode.BOTH,(4.0,30.0),ReferenceStrategy.CONDITION_MATCHED),"page":"Dashboard","single_cell_import_result":None,"single_cell_validation_report":None,"single_cell_qc":None,"marker_assignments":None,"reference_assignments":None,"combined_annotation_confidence":None,"annotation_review_state":None,"reference_diagnostics":None,"multi_anchor_manager":MultiAnchorManager(),"bulk_expression_matrix":None,"bulk_sample_metadata":None,"bulk_validation_report":None,"anchor_annotation_results":{}, "anchor_consensus_summaries":{}, "reference_build_result":None, "state_deconvolution_run_result":None, "last_export_path":None}
    for k,v in defaults.items():
        if k not in st.session_state: st.session_state[k]=v

def readiness():
    mgr=st.session_state.multi_anchor_manager
    return {"SC imported":"complete" if st.session_state.single_cell_import_result is not None else "missing","QC":"complete" if st.session_state.single_cell_qc is not None else "missing","Annotation":"complete" if st.session_state.annotation_review_state is not None else "missing","Diagnostics":"complete" if st.session_state.reference_diagnostics is not None else "missing","Anchors":"complete" if len(mgr.get_anchor_panel().anchor_records)>0 else "missing","Bulk":"complete" if st.session_state.bulk_expression_matrix is not None else "missing"}

def render_dashboard():
    render_nature_methods_dashboard(st.session_state.dataset_config, readiness())
    st.write("This build hardens real-data execution for large Chow/TRF MEX archives.")


def render_single_cell_input():
    st.header("Single-Cell Input")
    importer = SingleCellImporter()
    mode = st.radio(
        "Import mode",
        ["Matrix + metadata", "10x MEX directory", "10x MEX tar.gz archive", "Batch multi-archive import"],
        horizontal=True,
    )

    if mode == "Matrix + metadata":
        tab_upload, tab_path = st.tabs(["Browser upload", "Local file paths"])
        with tab_upload:
            expr = st.file_uploader("Expression matrix (CSV/TSV/XLSX)", type=["csv", "tsv", "xlsx"])
            meta = st.file_uploader("Cell metadata (CSV/TSV/XLSX)", type=["csv", "tsv", "xlsx"], key="meta")
            if expr and meta and st.button("Import matrix bundle"):
                result = importer.import_bundle(save_uploaded_file(expr), save_uploaded_file(meta), dataset_metadata={"mode": "matrix_bundle"})
                st.session_state.single_cell_import_result = result
                st.session_state.single_cell_validation_report = validate_single_cell_input(result.single_cell)
                st.success("Imported matrix bundle.")
        with tab_path:
            expr_path = st.text_input("Local path to expression matrix", placeholder="/Users/yourname/Downloads/expression.csv")
            meta_path = st.text_input("Local path to cell metadata", placeholder="/Users/yourname/Downloads/metadata.csv")
            if st.button("Import matrix bundle from local paths"):
                p1, p2 = Path(expr_path).expanduser(), Path(meta_path).expanduser()
                if not p1.exists() or not p2.exists():
                    st.error("One or both local file paths do not exist.")
                else:
                    result = importer.import_bundle(p1, p2, dataset_metadata={"mode": "matrix_bundle_local_path"})
                    st.session_state.single_cell_import_result = result
                    st.session_state.single_cell_validation_report = validate_single_cell_input(result.single_cell)
                    st.success("Imported matrix bundle from local paths.")

    elif mode == "10x MEX directory":
        st.info("Use the archive mode for your uploaded `.tar.gz` files.")

    elif mode == "10x MEX tar.gz archive":
        tab_upload, tab_path = st.tabs(["Browser upload", "Local file path"])
        with tab_upload:
            archive = st.file_uploader("10x MEX `.tar.gz` archive", type=["gz"])
            inferred_condition, inferred_zt = infer_condition_zt_from_name(archive.name if archive is not None else "")
            if archive is not None:
                st.caption(f"Inferred from filename: condition={inferred_condition or 'unknown'}, ZT={inferred_zt if inferred_zt is not None else 'unknown'}")
            condition_default = (["Chow","TRF"].index(inferred_condition) if inferred_condition in ["Chow","TRF"] else 0)
            condition = st.selectbox("Condition to assign after import", ["Chow","TRF"], index=condition_default, key="cond_upload")
            zeitgeber_time = st.number_input("ZT to assign after import", min_value=0, max_value=24, value=(inferred_zt if inferred_zt is not None else 8), key="zt_upload")
            if archive and st.button("Import MEX archive"):
                path = save_uploaded_file(archive)
                result = importer.import_mex_archive(path, dataset_metadata={"mode": "mex_archive"})
                result.single_cell.cell_metadata["condition"] = condition
                result.single_cell.cell_metadata["zeitgeber_time"] = zeitgeber_time
                st.session_state.single_cell_import_result = result
                st.session_state.single_cell_validation_report = validate_single_cell_input(result.single_cell)
                st.success("Imported MEX archive and assigned anchor metadata.")
        with tab_path:
            archive_path = st.text_input("Local path to 10x MEX `.tar.gz` archive", placeholder="/Users/yourname/Downloads/Chow_ZT8_filtered_feature_bc_matrix.tar.gz")
            inferred_condition, inferred_zt = infer_condition_zt_from_name(archive_path)
            if archive_path:
                st.caption(f"Inferred from filename: condition={inferred_condition or 'unknown'}, ZT={inferred_zt if inferred_zt is not None else 'unknown'}")
            condition_default = (["Chow","TRF"].index(inferred_condition) if inferred_condition in ["Chow","TRF"] else 0)
            condition = st.selectbox("Condition to assign after import", ["Chow","TRF"], index=condition_default, key="cond_path")
            zeitgeber_time = st.number_input("ZT to assign after import", min_value=0, max_value=24, value=(inferred_zt if inferred_zt is not None else 8), key="zt_path")
            if st.button("Import MEX archive from local path"):
                p = Path(archive_path).expanduser()
                if not p.exists():
                    st.error("Local archive path does not exist.")
                else:
                    result = importer.import_mex_archive(p, dataset_metadata={"mode": "mex_archive_local_path"})
                    result.single_cell.cell_metadata["condition"] = condition
                    result.single_cell.cell_metadata["zeitgeber_time"] = zeitgeber_time
                    st.session_state.single_cell_import_result = result
                    st.session_state.single_cell_validation_report = validate_single_cell_input(result.single_cell)
                    st.success("Imported MEX archive from local path and assigned anchor metadata.")

    else:
        st.subheader("Batch multi-archive import")
        st.caption("Paste one local `.tar.gz` path per line. The app will infer Chow/TRF and ZT from each filename, import each archive, validate it, and auto-register it as an anchor.")
        paths_text = st.text_area(
            "Local archive paths",
            placeholder="/Users/yourname/Downloads/Chow_ZT8_filtered_feature_bc_matrix.tar.gz\n/Users/yourname/Downloads/TRF_ZT8_filtered_feature_bc_matrix.tar.gz",
            height=160,
        )
        override_unknown = st.checkbox("Stop if a filename does not clearly contain condition/ZT", value=True)
        if st.button("Batch import and auto-register anchors"):
            mgr = st.session_state.multi_anchor_manager
            lines = [line.strip() for line in paths_text.splitlines() if line.strip()]
            if not lines:
                st.error("Please provide at least one local archive path.")
            else:
                batch_rows = []
                any_error = False
                for line in lines:
                    p = Path(line).expanduser()
                    if not p.exists():
                        batch_rows.append({"path": str(p), "dataset_id": "", "condition": "", "zt": "", "status": "error", "message": "Path does not exist"})
                        any_error = True
                        continue
                    condition, zt = infer_condition_zt_from_name(p.name)
                    dataset_id = infer_dataset_id_from_name(p.name)
                    if override_unknown and (condition is None or zt is None):
                        batch_rows.append({"path": str(p), "dataset_id": dataset_id, "condition": condition or "", "zt": zt if zt is not None else "", "status": "error", "message": "Could not infer condition and/or ZT from filename"})
                        any_error = True
                        continue
                    try:
                        result = importer.import_mex_archive(p, dataset_metadata={"mode": "mex_archive_batch_local_path"})
                        if condition is not None:
                            result.single_cell.cell_metadata["condition"] = condition
                        if zt is not None:
                            result.single_cell.cell_metadata["zeitgeber_time"] = int(zt)
                        validation_report = validate_single_cell_input(result.single_cell)
                        mgr.register_anchor(
                            dataset_id=dataset_id,
                            sc_input=result.single_cell,
                            source_file=str(p),
                            import_mode="mex_archive_batch_local_path",
                            dataset_label=dataset_id,
                            validation_report=validation_report,
                            review_table=None,
                            keep_dataset=False,
                        )
                        batch_rows.append({
                            "path": str(p),
                            "dataset_id": dataset_id,
                            "condition": condition or "",
                            "zt": int(zt) if zt is not None else "",
                            "status": "imported",
                            "message": f"{result.single_cell.expression_matrix.shape[0]} cells × {result.single_cell.expression_matrix.shape[1]} genes",
                        })
                        st.session_state.single_cell_import_result = result
                        st.session_state.single_cell_validation_report = validation_report
                    except Exception as e:
                        batch_rows.append({"path": str(p), "dataset_id": dataset_id, "condition": condition or "", "zt": zt if zt is not None else "", "status": "error", "message": str(e)})
                        any_error = True
                st.session_state.multi_anchor_manager = mgr
                batch_df = pd.DataFrame(batch_rows)
                st.dataframe(batch_df, use_container_width=True)
                if any_error:
                    st.warning("Batch import finished with one or more errors.")
                else:
                    st.success("Batch import and auto-registration completed successfully.")

    if st.session_state.single_cell_import_result is not None:
        sc = st.session_state.single_cell_import_result.single_cell
        st.write(f"Imported shape: {sc.expression_matrix.shape[0]} cells × {sc.expression_matrix.shape[1]} genes")
        if st.session_state.single_cell_validation_report is not None:
            st.dataframe(st.session_state.single_cell_validation_report, use_container_width=True)

def render_single_cell_qc():
    st.header("Single-Cell QC"); result=st.session_state.single_cell_import_result
    if result is None: st.info("Import single-cell data first."); return
    if st.button("Compute sparse-safe QC"): st.session_state.single_cell_qc=compute_single_cell_qc(result.single_cell); st.success("QC computed.")
    qc=st.session_state.single_cell_qc
    if qc is not None:
        st.dataframe(qc["dataset_summary"], use_container_width=True); st.dataframe(qc["condition_counts"], use_container_width=True)
        fig=plot_histogram(qc["per_cell_qc"], "detected_genes", "Detected genes per cell")
        if fig is not None: st.pyplot(fig)

def render_annotation_engine():
    st.header("Annotation Engine"); result=st.session_state.single_cell_import_result
    if result is None: st.info("Import single-cell data first."); return
    if st.button("Run sparse-safe marker annotation"):
        scores=score_markers_per_cell(result.single_cell, MARKER_LIBRARY); st.session_state.marker_assignments=assign_marker_labels(scores); st.success("Marker annotation completed.")
    if st.session_state.marker_assignments is not None: st.dataframe(st.session_state.marker_assignments.head(50), use_container_width=True)
    ref_expr=st.file_uploader("Reference expression matrix", type=["csv","tsv","xlsx"], key="ref_expr"); ref_meta=st.file_uploader("Reference labels metadata", type=["csv","tsv","xlsx"], key="ref_meta")
    if ref_expr and ref_meta and st.button("Run sparse-safe reference annotation"):
        ref_expr_df=read_table(save_uploaded_file(ref_expr)); ref_meta_df=read_table(save_uploaded_file(ref_meta))
        if "cell_id" not in ref_expr_df.columns: ref_expr_df=ref_expr_df.rename(columns={ref_expr_df.columns[0]:"cell_id"})
        ref_expr_df["cell_id"]=ref_expr_df["cell_id"].astype(str); ref_expr_df=ref_expr_df.set_index("cell_id")
        if "cell_id" not in ref_meta_df.columns: ref_meta_df=ref_meta_df.rename(columns={ref_meta_df.columns[0]:"cell_id"})
        labels=ref_meta_df.set_index("cell_id")["cell_type"]; centroids=build_reference_centroids(ref_expr_df, labels); score_df=score_cells_against_reference(result.single_cell, centroids); st.session_state.reference_assignments=assign_reference_labels(score_df); st.success("Reference annotation completed.")
    if st.session_state.reference_assignments is not None: st.dataframe(st.session_state.reference_assignments.head(50), use_container_width=True)
    if st.button("Build combined annotation confidence"):
        st.session_state.combined_annotation_confidence=combine_annotation_confidence(st.session_state.marker_assignments, st.session_state.reference_assignments); review_table=merge_annotation_sources(result.single_cell.cell_metadata, st.session_state.marker_assignments, st.session_state.reference_assignments, st.session_state.combined_annotation_confidence); st.session_state.annotation_review_state=initialize_annotation_review(review_table); st.success("Combined annotation confidence built.")
    if st.session_state.combined_annotation_confidence is not None: st.dataframe(st.session_state.combined_annotation_confidence.head(50), use_container_width=True)

def render_annotation_review():
    st.header("Annotation Review")
    if st.session_state.annotation_review_state is None: st.info("Build combined annotation confidence first."); return
    st.session_state.annotation_review_state=render_annotation_review_panel(st.session_state.annotation_review_state)
    if st.button("Apply reviewed labels to dataset"):
        result=st.session_state.single_cell_import_result; result.single_cell.cell_metadata=apply_reviewed_labels(result.single_cell.cell_metadata, st.session_state.annotation_review_state); st.session_state.single_cell_import_result=result; st.success("Reviewed labels applied.")

def render_reference_diagnostics():
    st.header("Reference Diagnostics"); result=st.session_state.single_cell_import_result
    if result is None: st.info("Import single-cell data first."); return
    if st.button("Compute diagnostics"):
        st.session_state.reference_diagnostics=compute_reference_diagnostics(result.single_cell, st.session_state.marker_assignments, st.session_state.reference_assignments); st.success("Diagnostics computed.")
    if st.session_state.reference_diagnostics is not None:
        for name,df in st.session_state.reference_diagnostics.items(): st.markdown(f"**{name}**"); st.dataframe(df, use_container_width=True)


def render_anchor_preflight_qc():
    st.header("Anchor Preflight QC")
    result = st.session_state.single_cell_import_result
    review_state = st.session_state.annotation_review_state
    validation = st.session_state.single_cell_validation_report
    if result is None:
        st.info("Import single-cell data first.")
        return
    if st.button("Run anchor preflight QC"):
        review_table = review_state.review_table if review_state is not None else None
        preflight = compute_anchor_preflight(result.single_cell, review_table=review_table, validation_report=validation)
        st.session_state["anchor_preflight_table"] = preflight
        st.success("Anchor preflight QC computed.")
    if "anchor_preflight_table" in st.session_state and st.session_state["anchor_preflight_table"] is not None:
        st.dataframe(st.session_state["anchor_preflight_table"], use_container_width=True)



def run_marker_annotation_on_sc_input(sc_input):
    scores = score_markers_per_cell(sc_input, MARKER_LIBRARY)
    assigns = assign_marker_labels(scores)
    combined = combine_annotation_confidence(assigns, None)
    review_table = merge_annotation_sources(sc_input.cell_metadata, marker_df=assigns, reference_df=None, combined_conf_df=combined)
    preflight = compute_anchor_preflight(sc_input, review_table=review_table, validation_report=validate_single_cell_input(sc_input))
    return assigns, combined, review_table, preflight

def recompute_anchor_panel_summary(mgr):
    return compute_panel_preflight_summary(mgr)

def apply_consensus_labels_to_sc_input(sc_input, consensus_df):
    if consensus_df is None or consensus_df.empty:
        return sc_input
    labels = consensus_df[["cell_id", "consensus_cell_type"]].rename(columns={"consensus_cell_type":"cell_type"})
    meta = sc_input.cell_metadata.drop(columns=["cell_type"], errors="ignore").merge(labels, on="cell_id", how="left")
    meta["cell_type"] = meta["cell_type"].fillna("unknown")
    sc_input.cell_metadata = meta
    return sc_input


def render_anchor_panel_summary():
    st.header("Anchor Panel Summary")
    mgr = st.session_state.multi_anchor_manager
    summary_df = compute_panel_preflight_summary(mgr)
    if summary_df.empty:
        st.info("No anchors are registered yet.")
        return
    counts_df = compute_panel_preflight_status_counts(summary_df)
    c1, c2, c3 = st.columns(3)
    ready_n = int((summary_df["readiness_status"] == "ready").sum()) if "readiness_status" in summary_df.columns else 0
    warning_n = int((summary_df["readiness_status"] == "warning").sum()) if "readiness_status" in summary_df.columns else 0
    blocked_n = int((summary_df["readiness_status"] == "blocked").sum()) if "readiness_status" in summary_df.columns else 0
    c1.metric("Ready anchors", ready_n)
    c2.metric("Warning anchors", warning_n)
    c3.metric("Blocked anchors", blocked_n)
    st.dataframe(summary_df, use_container_width=True)
    if not counts_df.empty:
        st.markdown("**Preflight status counts**")
        st.dataframe(counts_df, use_container_width=True)
    st.session_state["anchor_panel_preflight_summary"] = summary_df


def render_batch_anchor_actions():
    st.header("Batch Anchor Actions")
    mgr = st.session_state.multi_anchor_manager
    summary = mgr.to_anchor_table()
    if summary.empty:
        st.info("No anchors are registered yet.")
        return
    st.dataframe(summary, use_container_width=True)

    if st.button("Run batch marker annotation + preflight recomputation"):
        batch_rows = []
        for record in mgr.get_anchor_panel().anchor_records:
            sc_input = mgr.get_dataset(record.dataset_id, cache=False)
            try:
                assigns, combined, review_table, preflight = run_marker_annotation_on_sc_input(sc_input)
                # apply final labels from marker-only review fallback
                sc_input.cell_metadata = apply_reviewed_labels(sc_input.cell_metadata, initialize_annotation_review(review_table))
                preflight_status = preflight.iloc[0]["preflight_status"] if not preflight.empty else "warning"
                unresolved_fraction = float(review_table["final_cell_type"].astype(str).str.lower().eq("unknown").mean()) if "final_cell_type" in review_table.columns else 1.0
                low_conf_fraction = float(combined["combined_annotation_confidence_label"].astype(str).str.lower().eq("low").mean()) if not combined.empty else 1.0
                disagreement_fraction = float((~combined["methods_agree"].fillna(False)).mean()) if not combined.empty and "methods_agree" in combined.columns else 1.0
                annotation_status = "reviewed" if unresolved_fraction <= 0.3 else "partially_resolved"
                mgr.update_anchor_dataset(
                    record.dataset_id,
                    sc_input=sc_input,
                    annotation_status=annotation_status,
                    readiness_status=preflight_status,
                    low_confidence_fraction=low_conf_fraction,
                    unresolved_fraction=unresolved_fraction,
                    disagreement_fraction=disagreement_fraction,
                )
                batch_rows.append({
                    "dataset_id": record.dataset_id,
                    "status": "processed",
                    "annotation_status": annotation_status,
                    "preflight_status": preflight_status,
                    "unresolved_fraction": unresolved_fraction,
                    "low_confidence_fraction": low_conf_fraction,
                    "disagreement_fraction": disagreement_fraction,
                })
            except Exception as e:
                batch_rows.append({
                    "dataset_id": record.dataset_id,
                    "status": "error",
                    "annotation_status": "",
                    "preflight_status": "blocked",
                    "unresolved_fraction": "",
                    "low_confidence_fraction": "",
                    "disagreement_fraction": "",
                    "message": str(e),
                })
        st.session_state.multi_anchor_manager = mgr
        batch_df = pd.DataFrame(batch_rows)
        st.session_state["batch_anchor_actions_table"] = batch_df
        st.session_state["anchor_panel_preflight_summary"] = recompute_anchor_panel_summary(mgr)
        st.success("Batch marker annotation and preflight recomputation completed.")
        st.dataframe(batch_df, use_container_width=True)

    if "batch_anchor_actions_table" in st.session_state and isinstance(st.session_state["batch_anchor_actions_table"], pd.DataFrame):
        st.markdown("**Last batch action summary**")
        st.dataframe(st.session_state["batch_anchor_actions_table"], use_container_width=True)


def render_batch_reference_consensus():
    st.header("Batch Reference & Consensus")
    mgr = st.session_state.multi_anchor_manager
    summary = mgr.to_anchor_table()
    if summary.empty:
        st.info("No anchors are registered yet.")
        return

    st.dataframe(summary, use_container_width=True)
    tab_upload, tab_path = st.tabs(["Browser upload", "Local file paths"])

    ref_expr_df = None
    ref_meta_df = None

    with tab_upload:
        ref_expr = st.file_uploader("Reference expression matrix (cells x genes)", type=["csv", "tsv", "xlsx"], key="batch_ref_expr")
        ref_meta = st.file_uploader("Reference labels metadata with cell_id and cell_type", type=["csv", "tsv", "xlsx"], key="batch_ref_meta")
        if ref_expr is not None and ref_meta is not None:
            ref_expr_df = read_table(save_uploaded_file(ref_expr))
            ref_meta_df = read_table(save_uploaded_file(ref_meta))

    with tab_path:
        ref_expr_path = st.text_input("Local path to reference expression matrix", placeholder="/Users/yourname/Downloads/reference_expr.csv")
        ref_meta_path = st.text_input("Local path to reference metadata", placeholder="/Users/yourname/Downloads/reference_meta.csv")
        if ref_expr_path and ref_meta_path:
            p1, p2 = Path(ref_expr_path).expanduser(), Path(ref_meta_path).expanduser()
            if p1.exists() and p2.exists():
                ref_expr_df = read_table(p1)
                ref_meta_df = read_table(p2)

    if ref_expr_df is None or ref_meta_df is None:
        st.info("Provide a reference expression matrix and metadata to run batch reference annotation and consensus labeling.")
        return

    if st.button("Run batch reference annotation + consensus relabeling"):
        if "cell_id" not in ref_expr_df.columns:
            ref_expr_df = ref_expr_df.rename(columns={ref_expr_df.columns[0]: "cell_id"})
        ref_expr_df["cell_id"] = ref_expr_df["cell_id"].astype(str)
        ref_expr_df = ref_expr_df.set_index("cell_id")

        if "cell_id" not in ref_meta_df.columns:
            ref_meta_df = ref_meta_df.rename(columns={ref_meta_df.columns[0]: "cell_id"})
        label_col = "cell_type" if "cell_type" in ref_meta_df.columns else ref_meta_df.columns[1]
        ref_labels = ref_meta_df.set_index("cell_id")[label_col]
        centroids = build_reference_centroids(ref_expr_df, ref_labels)

        batch_rows = []
        annotation_results = st.session_state.get("anchor_annotation_results", {})
        consensus_summaries = st.session_state.get("anchor_consensus_summaries", {})

        panel = mgr.get_anchor_panel()
        for record in panel.anchor_records:
            sc_input = mgr.get_dataset(record.dataset_id)
            try:
                marker_assigns = None
                if record.dataset_id in annotation_results and "marker" in annotation_results[record.dataset_id]:
                    marker_assigns = annotation_results[record.dataset_id]["marker"]
                else:
                    scores = score_markers_per_cell(sc_input, MARKER_LIBRARY)
                    marker_assigns = assign_marker_labels(scores)

                score_df = score_cells_against_reference(sc_input, centroids)
                reference_assigns = assign_reference_labels(score_df)
                combined = combine_annotation_confidence(marker_assigns, reference_assigns)
                consensus = build_consensus_labels(marker_assigns, reference_assigns)
                consensus_summary = summarize_consensus(consensus)

                sc_input = apply_consensus_labels_to_sc_input(sc_input, consensus)
                review_table = merge_annotation_sources(sc_input.cell_metadata, marker_assigns, reference_assigns, combined_conf_df=combined)
                validation = validate_single_cell_input(sc_input)
                preflight = compute_anchor_preflight(sc_input, review_table=review_table, validation_report=validation)

                preflight_status = preflight.iloc[0]["preflight_status"] if not preflight.empty else "warning"
                unresolved_fraction = float(consensus["consensus_cell_type"].astype(str).str.lower().eq("unknown").mean()) if not consensus.empty else 1.0
                low_conf_fraction = float(consensus["consensus_confidence_label"].astype(str).str.lower().eq("low").mean()) if not consensus.empty else 1.0
                disagreement_fraction = float((consensus["marker_predicted_label"].astype(str) != consensus["reference_predicted_label"].astype(str)).mean()) if not consensus.empty else 1.0
                annotation_status = "reviewed" if unresolved_fraction <= 0.3 else "partially_resolved"

                mgr.update_anchor_dataset(
                    record.dataset_id,
                    sc_input=sc_input,
                    annotation_status=annotation_status,
                    readiness_status=preflight_status,
                    low_confidence_fraction=low_conf_fraction,
                    unresolved_fraction=unresolved_fraction,
                    disagreement_fraction=disagreement_fraction,
                )

                annotation_results[record.dataset_id] = {
                    "marker": marker_assigns,
                    "reference": reference_assigns,
                    "combined": combined,
                    "consensus": consensus,
                }
                consensus_summaries[record.dataset_id] = consensus_summary

                batch_rows.append({
                    "dataset_id": record.dataset_id,
                    "status": "processed",
                    "preflight_status": preflight_status,
                    "annotation_status": annotation_status,
                    "unresolved_fraction": unresolved_fraction,
                    "low_confidence_fraction": low_conf_fraction,
                    "disagreement_fraction": disagreement_fraction,
                    "n_consensus_labels": int(consensus_summary["n_cells"].sum()) if not consensus_summary.empty else 0,
                })
            except Exception as e:
                batch_rows.append({
                    "dataset_id": record.dataset_id,
                    "status": "error",
                    "preflight_status": "blocked",
                    "annotation_status": "",
                    "unresolved_fraction": "",
                    "low_confidence_fraction": "",
                    "disagreement_fraction": "",
                    "n_consensus_labels": 0,
                    "message": str(e),
                })

        st.session_state.multi_anchor_manager = mgr
        st.session_state["anchor_annotation_results"] = annotation_results
        st.session_state["anchor_consensus_summaries"] = consensus_summaries
        st.session_state["batch_reference_consensus_table"] = pd.DataFrame(batch_rows)
        st.session_state["anchor_panel_preflight_summary"] = recompute_anchor_panel_summary(mgr)
        st.success("Batch reference annotation and consensus relabeling completed.")
        st.dataframe(st.session_state["batch_reference_consensus_table"], use_container_width=True)

    if "batch_reference_consensus_table" in st.session_state and isinstance(st.session_state["batch_reference_consensus_table"], pd.DataFrame):
        st.markdown("**Last batch reference + consensus summary**")
        st.dataframe(st.session_state["batch_reference_consensus_table"], use_container_width=True)

    consensus_summaries = st.session_state.get("anchor_consensus_summaries", {})
    if consensus_summaries:
        dataset_id = st.selectbox("View consensus summary for anchor", options=sorted(consensus_summaries.keys()))
        st.dataframe(consensus_summaries[dataset_id], use_container_width=True)


def render_anchor_panel_manager():
    st.header("Anchor Panel Manager"); mgr=st.session_state.multi_anchor_manager; result=st.session_state.single_cell_import_result; review_state=st.session_state.annotation_review_state; validation=st.session_state.single_cell_validation_report
    if result is not None:
        default_id=""; meta=result.single_cell.cell_metadata
        if {"condition","zeitgeber_time"}.issubset(meta.columns): summary=meta[["condition","zeitgeber_time"]].drop_duplicates(); 
        else: summary=pd.DataFrame()
        if not summary.empty and len(summary)==1: default_id=f"{summary.iloc[0]['condition']}_ZT{int(float(summary.iloc[0]['zeitgeber_time']))}"
        dataset_id=st.text_input("Dataset ID", value=default_id); dataset_label=st.text_input("Dataset label", value=default_id); st.dataframe(validate_anchor_metadata(result.single_cell.cell_metadata), use_container_width=True)
        if st.button("Register current dataset as anchor"):
            source_file=(result.single_cell.audit_trail.source_files[0] if result.single_cell.audit_trail and result.single_cell.audit_trail.source_files else "uploaded_bundle")
            try:
                rec=mgr.register_anchor(dataset_id, result.single_cell, source_file, result.single_cell.dataset_metadata.get("mode","upload"), dataset_label=dataset_label, validation_report=validation, review_table=(review_state.review_table if review_state is not None else None)); st.session_state.multi_anchor_manager=mgr; st.success(f"Registered anchor: {rec.dataset_id}")
            except Exception as e: st.error(str(e))
    v=mgr.validate_completeness(); st.dataframe(v.completeness_table, use_container_width=True); summary=mgr.to_anchor_table()
    if not summary.empty: st.dataframe(summary, use_container_width=True)
    for w in v.warnings: st.warning(w)


def render_bulk_input():
    st.header("Bulk Input")
    mode = st.radio("Bulk import mode", ["Matrix + metadata", "Joined bulk Excel auto-parse"], horizontal=True)
    if mode == "Matrix + metadata":
        expr = st.file_uploader("Bulk expression matrix", type=["csv", "tsv", "xlsx"])
        meta = st.file_uploader("Bulk metadata", type=["csv", "tsv", "xlsx"], key="bulk_meta")
        if expr and meta and st.button("Import bulk data"):
            profiler = StepProfiler()
            token = profiler.start("bulk_matrix_plus_metadata_import")
            expr_df = read_table(save_uploaded_file(expr))
            meta_df = read_table(save_uploaded_file(meta))
            gene_col = expr_df.columns[0]
            expr_df = expr_df.rename(columns={gene_col: "gene_id"}).set_index("gene_id")
            st.session_state.bulk_expression_matrix = expr_df
            st.session_state.bulk_sample_metadata = meta_df
            st.session_state.bulk_validation_report = validate_bulk_input(expr_df, meta_df)
            profiler.stop(token, {"n_genes": int(expr_df.shape[0]), "n_samples": int(expr_df.shape[1])})
            st.session_state.profiling_table = profiler.to_dataframe()
            st.success("Imported bulk data.")
    else:
        bulk_joined = st.file_uploader("Joined bulk Excel file", type=["xlsx"], key="bulk_joined")
        overlay = st.file_uploader("Optional bulk metadata overlay (CSV/TSV/XLSX with sample_id or accession, plus condition and zeitgeber_time)", type=["csv","tsv","xlsx"], key="bulk_overlay")
        if bulk_joined and st.button("Auto-parse joined bulk Excel"):
            profiler = StepProfiler()
            token = profiler.start("bulk_joined_excel_autoparse")
            df = read_table(save_uploaded_file(bulk_joined))
            expr_df, meta_df = parse_joined_bulk_excel(df)
            if overlay is not None:
                overlay_df = read_table(save_uploaded_file(overlay))
                meta_df = merge_bulk_metadata_overlay(meta_df, overlay_df)
            st.session_state.bulk_expression_matrix = expr_df
            st.session_state.bulk_sample_metadata = meta_df
            st.session_state.bulk_validation_report = validate_bulk_input(expr_df, meta_df)
            profiler.stop(token, {"n_genes": int(expr_df.shape[0]), "n_samples": int(expr_df.shape[1])})
            st.session_state.profiling_table = profiler.to_dataframe()
            st.success("Parsed joined bulk Excel file.")
            st.markdown("**Auto-parsed bulk metadata summary**")
            st.dataframe(summarize_bulk_metadata(meta_df), use_container_width=True)
            if "condition" in meta_df.columns and meta_df["condition"].astype(str).eq("unknown").all():
                st.warning("Condition labels could not be inferred from the bulk file alone. Upload a metadata overlay to enable condition-aware downstream comparisons.")
    if st.session_state.bulk_validation_report is not None:
        st.dataframe(st.session_state.bulk_validation_report, use_container_width=True)


def render_reference_builder():
    st.header("Reference Builder")
    mgr = st.session_state.multi_anchor_manager
    summary = mgr.to_anchor_table()
    if summary.empty:
        st.info("No anchors are registered yet.")
        return
    st.dataframe(summary, use_container_width=True)
    if st.button("Build multi-anchor reference atlas"):
        profiler = StepProfiler()
        token = profiler.start("reference_builder")
        result = MultiAnchorReferenceBuilder().run(mgr)
        profiler.stop(token, {
            "n_anchors": int(summary.shape[0]),
            "n_reference_rows": int(result.atlas.reference_table.shape[0]) if result.atlas.reference_table is not None else 0,
            "n_shared_genes": int(len(result.atlas.shared_genes)),
        })
        st.session_state.profiling_table = profiler.to_dataframe()
        st.session_state.reference_build_result = result
        st.success("Reference atlas built.")
    if st.session_state.reference_build_result is not None:
        for warning in st.session_state.reference_build_result.warnings:
            st.warning(warning)
        st.markdown("**Reference support table**")
        st.dataframe(st.session_state.reference_build_result.previews["support_table"], use_container_width=True)
        st.markdown("**Reference preview**")
        st.dataframe(st.session_state.reference_build_result.previews["reference_preview"], use_container_width=True)

def render_state_aware_deconvolution():
    st.header("State-Aware Deconvolution")
    bulk = st.session_state.bulk_expression_matrix
    bulk_md = st.session_state.bulk_sample_metadata
    ref_result = st.session_state.reference_build_result
    if ref_result is None:
        st.info("Build the reference atlas first.")
        return
    if bulk is None or bulk_md is None:
        st.info("Import bulk data first.")
        return
    st.markdown("**Bulk metadata summary**")
    st.dataframe(summarize_bulk_metadata(bulk_md), use_container_width=True)
    if st.button("Run state-aware deconvolution"):
        profiler = StepProfiler()
        token = profiler.start("state_aware_deconvolution")
        result = StateAwareDeconvolutionEngine().run(bulk, bulk_md, ref_result.atlas)
        profiler.stop(token, {
            "n_bulk_samples": int(bulk.shape[1]),
            "n_bulk_genes": int(bulk.shape[0]),
            "n_reference_rows": int(ref_result.atlas.reference_table.shape[0]),
            "n_shared_genes": int(len(ref_result.atlas.shared_genes)),
        })
        st.session_state.profiling_table = profiler.to_dataframe()
        st.session_state.state_deconvolution_run_result = result
        st.success("State-aware deconvolution completed.")
    if st.session_state.state_deconvolution_run_result is not None:
        for warning in st.session_state.state_deconvolution_run_result.warnings:
            st.warning(warning)
        st.markdown("**Fit metrics**")
        st.dataframe(st.session_state.state_deconvolution_run_result.result.fit_metrics, use_container_width=True)
        st.markdown("**State program scores**")
        st.dataframe(st.session_state.state_deconvolution_run_result.result.state_program_scores, use_container_width=True)
        st.markdown("**Cell fractions preview**")
        st.dataframe(st.session_state.state_deconvolution_run_result.result.cell_fractions.head(20), use_container_width=True)



def render_rhythm_discovery():
    st.header("Rhythm Discovery")
    bulk = st.session_state.bulk_expression_matrix
    bulk_md = st.session_state.bulk_sample_metadata
    if bulk is None or bulk_md is None:
        st.info("Import bulk data first.")
        return
    if st.button("Run rhythm discovery"):
        profiler = StepProfiler()
        token = profiler.start("rhythm_discovery")
        result = RhythmDiscoveryEngine().run(bulk, bulk_md)
        profiler.stop(token, {"n_bulk_samples": int(bulk.shape[1]), "n_bulk_genes": int(bulk.shape[0])})
        st.session_state.profiling_table = profiler.to_dataframe()
        st.session_state.rhythm_run_result = result
        st.success("Rhythm discovery completed.")
    if st.session_state.rhythm_run_result is not None:
        st.dataframe(st.session_state.rhythm_run_result.result.gene_results.head(50), use_container_width=True)
        summary = summarize_rhythm_results(simple_df=st.session_state.rhythm_run_result.result.gene_results)
        st.markdown("**Rhythm summary**")
        st.dataframe(summary.summary_table, use_container_width=True)

def render_rhythm_comparison():
    st.header("Rhythm Comparison")
    bulk = st.session_state.bulk_expression_matrix
    bulk_md = st.session_state.bulk_sample_metadata
    if bulk is None or bulk_md is None:
        st.info("Import bulk data first.")
        return
    if st.button("Run rhythm comparison"):
        profiler = StepProfiler()
        token = profiler.start("rhythm_comparison")
        result = RhythmComparisonEngine().run(bulk, bulk_md)
        profiler.stop(token, {"n_bulk_samples": int(bulk.shape[1]), "n_bulk_genes": int(bulk.shape[0])})
        st.session_state.profiling_table = profiler.to_dataframe()
        st.session_state.rhythm_comparison_run_result = result
        st.success("Rhythm comparison completed.")
    if st.session_state.rhythm_comparison_run_result is not None:
        st.dataframe(st.session_state.rhythm_comparison_run_result.result.gene_results.head(50), use_container_width=True)
        summary = summarize_rhythm_results(comparison_df=st.session_state.rhythm_comparison_run_result.result.gene_results)
        st.markdown("**Comparison summary**")
        st.dataframe(summary.summary_table, use_container_width=True)

def render_biological_decomposition():
    st.header("Biological Decomposition")
    bulk = st.session_state.bulk_expression_matrix
    state_deconv = st.session_state.state_deconvolution_run_result
    if bulk is None or state_deconv is None:
        st.info("Run state-aware deconvolution first.")
        return
    rhythm_df = st.session_state.rhythm_run_result.result.gene_results if st.session_state.rhythm_run_result is not None else None
    if st.button("Run biological decomposition"):
        profiler = StepProfiler()
        token = profiler.start("biological_decomposition")
        basic = BiologicalDecompositionEngine().run(bulk, state_deconv.result, rhythm_results=rhythm_df)
        state_prog = StateProgramDecompositionEngine().run(bulk, state_deconv.result, rhythm_results=rhythm_df)
        unified = build_unified_decomposition_table(
            basic.result.gene_results if basic is not None else None,
            state_prog.result.gene_results if state_prog is not None else None,
        )
        profiler.stop(token, {"n_bulk_samples": int(bulk.shape[1]), "n_bulk_genes": int(bulk.shape[0])})
        st.session_state.profiling_table = profiler.to_dataframe()
        st.session_state.decomposition_run_result = basic
        st.session_state.state_program_decomposition_run_result = state_prog
        st.session_state.unified_decomposition_result = unified
        st.success("Biological decomposition completed.")
    if st.session_state.unified_decomposition_result is not None:
        st.markdown("**Unified decomposition table**")
        st.dataframe(st.session_state.unified_decomposition_result.unified_gene_table.head(50), use_container_width=True)

def render_figure_studio():
    st.header("Figure Studio")
    choices = []
    if st.session_state.rhythm_run_result is not None:
        choices.append("rhythm_fit_r2")
    if st.session_state.rhythm_comparison_run_result is not None:
        choices.append("rhythm_comparison_delta")
    if st.session_state.unified_decomposition_result is not None and not st.session_state.unified_decomposition_result.unified_gene_table.empty:
        choices.append("decomposition_two_panel")
    if not choices:
        st.info("Run rhythm and/or decomposition analyses first.")
        return
    choice = st.selectbox("Figure preset", options=choices)
    if st.button("Create figure bundle"):
        fig = None
        if choice == "rhythm_fit_r2":
            df = st.session_state.rhythm_run_result.result.gene_results.head(30)
            fig = make_bar_figure(df, "gene_id", "fit_r2", "Top rhythm fit R2", rotation=60)
        elif choice == "rhythm_comparison_delta":
            df = st.session_state.rhythm_comparison_run_result.result.gene_results.head(20)
            fig = make_bar_figure(df, "gene_id", "amplitude_delta", "Rhythm comparison amplitude delta", rotation=60)
        elif choice == "decomposition_two_panel":
            left = st.session_state.rhythm_comparison_run_result.result.gene_results.head(15) if st.session_state.rhythm_comparison_run_result is not None else None
            right = st.session_state.unified_decomposition_result.unified_gene_table.head(15)
            if left is not None and "amplitude_delta" in left.columns:
                y2 = "composition_score" if "composition_score" in right.columns else ("state_program_r2" if "state_program_r2" in right.columns else right.select_dtypes(include="number").columns[0])
                fig = make_two_panel(left, "gene_id", "amplitude_delta", "Rhythm amplitude delta", right, "gene_id", y2, "Decomposition score")
        if fig is not None:
            st.session_state.last_figure_paths = save_figure_bundle(fig, str(Path.cwd() / "exports" / choice))
            st.success("Figure bundle created.")
    if st.session_state.last_figure_paths:
        st.json(st.session_state.last_figure_paths)

def render_report_generator():
    st.header("Report Generator")
    if st.button("Generate markdown report"):
        sections = {}
        if "anchor_panel_preflight_summary" in st.session_state and st.session_state["anchor_panel_preflight_summary"] is not None:
            sections["Anchor Panel Summary"] = st.session_state["anchor_panel_preflight_summary"].head(20).to_markdown(index=False)
        if st.session_state.reference_build_result is not None:
            sections["Reference Support"] = st.session_state.reference_build_result.previews["support_table"].head(20).to_markdown(index=False)
        if st.session_state.state_deconvolution_run_result is not None:
            sections["State-Aware Deconvolution"] = st.session_state.state_deconvolution_run_result.result.fit_metrics.head(20).to_markdown(index=False)
        if st.session_state.rhythm_run_result is not None:
            sections["Rhythm Discovery"] = st.session_state.rhythm_run_result.result.gene_results.head(20).to_markdown(index=False)
        if st.session_state.rhythm_comparison_run_result is not None:
            sections["Rhythm Comparison"] = st.session_state.rhythm_comparison_run_result.result.gene_results.head(20).to_markdown(index=False)
        if st.session_state.unified_decomposition_result is not None and not st.session_state.unified_decomposition_result.unified_gene_table.empty:
            sections["Unified Decomposition"] = st.session_state.unified_decomposition_result.unified_gene_table.head(20).to_markdown(index=False)
        if not sections:
            st.info("Run one or more downstream analyses first.")
        else:
            figs = st.session_state.last_figure_paths if st.session_state.last_figure_paths is not None else None
            out = generate_markdown_report(Path.cwd() / "exports" / "bayesrhythm_report.md", sections, figure_paths=figs)
            st.session_state.last_report_path = str(out)
            st.success("Report generated.")
    if st.session_state.last_report_path:
        st.code(st.session_state.last_report_path)


def render_validation_profiling():
    st.header("Validation & Profiling")
    note_path = Path.cwd() / "REAL_DATA_VALIDATION_NOTE.md"
    if note_path.exists():
        st.markdown("**Real-data validation note**")
        st.code(note_path.read_text(encoding="utf-8"))
    st.markdown("**Execution-path notes**")
    st.write("This build exposes the stabilized reference-builder and state-aware deconvolution workflow in the UI and supports optional bulk metadata overlays for joined bulk files.")
    st.markdown("**Execution-pass fixes included in this build**")
    st.write("This build uses lazy archive-backed anchor loading in batch pages, supports optional bulk metadata overlays, and exposes the reference-builder and state-aware deconvolution workflow directly in the app.")
    if st.session_state.profiling_table is not None:
        st.markdown("**Latest profiling table**")
        st.dataframe(st.session_state.profiling_table, use_container_width=True)
    else:
        st.info("No profiling data has been recorded yet. Import bulk data or run workflow steps that collect timing information.")

def render_export_results():
    st.header("Export Results"); tables={}
    for name,key in [("single_cell_validation","single_cell_validation_report"),("marker_assignments","marker_assignments"),("reference_assignments","reference_assignments"),("combined_annotation_confidence","combined_annotation_confidence"),("bulk_validation","bulk_validation_report")]:
        val=st.session_state.get(key)
        if isinstance(val,pd.DataFrame) and not val.empty: tables[name]=val
    if st.session_state.single_cell_qc is not None:
        for k,df in st.session_state.single_cell_qc.items(): tables[f"qc_{k}"]=df
    if st.session_state.reference_diagnostics is not None:
        for k,df in st.session_state.reference_diagnostics.items():
            if isinstance(df,pd.DataFrame) and not df.empty: tables[f"diag_{k}"]=df
    if "anchor_preflight_table" in st.session_state and isinstance(st.session_state["anchor_preflight_table"], pd.DataFrame) and not st.session_state["anchor_preflight_table"].empty:
        tables["anchor_preflight_qc"] = st.session_state["anchor_preflight_table"]
    if "anchor_panel_preflight_summary" in st.session_state and isinstance(st.session_state["anchor_panel_preflight_summary"], pd.DataFrame) and not st.session_state["anchor_panel_preflight_summary"].empty:
        tables["anchor_panel_preflight_summary"] = st.session_state["anchor_panel_preflight_summary"]
    if "batch_anchor_actions_table" in st.session_state and isinstance(st.session_state["batch_anchor_actions_table"], pd.DataFrame) and not st.session_state["batch_anchor_actions_table"].empty:
        tables["batch_anchor_actions_table"] = st.session_state["batch_anchor_actions_table"]
    if "batch_reference_consensus_table" in st.session_state and isinstance(st.session_state["batch_reference_consensus_table"], pd.DataFrame) and not st.session_state["batch_reference_consensus_table"].empty:
        tables["batch_reference_consensus_table"] = st.session_state["batch_reference_consensus_table"]
    if st.session_state.reference_build_result is not None:
        tables["reference_support_table"] = st.session_state.reference_build_result.previews["support_table"]
    if st.session_state.state_deconvolution_run_result is not None:
        tables["state_deconvolution_fit_metrics"] = st.session_state.state_deconvolution_run_result.result.fit_metrics
        tables["state_program_scores"] = st.session_state.state_deconvolution_run_result.result.state_program_scores
        tables["state_deconvolution_fractions"] = st.session_state.state_deconvolution_run_result.result.cell_fractions
    if st.session_state.rhythm_run_result is not None:
        tables["rhythm_gene_results"] = st.session_state.rhythm_run_result.result.gene_results
    if st.session_state.rhythm_comparison_run_result is not None:
        tables["rhythm_comparison_results"] = st.session_state.rhythm_comparison_run_result.result.gene_results
    if st.session_state.unified_decomposition_result is not None and not st.session_state.unified_decomposition_result.unified_gene_table.empty:
        tables["unified_decomposition_results"] = st.session_state.unified_decomposition_result.unified_gene_table
    summary=st.session_state.multi_anchor_manager.to_anchor_table()
    if not summary.empty: tables["anchor_panel_summary"]=summary
    if not tables: st.info("No exportable tables available yet."); return
    if st.button("Create export bundle"): export_path=export_results_bundle(Path.cwd()/"exports", tables); st.session_state.last_export_path=str(export_path); st.success(f"Created export bundle: {export_path.name}")
    if st.session_state.last_export_path: st.code(st.session_state.last_export_path)

def main():
    st.set_page_config(page_title=APP_TITLE, layout="wide"); init_state(); page=st.sidebar.radio("Navigate", NAV)
    if page=="Dashboard": render_dashboard()
    elif page=="Single-Cell Input": render_single_cell_input()
    elif page=="Single-Cell QC": render_single_cell_qc()
    elif page=="Annotation Engine": render_annotation_engine()
    elif page=="Annotation Review": render_annotation_review()
    elif page=="Reference Diagnostics": render_reference_diagnostics()
    elif page=="Anchor Preflight QC": render_anchor_preflight_qc()
    elif page=="Anchor Panel Summary": render_anchor_panel_summary()
    elif page=="Batch Anchor Actions": render_batch_anchor_actions()
    elif page=="Batch Reference & Consensus": render_batch_reference_consensus()
    elif page=="Anchor Panel Manager": render_anchor_panel_manager()
    elif page=="Bulk Input": render_bulk_input()
    elif page=="Reference Builder": render_reference_builder()
    elif page=="State-Aware Deconvolution": render_state_aware_deconvolution()
    elif page=="Rhythm Discovery": render_rhythm_discovery()
    elif page=="Rhythm Comparison": render_rhythm_comparison()
    elif page=="Biological Decomposition": render_biological_decomposition()
    elif page=="Figure Studio": render_figure_studio()
    elif page=="Report Generator": render_report_generator()
    elif page=="Validation & Profiling": render_validation_profiling()
    elif page=="Export Results": render_export_results()

if __name__=="__main__": main()
