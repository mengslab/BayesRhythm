mods=["streamlit","pandas","numpy","scipy","matplotlib","bayesrhythm_core_models","bayesrhythm_dashboard_nature_methods_style","br_io.single_cell_importers","exports.profiling","br_io.bulk_auto_parser","br_io.mex_loader","br_io.validators","annotation.qc_engine","annotation.marker_annotation_engine","annotation.reference_annotation_engine","annotation.confidence_engine","annotation.consensus_engine","annotation.annotation_review","anchors.multi_anchor_manager","anchors.preflight_qc","references.reference_diagnostics","references.multi_anchor_reference_builder","trajectories.state_aware_deconvolution_engine","exports.export_bundle","plotting.figure_studio","decomposition.decomposition_refactor","rhythms.rhythm_workflow"]
failed=[]
for m in mods:
    try: __import__(m)
    except Exception as e: failed.append((m,repr(e)))
if failed:
    print("Startup diagnostics failed:")
    for m,e in failed: print(f" - {m}: {e}")
    raise SystemExit(1)
print("Startup diagnostics passed.")
