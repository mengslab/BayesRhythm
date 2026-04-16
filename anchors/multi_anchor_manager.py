import pandas as pd
from bayesrhythm_core_models import AnchorDatasetRecord, AnchorPanel, AnchorPanelValidationResult
from anchors.preflight_qc import compute_anchor_preflight

EXPECTED_ANCHORS=[("Chow",8.0),("Chow",20.0),("TRF",8.0),("TRF",20.0)]

class MultiAnchorManager:
    def _normalize_source_file(self, source_file):
        src = str(source_file)
        parts = [p.strip() for p in src.split(',') if p.strip()]
        for p in parts:
            low = p.lower()
            if low.endswith('.tar.gz') or low.endswith('.tgz'):
                return p
        return parts[0] if parts else src

    def __init__(self):
        self._datasets={}
        self._records={}
        self._dataset_sources={}

    def register_anchor(self, dataset_id, sc_input, source_file, import_mode, dataset_label=None, validation_report=None, review_table=None, keep_dataset=False):
        conds=sorted(sc_input.cell_metadata["condition"].astype(str).replace("", pd.NA).dropna().unique().tolist())
        zts=sorted(pd.to_numeric(sc_input.cell_metadata["zeitgeber_time"], errors="coerce").dropna().unique().tolist())
        if len(conds)!=1 or len(zts)!=1:
            raise ValueError(f"Each anchor dataset must contain exactly one condition and one zeitgeber time. Found conditions={conds}, zt={zts}")
        validation_status="passed"
        if validation_report is not None and not validation_report.empty:
            if (validation_report["severity"]=="error").any(): validation_status="failed"
            elif (validation_report["severity"]=="warning").any(): validation_status="warning"
        annotation_status="annotation_required"; readiness_status="not_ready"; low_conf=1.0; unresolved=1.0; disagreement=1.0
        preflight = compute_anchor_preflight(sc_input, review_table=review_table, validation_report=validation_report)
        preflight_status = preflight.iloc[0]["preflight_status"] if not preflight.empty else "warning"
        if review_table is not None and not review_table.empty:
            final=review_table.get("final_cell_type", pd.Series(dtype=str)).astype(str)
            unresolved=float(final.str.lower().eq("unknown").mean()) if len(final) else 1.0
            low_conf=float(review_table["combined_annotation_confidence_label"].astype(str).str.lower().eq("low").mean()) if "combined_annotation_confidence_label" in review_table.columns else 0.0
            disagreement=float((~review_table["methods_agree"].fillna(False)).mean()) if "methods_agree" in review_table.columns else 0.0
            annotation_status="reviewed" if unresolved<=0.3 else "partially_resolved"
        readiness_status = preflight_status
        record=AnchorDatasetRecord(dataset_id, conds[0], float(zts[0]), source_file, int(sc_input.expression_matrix.shape[0]), int(sc_input.expression_matrix.shape[1]), annotation_status, validation_status, readiness_status, import_mode, dataset_label or dataset_id, low_confidence_fraction=low_conf, unresolved_fraction=unresolved, disagreement_fraction=disagreement)
        self._records[dataset_id]=record
        if keep_dataset:
            self._datasets[dataset_id]=sc_input
        else:
            self._datasets[dataset_id]=None
            self._dataset_sources[dataset_id]={"source_file":source_file,"import_mode":import_mode,"condition":conds[0],"zeitgeber_time":float(zts[0])}
        return record

    def get_dataset(self, dataset_id, cache=False):
        if dataset_id in self._datasets and self._datasets[dataset_id] is not None:
            return self._datasets[dataset_id]
        if dataset_id not in self._dataset_sources:
            raise KeyError(f"No dataset or source registered for {dataset_id}")
        src=self._dataset_sources[dataset_id]
        from br_io.single_cell_importers import SingleCellImporter
        imp=SingleCellImporter()
        sc=imp.import_mex_archive(src["source_file"], dataset_metadata={"mode":src["import_mode"]}).single_cell
        sc.cell_metadata["condition"]=src["condition"]
        sc.cell_metadata["zeitgeber_time"]=src["zeitgeber_time"]
        if cache:
            self._datasets[dataset_id]=sc
        return sc

    def update_anchor_dataset(self, dataset_id, sc_input=None, annotation_status=None, readiness_status=None,
                             low_confidence_fraction=None, unresolved_fraction=None, disagreement_fraction=None,
                             keep_dataset=False):
        if dataset_id not in self._records:
            raise KeyError(f"Unknown anchor dataset_id: {dataset_id}")
        rec=self._records[dataset_id]
        if sc_input is not None:
            rec.n_cells=int(sc_input.expression_matrix.shape[0]); rec.n_genes=int(sc_input.expression_matrix.shape[1])
            if keep_dataset:
                self._datasets[dataset_id]=sc_input
        if annotation_status is not None: rec.annotation_status=annotation_status
        if readiness_status is not None: rec.readiness_status=readiness_status
        if low_confidence_fraction is not None: rec.low_confidence_fraction=float(low_confidence_fraction)
        if unresolved_fraction is not None: rec.unresolved_fraction=float(unresolved_fraction)
        if disagreement_fraction is not None: rec.disagreement_fraction=float(disagreement_fraction)
        self._records[dataset_id]=rec

    def remove_anchor(self, dataset_id):
        self._datasets.pop(dataset_id,None); self._records.pop(dataset_id,None); self._dataset_sources.pop(dataset_id,None)

    def get_anchor_panel(self):
        recs=sorted(self._records.values(), key=lambda r:(r.condition,r.zeitgeber_time,r.dataset_id))
        return AnchorPanel(recs, {}, {"n_registered_anchors":len(recs)})

    def validate_completeness(self):
        panel=self.get_anchor_panel(); observed={(r.condition,float(r.zeitgeber_time)):r.dataset_id for r in panel.anchor_records}; rows=[]; missing=[]; warnings=[]
        for c,z in EXPECTED_ANCHORS:
            key=f"{c}_ZT{int(z)}"; ds=observed.get((c,z)); present=ds is not None; rows.append({"anchor_key":key,"condition":c,"zeitgeber_time":z,"present":present,"dataset_id":ds})
            if not present: missing.append(key)
        summary=panel.to_summary_table()
        if not summary.empty:
            if (summary["validation_status"]=="failed").any(): warnings.append("One or more anchors failed validation.")
            if (summary["annotation_status"].isin(["annotation_required","partially_resolved"])).any(): warnings.append("One or more anchors still need annotation cleanup before reference construction.")
            if (summary["readiness_status"]!="ready").any(): warnings.append("One or more anchors are not yet ready for real downstream runs.")
        if missing: warnings.append("Anchor panel is incomplete. Missing anchors: " + ", ".join(missing))
        return AnchorPanelValidationResult(len(missing)==0, pd.DataFrame(rows), missing, warnings)

    def to_anchor_table(self):
        return self.get_anchor_panel().to_summary_table()
