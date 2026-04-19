BayesRhythm

BayesRhythm is a unified computational framework for state-aware rhythmic analysis of transcriptomes, integrating single-cell and bulk RNA-seq data through a novel multi-anchor reference strategy.
Biological rhythms—ranging from circadian (24-hour) to ultradian (e.g., 12-hour) cycles—are fundamentally shaped by cellular composition, cell-state programs, and environmental conditions. However, existing approaches typically analyze bulk and single-cell data in isolation, and lack the ability to jointly resolve rhythmicity, composition, and state-dependent regulation.
BayesRhythm addresses this gap by introducing a multi-anchor, state-aware analytical pipeline that connects single-cell references with bulk transcriptomic observations to enable mechanistic interpretation of rhythmic gene regulation.

Key Concepts

* Multi-anchor reference framework Instead of relying on a single reference dataset, BayesRhythm integrates multiple condition- and time-specific single-cell datasets (“anchors”) to construct a more biologically faithful reference atlas.
* State-aware deconvolution Bulk transcriptomes are decomposed not only into cell-type proportions, but also into cell-state programs, enabling resolution of dynamic biological processes.
* Integrated rhythmic analysis The framework supports:
    * rhythmic gene discovery
    * cross-condition rhythm comparison
    * state-dependent rhythmic decomposition
* End-to-end workflow From raw input to publication outputs, BayesRhythm provides:
    * anchor validation and consensus building
    * reference atlas construction
    * deconvolution and rhythm modeling
    * figure generation and report export

What BayesRhythm Enables

BayesRhythm allows researchers to:
* Dissect cell-state–specific rhythmic gene programs
* Compare rhythmic regulation across conditions and timepoints
* Separate true transcriptional rhythms from composition-driven effects
* Integrate single-cell resolution with bulk-scale statistical power
* Generate publication-ready figures and reports directly from analysis workflows

Scope and Applications

BayesRhythm is designed for studies involving:
* Metabolic and physiological rhythms (e.g., circadian, ultradian)
* Nutrition and environmental perturbations (e.g., feeding regimens)
* Disease-associated transcriptomic remodeling
* Multi-condition and multi-timepoint RNA-seq datasets

Corresponding Author & Copyright (c) 2026 Dr. Huan Meng, huan.meng@gmail.com


# BayesRhythm v3.2.0

This build hardens the single-cell preparation layer using real uploaded Chow/TRF MEX archives.

## Real issues addressed
- direct `.tar.gz` MEX archive ingestion
- sparse-safe QC on large matrices
- sparse-safe marker annotation using only marker genes
- sparse-safe reference scoring on shared genes only
- explicit condition / ZT assignment after archive import
- preflight summaries for real anchor readiness

## Run
```bash
chmod +x run.sh
./run.sh
```


## Upload limit
This package sets Streamlit `server.maxUploadSize = 500`, so app uploads can reach 500 MB when the local environment allows it.


## Local path import
This package also supports importing large files by local path, which bypasses browser upload limits.

Examples:
- `/Users/yourname/Downloads/Chow_ZT8_filtered_feature_bc_matrix.tar.gz`
- `/Users/yourname/Downloads/TRF_ZT8_filtered_feature_bc_matrix.tar.gz`


## Automatic anchor metadata inference
For archive filenames containing patterns like `Chow_ZT8` or `TRF-ZT20`, the app now auto-infers:
- condition
- zeitgeber time

You can still override the inferred values manually before import.


## Batch multi-archive import
This package adds a batch import mode for multiple local `.tar.gz` MEX archives.

Paste one local file path per line, for example:
- `/Users/yourname/Downloads/Chow_ZT8_filtered_feature_bc_matrix.tar.gz`
- `/Users/yourname/Downloads/Chow_ZT20_filtered_feature_bc_matrix.tar.gz`
- `/Users/yourname/Downloads/TRF_ZT8_filtered_feature_bc_matrix.tar.gz`
- `/Users/yourname/Downloads/TRF_ZT20_filtered_feature_bc_matrix.tar.gz`

The app will:
- infer condition and ZT from each filename
- import each archive
- validate each dataset
- auto-register each one into the anchor panel


## Batch import helper fix
This package fixes a missing helper-function bug in batch multi-archive import.
If you saw a `NameError` involving filename inference, this build resolves it.


## Anchor preflight QC + gating
This package adds an anchor preflight QC layer that classifies an imported anchor as:
- ready
- warning
- blocked

based on:
- validation status
- cell count
- median detected genes per cell
- unresolved label fraction
- low-confidence annotation fraction
- marker/reference disagreement fraction


## Anchor panel summary
This package adds a panel-wide preflight summary page for all registered anchors.

It shows:
- ready anchors
- warning anchors
- blocked anchors
- a table of anchor-level readiness and annotation risk fields


## Batch anchor actions
This package adds a batch anchor execution page that can:
- run marker annotation across all registered anchors
- recompute preflight readiness for all anchors
- refresh the panel-wide summary in one pass


## Batch reference annotation + consensus labeling
This package adds a batch page that can:
- run reference annotation across all registered anchors
- combine marker and reference outputs
- compute consensus labels
- apply consensus labels back into each anchor
- recompute readiness and panel summary


## Finish pass additions
This package adds:
- condition-aware rhythm comparison
- figure export to PNG/PDF/SVG
- decomposition that can use rhythm support
- smoke-test files for core modules

It is closer to a finished end-to-end research app, though still conservative rather than fully optimized.


## Completion pass additions
This package adds:
- deeper state-aware deconvolution
- more advanced rhythm fitting
- stronger decomposition tied to state programs
- richer publication-style figure layouts
- a more realistic smoke/regression test set
- more defensive data validation for future single-cell and bulk inputs


## Publication release additions
This package adds:
- direct auto-parsing of the uploaded joined bulk Excel format
- simple profiling hooks for major workflow steps
- fuller smoke/regression tests across connected workflows
- a real-data validation note for the files actually present during packaging

This is the final publication-release candidate in this package lineage.


## Four-anchor validation update
This package revision reflects validation of all four single-cell anchors:
- Chow ZT8
- Chow ZT20
- TRF ZT8
- TRF ZT20


## Execution-pass fixes
This package patches concrete failures exposed by real four-anchor execution:
- anchor registration now correctly deduplicates repeated per-cell zeitgeber values
- reference building caps shared genes for memory safety on large four-anchor runs
- joined bulk Excel parsing now supports an optional metadata overlay when sample columns lack explicit condition labels

## Lazy archive-backed anchor loading
This package avoids holding all four giant anchor matrices in memory at once. Anchors registered from archive files are now reloaded on demand for batch actions and reference building.


## Stable execution workflow patch
This package stabilizes the real four-anchor execution path by:
- using lazy anchor loading in batch marker and batch reference/consensus actions
- exposing reference building and state-aware deconvolution pages in the app
- supporting an optional bulk metadata overlay for joined bulk Excel files
- exporting reference and state-aware deconvolution outputs


## Final candidate stability patch
This build resolves a packaging inconsistency where stabilized workflow functions existed in code but were not exposed through the visible UI navigation and routes.
It also includes:
- optional bulk metadata overlay support in the visible Bulk Input page
- visible Reference Builder and State-Aware Deconvolution pages
- startup checks aligned to the packaged modules
- smoke tests for overlay merge and unique-ZT anchor registration


## Source-file normalization fix
This build fixes lazy anchor reloading when earlier records stored a comma-joined source path string containing both the archive and extracted temporary files.
Anchor records now keep a clean archive path, and old malformed records are normalized on reload.


## Downstream finish patch
This build adds and exposes:
- rhythm discovery
- rhythm comparison
- biological decomposition refactor
- figure studio
- report generator

The package now includes a clearer downstream analysis and reporting flow inside the app.
