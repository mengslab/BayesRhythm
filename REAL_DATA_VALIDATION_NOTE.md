# BayesRhythm v5.0.1 Four-Anchor Validation Note

## Single-cell anchors validated
All four anchor archives are present and recognized as valid 10x MEX tarballs containing:
- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

### Anchor dimensions
| Anchor | Cells | Genes | Nonzero entries |
|---|---:|---:|---:|
| Chow ZT8 | 13,137 | 114,860 | 61,121,606 |
| Chow ZT20 | 10,853 | 108,503 | 51,826,943 |
| TRF ZT8 | 12,007 | 105,151 | 48,955,622 |
| TRF ZT20 | 11,929 | 118,527 | 64,520,416 |

## Bulk files validated
The joined hypothalamus bulk Excel files are supported through the joined-bulk auto-parse route.

## Validation conclusion
The package lineage is now validated for a full four-anchor single-cell panel at the file-structure and dimension level:
- Chow ZT8
- Chow ZT20
- TRF ZT8
- TRF ZT20

This removes the previous blocker caused by missing ZT20 anchor uploads.

## Remaining reality check
This confirms:
- all four anchors are present
- all four anchors have valid MEX structure
- the software can target a real four-anchor workflow

This does not, by itself, prove that every downstream analysis will run perfectly on first pass under all real inputs. The remaining work is execution-level validation inside the app workflow:
- batch import all four anchors
- batch annotation / consensus
- reference building
- deconvolution
- rhythm analysis
- decomposition
- report / figure export


## Execution-pass fixes applied
Concrete issues found during true execution pass:
1. Anchor registration originally failed because repeated per-cell ZT values were not deduplicated before uniqueness checks.
2. Multi-anchor reference building needed a memory-safety cap on shared genes for very large four-anchor runs.
3. The uploaded joined bulk Excel file does not encode Chow/TRF per sample in the sample column names, so condition-aware comparison requires an external metadata overlay.

## Additional execution fix applied
4. Multi-anchor execution previously became unstable because multiple full sparse anchor matrices were retained simultaneously. The package now uses lazy archive-backed anchor loading for registered anchors.


## Stable execution patch
This revision also patches execution-path stability issues:
- batch anchor pages now reload archive-backed anchors on demand
- joined bulk parsing supports a metadata overlay when condition labels are not embedded in sample names
- reference building and state-aware deconvolution are directly exposed in the UI for real workflow execution


## Final candidate stability patch
This revision aligns the packaged UI with the stabilized workflow code so the reference-builder and state-aware deconvolution path can actually be run from the app.
