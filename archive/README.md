# archive/

These files were part of an earlier version of the pipeline that used a
**per-plate normalization → merge** approach. The pipeline was subsequently
refactored to use **combined normalization** (all plates normalized together in
a single meffil functional normalization call), which supersedes this approach.

None of these files are referenced in the current `main_meffil.nf`.

| File | Reason archived |
|------|-----------------|
| `merge_plates.r` | Merged per-plate beta CSVs and applied ComBat; replaced by `combined_normalize_meffil.r` which handles all plates in one pass |
| `plate_normalize_meffil.r` | Per-plate normalization script; made redundant by combined normalization |
| `multi_plate_qc_report.r` | Cross-plate QC report; the same functionality is now embedded in Step 6 of `combined_normalize_meffil.r` |
| `diagnose_preprocesscore.r` | One-off diagnostic for the `pthread_create` preprocessCore bug; superseded by the fix documented in the README |
| `fix_meffil_threading.sh` | Workaround shell script for the threading bug; the fix (`--disable-threading` reinstall) is now documented in the README setup section and should be applied once during environment creation |

To restore any of these files: `cp archive/<file> ./`
