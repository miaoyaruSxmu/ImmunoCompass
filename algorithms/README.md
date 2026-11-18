# ImmunoCompass analysis algorithms

This folder contains the key analysis scripts used in the ImmunoCompass
project. The installable R package is stored separately in the
`ImmunoCompass` package directory. The scripts here preserve the analysis
workflow and are provided for reproducibility, audit and future extension.

## Contents

- `01_IBI_calculation`: fixed-panel feature construction, IBI prediction and
  development-cohort validation.
- `02_feature_selection`: gene-pair screening, Top390 ranking and repeated
  feature-selection stability analyses.
- `03_model_comparison`: nested or cohort-grouped algorithm comparisons and
  metric summaries. `train_final_svr_from_cohort_table.R` is the authoritative
  package-model build script; it verifies the development allocation against
  Table S2 before fitting the locked 390-pair SVR.
- `04_BAS_CAS_IAS`: standardized gene-effect decomposition, dominance and
  partial-R2 analyses.
- `05_external_validation`: GSE65391, immunodeficiency and irAE validation
  workflows.
- `06_single_cell`: single-cell sparsity, dropout and pseudobulk sensitivity
  analyses.
- `07_cross_platform`: absolute-expression, continuous log-ratio and binary
  gene-pair representation comparisons.

## Reproducibility notes

Most scripts use project-specific input and output paths. Before running a
script, inspect and update its `root`, `tool_root`, `input_dir` or `output_dir`
variables. Scripts that use the `IBI_TOOL_ROOT` environment variable can be
run after setting that variable to the local project root.

The analysis scripts require R packages beyond the base R installation,
including `e1071`, `randomForest`, `xgboost`, `limma`, `pROC`, `ggplot2`,
`data.table`, `GSVA`, `msigdbr` and other packages documented in the resource
table. The scripts are not intended to download restricted data automatically.

The package implementation provides the user-facing functions:
`build_pair_features()`, `calculate_ibi()` and `decompose_gene_effect()`.

## Current bundled model

The current package uses the latest Top390 panel from the 781-sample
development configuration, matching Table S6. The final SVR model uses 807
samples from the 12 cohort-level
allocations marked `Model development` in Table S2. The class counts are 179
immunodeficiency or immunodeficiency-associated, 285 healthy and 343
autoimmune samples. The 390-pair feature panel remains fixed; GSE122459 adds
20 SLE and 6 healthy samples to SVR fitting. Run the build script with:

```text
Rscript 03_model_comparison/train_final_svr_from_cohort_table.R \
  fixed_feature_data.rda GSE122459.xlsx TableS2.xlsx \
  output_dir package_data_dir package_extdata_dir
```

Some scripts reproduce historical sensitivity analyses based on earlier
development allocations and are retained for manuscript auditing. They do not
rebuild the current bundled model. New-data scoring should use the installed R
package or `01_IBI_calculation/new_data_prediction.R`.
