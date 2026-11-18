# ImmunoCompass

ImmunoCompass provides a reproducible R implementation of the transcriptomic
Immune Balance Index (IBI) and immune-context-aware gene effect decomposition.
The package contains the latest 781-sample Top390 gene-pair panel reported in
Table S6 and a radial SVR model refitted to the 807 model-development samples specified in Table S2
(285 healthy, 343 autoimmune, and 179 immunodeficiency or
immunodeficiency-associated samples). The panel is fixed after selection in
the 781-sample configuration; no feature re-selection is performed when
scoring new data.

## Installation

Install the local package from its root directory with:

```r
install.packages(".", repos = NULL, type = "source")
```

The package depends on `e1071` for the bundled SVR model. To install the
development version from GitHub after uploading this directory:

```r
remotes::install_github("lydiaMyr/ImmunoCompass", subdir = "ImmunoCompass")
```

## Example

The bundled example contains six samples and the gene rows required by the
locked panel.

```r
library(ImmunoCompass)
data(example_expression)
data(example_metadata)

ibi <- calculate_ibi(example_expression)
head(ibi)

decomposition <- decompose_gene_effect(
  expression = example_expression,
  ibi = ibi,
  condition = example_metadata$condition
)
head(decomposition)
```

`build_pair_features()` exposes the 390 binary pair features. By default,
GeneA > GeneB is encoded as 1 and all other ties or reverse orderings as 0,
matching the bundled SVR object. Missing pair values use the pair-specific
binary ordering modes estimated from the 807 model-development samples. Zero
encoding and strict failure remain available through `missing = "zero"` and
`missing = "error"`, respectively.

Model provenance, training-cohort counts, class counts, panel provenance,
missing-pair modes and fitted hyperparameters are available in
`model_metadata` and `gene_pair_modes`. The bundled SVR uses
`cost = 10`, `gamma = 0.0001`, `epsilon = 0.1`, a radial kernel and internal
feature scaling.

The `algorithms/` directory distributed with the project contains the key
feature-selection, model-comparison, IBI-validation, BAS/CAS/IAS and
single-cell analysis scripts. These scripts preserve the original analysis
workflow and may require local data paths and additional R packages.
