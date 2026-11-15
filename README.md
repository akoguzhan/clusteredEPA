# clusteredEPA

Testing Clustered Equal Predictive Ability with Unknown Clusters

## Installation

```r
devtools::install_github("akoguzhan/clusteredEPA")
```

or, from the RStudio project of the paper, `source("Simulation/00_setup.R")` (installs
PanelKmeansInference first).

## Main functions (version 0.2.3)

* `epa_clustered_selective()`: selective C-EPA tests (F_SI, joint C-EPA, merged and joint
  homogeneity, weighted Bonferroni, O-EPA); with `cluster_specific = TRUE`, also the selective
  tests that each selected cluster mean is zero, with Holm-adjusted p-values.
* `epa_unit_tests()`, `power_enhance()`: unit-level benchmarks for joint EPA.
* `epa_clustered_known()`, `epa_clustered_split()`, `overall_EPA_test()`.
* `epa_default_B()`: the common B for all tests on a sample.
* `generate_epa_simulation_data()`: extended Monte Carlo design of the paper.
