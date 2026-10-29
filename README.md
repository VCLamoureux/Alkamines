# Alkamines

R analysis of fecal drug-derived alkamines and sulfasalazine exposure in a longitudinal rheumatoid arthritis cohort: **17 participants, including 2 users and 15 non-users, sampled at T1, T2 and T3**.

## Main script

[**julius_standalone_revise_script_updated.R**](scripts/julius_standalone_revise_script_updated.R) prepares the data, performs participant-level analyses, and generates longitudinal plots and boxplots. The original analysis and supporting MZmine files are also retained in `scripts/`.

## Run

Install the required R packages:

```r
install.packages(c("dplyr", "tidyr", "purrr", "readr", "ggplot2", "patchwork"))
```

Set `project_dir` and check the paths in `input_files` at the top of the main script. Four inputs are required:

- Feature peak-area table, retaining the `6Mix` and `Blank` QC columns.
- GNPS annotation table.
- Clinical metadata.
- Supplementary drug metadata.

**Only the feature table is included in this repository; the other three inputs must be supplied separately.** Expected filenames are defined in `input_files`.

Click **Source** in RStudio, or run from the repository root:

```r
source("scripts/julius_standalone_revise_script_updated.R")
```

The script creates `table_julius` directly from these inputs.

## Statistics

Raw peak areas are averaged across the three visits within each participant; only complete profiles enter the tests. Between-group mean differences are assessed using two-sided exact participant-label permutation tests based on absolute differences, enumerating all 136 allocations. P values are unadjusted for multiple comparisons. The tests compare average abundance across visits, not temporal change. Findings are exploratory given the two exposed participants.

Secondary detection analyses use participant-level Fisher's exact tests with an operational peak-area threshold of `>1` at any visit.

## Outputs

Results are saved in `<project_dir>/julius_standalone_results/`: the reconstructed analysis table, participant summaries, statistical results, QC audits, and longitudinal and boxplot figures in **PDF, SVG and PNG**. PNG resolution is controlled by `figure_dpi`; set it to `900` for 900 dpi exports.

## Citation

Agongo, J. et al. **Alkamines reveal a hidden layer of steroid and drug metabolism.** *bioRxiv* (2026), preprint. [doi:10.64898/2026.05.13.724743](https://doi.org/10.64898/2026.05.13.724743).
