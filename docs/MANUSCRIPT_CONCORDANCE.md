# Manuscript-to-code concordance gate (requires approved data)

The command `Rscript scripts/10_check_manuscript_concordance.R` reads **private aggregates produced by the scripts**, compares them with the reported main text numbers, and **fails if any differs beyond printed rounding**. It does not expose any case-level rows. It also requires all 300 sampled-network refits to converge successfully.

## Hard checks (reported figures, NOT implementation constants)

| Metric | Manuscript | Tolerance |
|---|---:|---:|
| Cases in reconstruction | 1,329 | exact |
| Cases in regression | 1,280 | exact |
| Observed-source links | 720 | exact |
| Transmission-compatible clusters | 206 | exact |
| Cases in multi-case clusters | 926 | exact |
| Truly unconnected cases | 403 | exact |
| Largest cluster size | 42 | exact |
| Maximum generation depth | 12 | exact |
| Median residential link distance | 567.3 m | 0.11 m |
| Links within 500 m | 46.1% | 0.11 percentage points |
| Cluster diameter–duration OLS slope | 32.3 m/day | 0.11 m/day |
| Highest 20% of cases' share of inferred onward transmission | 61.5% | 0.11 percentage points |

**The script does NOT replace the paper's data.** If it fails, check the original person ordering, parent matrix and parameters, the late background-inclusive source selection, versions/seeds and sample restrictions. Do not overwrite the reported statistics or tune the code to produce the fingerprint without methodological review.

## Manually check in addition

- Original data QC and expected exclusions, case-type factor coding, matched complete-case 1,280 rows, standardized variables and model formula terms.
- Model 1/2/3 AIC, coefficient estimates, joint Wald tests and likelihood-ratio tests; all nested regressions must use the same cohort.
- Main Figure 2A/B/C/D and Figure 3A/B/C source data, interval definitions, labels and axes.
- The 300-network stability tables and Figure S4, including selection seed, coefficient sign, 2.5/97.5 quantiles, and exact significance proportions.
- SNP pair generation and sequence QC; the two-sided permutation P values and Holm-adjusted results for 1/2/3/5 SNP. The 9999 permutations are **at case-label level**, not pair level.
- All supplementary robustness analyses and final submitted figure numbering; graphic formatting and tables in this repository remain *components*, not a certified facsimile of the manuscript.

## Privacy boundary

Run only on a governed workstation. Keep `data/private/`, `results/private/`, per-case labels, pairwise SNP files, exact dates, residential coordinates, logs and model objects out of Git. `.gitignore` is necessary but not sufficient: check the actual staged changes and Git history, and obtain institutional approval.
