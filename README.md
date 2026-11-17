# Chikungunya transmission-compatible network reconstruction — Guangzhou, 2025

**Research-code release candidate v2 (2026-11-03) — restricted patient data are NOT included.**

Code supporting the computational analyses in the manuscript *“Transmission-compatible network reconstruction characterizes the organization of chikungunya transmission during the 2025 Guangzhou outbreak”* (working title). This repository consolidates the final analysis specifications from multiple R/R Markdown drafts; old and exploratory versions are intentionally excluded. **It has not yet been verified against the restricted original dataset, and an R runtime was not available for end-to-end execution during packaging.** Before labeling a release reproducible, the authorized research team must run the tests and compare all outputs with the accepted manuscript.

## What it does

| Script | Scope | Main private input |
|---|---|---|
| `scripts/00_knox.R` | Space–time Knox permutation tests (all and local cases) | `cases_network.csv` |
| `scripts/01_reconstruct.R` | 10,000 sampled transmission-compatible parent networks and case-level onward-transmission summaries | `cases_network.csv` |
| `scripts/02_cluster_summary.R` | Reconstructed network clusters, residential extents and durations | network output |
| `scripts/03_parameter_sensitivity.R` | Alternative spatial scales, cutoffs, background weights and serial intervals | `cases_network.csv` |
| `scripts/04_tobit_models.R` | Three temporally adjusted nested Tobit models; right-boundary and local-only checks; spatial residual diagnostics; marginal latent predictions | `Case_Data_Phase1_Ready.csv` |
| `scripts/05_sampled_network_uncertainty.R` | Refit Model 3 over 300 draws, matched by case ID to **the complete 10,000-network matrix** | prepared cases + private parent matrix |
| `scripts/06_figures.R` | Figure 2D, Figure 3 component curves, Figure S4 components | private derived results |
| `scripts/07_genomic_consistency_permutation.R` | **Manuscript case-level 9,999-permutation SNP test**, with Holm adjustment | governed clustered-sequence labels and pairwise SNP-distance inputs |
| `scripts/08_optional_tweedie.R` | Optional Tweedie Model 3 sensitivity analysis | prepared cases |
| `scripts/09_sensitivity_figure_S5.R` | Figure S5 aggregated sensitivity profile | finalized 21-scenario sensitivity output |
| `scripts/10_check_manuscript_concordance.R` | Fail-closed comparison against published aggregate fingerprint and 300-fit convergence | approved private output files |

The scripts regenerate computational results, **not a final typeset figure set**. Diagram positioning, manuscript figure numbering, Table 1, external POI acquisition, vector/ecology preprocessing and raw sequence alignment are out of scope of this release candidate. Their use depends on internal data-governance procedures.

## Final analysis specifications

- Gamma-distributed **onset-to-onset serial interval**: mean **8 days**, standard deviation **3 days**.
- Candidate parent occurs **strictly earlier** (1–30 days) and within **1,300 m** of recipient's residential coordinates.
- Spatial compatibility weight `exp(-d / 250)` with `d` in metres; background weight `1e-4`.
- Imported cases are assigned background state `0`; local cases may also have background source.
- **10,000** sampled parent networks, reconstruction random seed `2026`; 300 selected network columns for regression robustness, selection seed `20260708`.
- `mean_Ri` is the model-inferred expected **direct** offspring count of a case over sampled networks, **not** a directly observed reproduction number `R_t`.
- Regression Model 1: demographics, case classification, reporting-delay natural spline (3 df), symptom-onset natural spline (4 df), spatial second-degree polynomial. Model 2 adds ERI natural spline (3 df), Breteau Index and population outflow. Model 3 additionally adds parks, healthcare and public-transport variables, plus market-distance natural spline (3 df). All Tobit models are left-censored at zero.

**Important reconstruction caveat:** The **late corrective section** in the original reconstruction Rmd explicitly replaced the earlier selection of the most frequent *non-background* source. Version 2 selects each case's most frequent state **including background 0**; ties prefer 0 and otherwise the smallest numeric source ID. Clusters, links and unconnected cases now derive from this corrected state selection. An "unconnected case" has **neither an observed-source upstream nor an observed-source downstream**. Clusters are weakly connected components with **at least two cases**. This is a per-case *marginal posterior mode* assignment and **not** a jointly optimized maximum-posterior network; the manuscript's descriptive terminology should be reviewed accordingly.

**Simulation reproducibility:** Primary reconstruction uses the original simulation-major RNG draw order, with 10,000 networks and seed `2026`. The manuscript's late sensitivity section uses **case-major** draws, also 10,000 **per scenario** (21 scenarios total) and seed `2026`. These two orders consume the random stream differently: do not mix them or assume the sensitivity baseline will be bitwise identical to the primary baseline. The 300-refit sample is selected from the primary saved matrix using seed `20260708`.

**Genomics:** The paper's final approach compares proportions with SNP distances `<=1, <=2, <=3, <=5` within vs between reconstructed clusters using **9,999 case-level cluster-label permutations**, two-sided empirical P values, and Holm adjustment (seed `20260807`). The previous aggregate-count chi-square implementation was removed because it does **not** reproduce this test. The script runs only on local, governed pair-index/cluster-label inputs; sequence/name matching and raw data are not publicly included.

## Running locally on approved restricted data

Use an R environment with `readr`, `AER`, `spdep`, `car`, `ggplot2`, `scales` and `mgcv` (the plotting/alternative-model packages are optional). `stats` and `splines` are included with R. Install, where allowed:

```r
install.packages(c("readr", "AER", "spdep", "car", "ggplot2", "scales", "mgcv"))
```

Place permitted data **only on a secure machine** under `data/private/` (git-ignored). After extracting the ZIP, from the repository root create these empty local folders and run:

```bash
mkdir -p data/private results/private
Rscript tests/test_core.R
Rscript scripts/00_knox.R
Rscript scripts/01_reconstruct.R
Rscript scripts/02_cluster_summary.R
Rscript scripts/03_parameter_sensitivity.R
Rscript scripts/04_tobit_models.R
Rscript scripts/05_sampled_network_uncertainty.R
Rscript scripts/06_figures.R
Rscript scripts/09_sensitivity_figure_S5.R
# Only on the controlled computer, after approved per-case SNP/cluster inputs exist:
Rscript scripts/07_genomic_consistency_permutation.R
# Only after ALL analyses and their reviewed manuscript targets exist:
Rscript scripts/10_check_manuscript_concordance.R
```

The 21-scenario sensitivity run can be time-consuming; it is necessary for a complete scientific-concordance check, but may be skipped during quick smoke tests. `scripts/08_optional_tweedie.R` is optional. All generated files, including tables and graphics, go under `results/private/` and are ignored by Git. **Do not force-add them.** Input schema: `docs/DATA_DICTIONARY.md`.

### Input caution

`cases_network.csv` must be ordered consistently with the original reconstruction before tie sorting; code uses ascending `onset_date` with stable original row order. `case_id` must match across analysis and network files (no row-position fallback). `Case_Data_Phase1_Ready.csv` is the **already verified, restricted, final-model input**. Reproduce and quality-check it internally using the full original governed preprocessing workflow. The earlier proposed partial preprocessing template has been **removed** from this GitHub release because it was not sufficient to guarantee reproduction of the 1,280-case sample and standardized variables. The public repository excludes personally identifying linkage, source merging, address geocoding and POI API collection.

## Privacy, code review and validation

**The scripts are the only intended publication material.** No individual-level records, case IDs, addresses, residential coordinates, exact case dates, genomes, linkage maps, API credentials, raw monitoring data or case-level analysis outputs may be uploaded. Even deidentified coordinates, combined dates and genomic pair distances may permit reidentification. The `.gitignore` is defense in depth, **not** a substitute for authorization or examining Git history. See `docs/PRIVACY_AND_RELEASE.md`, `docs/CODE_AUDIT.md` and `docs/RELEASE_QA.md` before GitHub publication.

This candidate contains no claim of completed end-to-end results reproduction or public access to patient-level inputs. Data access remains governed by the originating public health authority. Do not ask GitHub users to contact the authors for unrestricted access or claim that the code alone reproduces tables without the restricted inputs.

## Maintainer release check

1. On a governed research computer, run `Rscript tests/test_core.R`, the parse command in `docs/RELEASE_QA.md`, and all applicable analysis scripts.
2. Run `scripts/10_check_manuscript_concordance.R`, then manually compare all nested model fits, tests, genomic Holm P values, 300-fit convergence, figure source data and supplementary tables to the **final submitted manuscript**. See `docs/MANUSCRIPT_CONCORDANCE.md`.
3. Review all output and inputs, screenshots, metadata, commit history and repository collaborators for reidentification or access concerns. Obtain data-owner permission for code and for aggregate results to be disclosed.
4. Add the agreed code license only after the research team/institution approves rights. **No open-source license is granted automatically by this draft repository.**

**Development status:** source-aligned v2 candidate. Static checks have been run; R execution on approved restricted data and manuscript numerical concordance are still pending. **Do not switch to a public repository yet.**
