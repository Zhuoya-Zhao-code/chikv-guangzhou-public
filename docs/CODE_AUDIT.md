# Scientific source audit — v2 candidate (2026-11-03)

This repository was derived from the **12 uploaded R and Rmd files**, not generated from private patient data. The earlier v1 ZIP contained deviations discovered during review. Version 2 corrects the **known code-level discrepancies** but has **not** passed R runtime or manuscript data-concordance tests.

## Source selection

| Original source | Status in v2 |
|---|---|
| `1 传播链重建.Rmd` | Canonical late corrective sections: source `0` included in per-case posterior modal parent, weak connected components, true isolates, generation depth; final sensitivity block uses 10,000 case-major samples/scenario, complete numeric summary and Wilson intervals; SNP block uses 9,999 label permutations |
| `可视化.Rmd` | Figure 2D cumulative onward-transmission method and associated manuscript plot semantics; layout drafts are not all retained |
| `2 多因素分析 3.Rmd` | Late covariate transformation descriptions only; identifiable raw joins omitted |
| `chikv_temporal_adjusted_reanalysis(1).R` | Primary nested Tobit M1–M3, onset spline, local-only and complete-window checks, prediction/diagnostics |
| `chikv_final_remaining_analysis(1).R` | Final 300 sampled-network direct-from-parent-matrix refitting, selection seed 20260708 |
| `updated_figureS4_tablesS4_S9(1).R` | Later labels and final summary figure context; display outputs remain to be checked |
| `code.R`, `code2.R`, `代码 1.R`, `2多因素分析*.Rmd`, `chikv_remaining_two_models(1).R` | Old or superseded exploratory implementations; not copied into public source tree |

## Confirmed v1-to-v2 fixes

1. **Descriptive source assignment:** v1 defaulted to the modal nonzero observed parent; late original block explicitly includes `0`. v2 uses `tabulate(parents + 1)` per case, chooses maximum frequency including `0`, and conservatively chooses `0` on ties. Source-support probability now corresponds to the chosen mode.
2. **Unconnected cases:** v1 reported `top_parent_id == 0` as unlinked. v2 defines isolated nodes as having **no incoming observed-source edge AND no outgoing observed-source edge**. A seed with descendants is not unconnected.
3. **Clusters:** weak components of a directed parent forest are counted as formal transmission-compatible clusters only when size ≥2. Singleton count, connected count, largest cluster, and **maximum generation depth** are all explicitly computed. Duration–diameter OLS includes all size ≥2 clusters, including zero-duration ones (as in the late source).
4. **Sensitivity:** v1 used 1,000 samples/scenario and a simplified output. v2 uses **21 scenarios × 10,000 samples**, preserving final source parameter grid, the sensitivity-specific case-major RNG order, Wilson intervals and the original numeric metrics for the sensitivity supplementary figure.
5. **Genomic inference:** v1 applied a Yates-corrected chi-square test on aggregate counts. v2 eliminates it and provides the late source's **case-level cluster-label permutation** with 9,999 draws, two-sided empirical P and Holm adjustment. Aggregate proportions do not contain enough information to reproduce the original permutation test.
6. **300-network summary fields:** v2 preserves linear-coefficient direction percentages and quantiles and spline joint-test summaries, needed by late supplementary figure and tables. Patient/network mapping remains strictly ID-based.
7. **Fail-closed concordance:** added a private-report script against 12 core manuscript-level targets and all 300 refit statuses; no value is hard-coded into the scientific estimator to force a match.

## Remaining uncertainties / intentional limitations

- Primary draw order **simulation-major** and late sensitivity draw order **case-major** differ in the original sources; v2 preserves both. Identical seeds cannot make the resulting matrices bitwise equal. Exact matching still requires the same input ordering, package versions and R RNG behavior.
- The source and paper often say “maximum-posterior network”. The algorithm is **per-case marginal modal source selection**, not necessarily a *joint* maximum-posterior sampled network. This terminology should be discussed with the authors rather than silently changing the method.
- `Case_Data_Phase1_Ready.csv` is derived from protected multi-source records. Full geocoding, personal linkage, POI harvesting, exclusions and complete standardization are not included. The previous *partial* upstream ETL script was removed to prevent accidentally regenerating a different analysis-ready cohort.
- Figures are regenerable computational components, not certified pixel-identical to the submitted artwork. Some multi-panel layouts remain outside the code release. Published figure/data concordance must be checked with approved inputs.
- No actual R parse/test or model fitting was possible in the packaging container (R runtime unavailable); no underlying patient/sequence data were provided. Static passing does **not** prove numerical reproducibility.
- Rights-holder approval, code license selection and disclosure clearance are not established. Individual-level case, coordinate, pair-SNP and genome linkage records must remain private.

**Status: static candidate only, NOT approved for public GitHub release.**
