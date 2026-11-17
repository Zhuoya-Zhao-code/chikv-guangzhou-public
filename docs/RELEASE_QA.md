# Release QA — scientific correction candidate v2 (2026-11-03)

This check concerns **code files only**. The original 12 uploaded R/Rmd scripts were used as source references but are **not** included in this deliverable.

## Changes requiring verification against the manuscript

- Primary descriptive source assignment now **includes background state 0**; `support_prob` refers to that selected state. The graph's edges, true isolates, cluster counts, weak components, generation depth, cluster diameter and OLS slope were rebuilt from this selection.
- All **21 network parameter sensitivity scenarios** use 10,000 networks each, and the **case-major** sampling order in the late Rmd. Primary reconstruction retains the original **simulation-major** order. Hence baseline sensitivity draws need not equal primary reconstruction draws; do not set their outputs equal by hand.
- Replaced Yates chi-square on pair-count aggregates with the source's **case-level SNP-cluster-label permutation** at ≤1/2/3/5 SNP, 9,999 draws, two-sided P and Holm adjustment.
- Extended 300 sampled-network outputs to include the final coefficient sign percentages, 2.5% and 97.5% quantiles and spline joint-Wald stability fields.
- Removed the incomplete upstream cleaning template; the prepared 1,280-case model input must come from the governed and independently verified original workflow. Added an automatic aggregate fingerprint comparison and expanded synthetic, privacy-preserving test cases. Figure 2A/B/C/D, 3A/B/C and S4 are analytic plotting *components*, not final manuscript layout replications.

## Gates already executed in the packaging environment

- `python3 scripts/static_audit.py` — PASS; checks file types, private paths, selected credential patterns, name-join patterns and delimiter balance. **Not an R syntax parser.**
- Git fresh-staging dry run with dummy private research files — PASS. Reviewed staged file paths; private sample files stayed ignored. `git diff --cached --check` — PASS.
- Reviewed the deliverable ZIP manifest for forbidden original datasets, raw Rmd notebooks, private files, FASTA, RDS and patient-level outputs — none included.

## Not executed, so cannot be claimed

- **R parse / unit tests:** no R executable is installed in this packaging container. All R script API calls, ggplot rendering and R package compatibility require a real R runtime.
- **Numerical concordance:** the authorized surveillance/sequence data, derived complete-case input and full parent matrix are not available here; sample counts, coefficients, cluster metrics, primary/sensitivity Monte Carlo draws, genomic P and graphs are **not verified**.
- Institutional authorization, intellectual-property rights, release approval, repository history review and independent scientific methods review are **not established**.

## Authorized-run commands (from repository root)

```bash
Rscript -e 'for(f in c(list.files("R", "\\.R$", full.names=TRUE),list.files("scripts", "\\.R$", full.names=TRUE),"tests/test_core.R")) parse(file=f); cat("R parser: PASS\n")'
Rscript tests/test_core.R
# Follow README for the primary pipeline; these need governed inputs.
Rscript scripts/10_check_manuscript_concordance.R
```

The concordance script checks the reported 12 key aggregate metrics and requires **300 successful sampled-network fits**. It does not substitute for manually checking AIC/coefficient values, uncertainty summaries, genomic permutation P and final figures. If it fails, stop before publication and trace the discrepancy to data, model or random-stream differences. Never tune code just to match published numbers.

**Release recommendation:** keep the repository PRIVATE until the authorized team passes R/runtime, result fingerprint, manual supplements, disclosure approval and Git history review. See `MANUSCRIPT_CONCORDANCE.md` and `PRIVACY_AND_RELEASE.md`.
