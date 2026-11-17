# Internal input contract (do not upload source datasets)

All paths are relative to the repository root. Input files below live under `data/private/` and **must not enter GitHub, attached issues, or logs**. These schemas are not an invitation to disclose real case-level values.

## `cases_network.csv` (governed individual-level, e.g. 1,329 analysis cases)

| Column | Format | Meaning |
|---|---|---|
| `case_id` | unique string | Consistent restricted pseudonymous case identifier; no names / original source IDs in Git |
| `onset_date` | ISO `YYYY-MM-DD` | Exact onset date, private |
| `case_type` | `Imported` or `Local` | Epidemiologic case classification |
| `utm_x` | numeric metres | Private residential UTM easting, EPSG:32649 |
| `utm_y` | numeric metres | Private residential UTM northing, EPSG:32649 |

**Coordinate system** matters: use metres, not raw longitude/latitude degrees or GCJ-02 longitude/latitude. Any conversion, geocoding, address joins and error checks must occur only under institutional authorization.

## Upstream data preparation (not included)

`Case_Data_Phase1_Ready.csv` must be produced by the institution's **full original governed preprocessing**. This project does not provide a standalone replacement for raw source linkage, residential geocoding, covariate-source merging, eligibility exclusions or restandardization. A previously provided partial upstream script was **removed** because it could silently yield a cohort or standardized variables different from those used by the paper.

## `Case_Data_Phase1_Ready.csv` (governed complete-case regression input, expected 1,280 rows)

The final-model script expects the already cleaned/standardized covariates below. It does not silently drop missing values.

| Column | Meaning / model definition |
|---|---|
| `case_id` | Exact mapping to restricted `cases_network.csv` (required for 300-refit analysis) |
| `mean_Ri` | From the **full 10,000 sampled networks** for the same case |
| `onset_date` | ISO date, converted to `onset_day = onset_date − minimum onset_date in analysis data` |
| `z_age` | Z-scored age |
| `sex` | `Female` or `Male` |
| `occ_recode` | `Indoor/Sedentary`, `Student/Childcare`, `High Exposure/Mobile`, `Home-based`, `Other/Unknown` |
| `case_type` | `Imported` or `Local` |
| `Reporting_Delay` | Reporting date minus onset date in days (non-negative) |
| `matched_P` | Environmental Risk Index matched to the case's exposure window |
| `z_BI` | Z-scored weekly Breteau Index |
| `z_outflow` | Z-scored population outflow |
| `z_dens_park` | Z-scored parks / green-space density or count as defined in internal analysis |
| `z_log_hosp` | Z-scored log healthcare-facility quantity |
| `z_log_trans` | Z-scored log distance to public transport |
| `log_dist_market` | **Unstandardized** log distance to the nearest agricultural market; verify `log(distance + 1)` in original preprocessing |
| `utm_x_km`, `utm_y_km` | Private residential UTM coordinates divided by 1,000 |

**Important:** The uploaded Rmd writes `Case_Data_Phase1_Ready.csv` after joins and normalization from personally identifiable surveillance sources. The public scripts do not rebuild those joins. Check factors and transformations (especially market log definition and regression case selection) with the original private pipeline and submitted Methods before claiming numerical agreement.

## Restricted inputs for the published genomic permutation analysis

- `genomic_cluster_labels_private.csv`: exactly one row per **clustered sequenced case**, in the same order as the pair-index file; columns `genomic_index` (1-based, ordered 1..N) and `cluster_id` (nonmissing transmission-compatible cluster assignment). No names, sequence identifiers or original case IDs are needed by this public algorithm.
- `genomic_pairwise_snp_private.csv`: columns `index_A`, `index_B`, `SNP_distance`, with integer indices into the previous file, strictly `index_A < index_B`, no duplicated pairs. Distances must be computed inside the authorized environment using the **same QC/sequence-filter choices** as the final manuscript, in particular the original Rmd's default `minimum_pairwise_comparable_fraction = NULL`. For the unfiltered analysis, there should be `choose(N,2)` unordered pairs; investigate any missing comparisons.

These inputs are **case-level / pair-level restricted data** even without names, and cannot be committed. The algorithm permutes case cluster labels without replacement with seed `20260807`, 9,999 draws, preserves label multiplicities, tests within-minus-between proportions at ≤1, ≤2, ≤3 and ≤5 SNP, and applies Holm adjustment to the four two-sided empirical P values. This matches the late Rmd's procedure and cannot be reconstructed from aggregate counts alone.

## Outputs

All outputs are written to `results/private/`. Some include case IDs, full precise coordinates, per-network case assignments or models whose objects may contain the original inputs; they are **not public-use derived data**. Neither renaming files nor removing a `name` column makes them safe to publish.
