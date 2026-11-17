# Privacy and release controls

## What must never be copied to this public repository

- Any original surveillance spreadsheet, case-line list, symptom-onset/report dates linked to a person, addresses, precise residential or activity coordinates, anonymized record-level tables, crosswalks, case identifiers or sample identifiers.
- Whole-genome FASTA/sequence files, raw paired SNP distances, name-to-sequence mappings, linkage and household/contact/medical information, or cohort membership identifying a person.
- Network matrices (`sim_parents_10000.rds`), reconstructed edges, per-case `Ri`, stored R model objects and geospatial maps derived from restricted inputs, absent specific permission.
- Authentication or API keys (especially former POI service credentials), personal machine paths, `.Rhistory`, `.RData`, GitHub Actions logs or run logs quoting restricted records.
- Dataset screenshots or figures with individual-case points. Suppressing column headers does not eliminate disclosure risk.

## Review requirements

- Use a new Git repository initialized **from the cleaned directory**. Never fork or push the original folder; ignored files can persist in historical commits, local branches and tags.
- Obtain explicit approval of **code disclosure and aggregate output disclosure** from the data owner/institution; patient data access remains restricted independently of code licensing.
- Apply `.gitignore` *before the first commit*; inspect `git status --short` and `git diff --cached --name-only` plus `git diff --cached --check` and inspect actual staged contents.
- Scan all staged content and commit history for keys, tokens, user home paths, identity fields, addresses and private file formats. Run `python3 scripts/static_audit.py`.
- Confirm collaborative ownership and determine the license before making the repository public.
- If a credential was ever committed, removing it from the working tree is insufficient: revoke/rotate it, clean history and evaluate exposure with the issuer.

## Public claims

"Code is publicly available; individual-level surveillance, residential location and genomic linkage data are not publicly shared due to privacy and data-governance restrictions. Reproduction of case-level results requires permission from the data custodian." Do not imply code publication guarantees unrestricted access to restricted patient data.
