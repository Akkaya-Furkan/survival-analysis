# P30 scope and acceptance criteria

Date: 2026-10-24. One-day methodological project; no clinical decision support.

## Question and fixed analysis

In the public NCCTG advanced lung cancer cohort, how does observed survival differ
by recorded sex, and how are age, sex and physician ECOG score associated with
mortality after mutual adjustment?

Use the lifelines v0.30.0 lung CSV (228 source records), downloaded from a versioned
URL with an expected SHA-256 digest. This copy encodes death as 1 and censoring as
0; do not apply the R survival dataset's 1/2 recoding. Use one complete-case cohort
for all estimates, dropping missing values only in time, status, age, sex, ph.ecog.
Age is modeled per 10 years, female versus male, ECOG as an ordinal linear term.
No selection of covariates based on significance. Cox uses Efron ties, Breslow
baseline hazard and no penalization. Report Wald 95% confidence intervals.

Produce sex-specific Kaplan–Meier estimates with pointwise 95% intervals and
at-risk counts; one unweighted two-sided log-rank test; three adjusted Cox hazard
ratios; rank-transformed scaled Schoenfeld tests for all three covariates at
unadjusted alpha 0.05; residual-versus-time plots with descriptive smoothing.
Flag violations rather than silently modifying the prespecified model.

## Numerical acceptance criteria

1. Download and validate exactly 228 rows, binary event coding and source SHA-256;
   save provenance. Report exact retained/excluded/event/censor counts.
2. Save a processed complete-case CSV, exactly 8 result CSV tables, 3 nonempty PNG
   figures (KM, hazard-ratio forest, Schoenfeld diagnostics), summary JSON,
   generated Markdown report and a SHA-256 artifact manifest.
3. Cox and PH tables each contain exactly 3 covariates; all estimates/statistics
   are finite, HRs positive, and each confidence interval contains its estimate.
   KM curves are bounded by [0,1] and nonincreasing. All reported p-values lie
   in [0,1]. Every flagged PH test is disclosed in the report and README.
4. At least 12 pytest cases pass, including hand-computable censoring/tied-event
   KM behavior, source/coding rejection, missing-data handling, regression
   reference values, and end-to-end artifact verification.
5. A pinned Python 3.11 environment reproduces the analysis offline from the
   checked-in input; a separate script fetches/validates the source. CI downloads
   the source, executes tests, regenerates outputs, verifies the manifest and
   uploads artifacts. Read the actual remote CI conclusion after push.
6. English README states question, approach, computed results, sources,
   reproduction commands and limitations. Unimplemented scientific extensions
   are explicitly marked "not done". Commit and push main, then verify the remote
   commit hash against the local hash.

## Excluded from this day

TCGA integration, external validation, causal inference, imputation, nonlinear
functional-form assessment, time-varying coefficients, treatment recommendations,
and deployment. No thesis-ligand or homelab work. No vault writes.
