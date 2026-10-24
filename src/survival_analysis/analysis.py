"""Deterministic complete-case analysis and auditable result artifacts."""
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, proportional_hazard_test

from .data import COVARIATES, REQUIRED, prepare, validate_source

TABLES = ["cohort_summary", "missingness", "km_curves", "at_risk", "logrank", "cox_model", "ph_test", "schoenfeld"]
FIGURES = ["kaplan_meier", "hazard_ratios", "schoenfeld"]
LABELS = {"age_decade": "Age (per 10 years)", "female": "Female vs male", "ecog": "ECOG (per point)"}


def fit_cox(data):
    return CoxPHFitter(penalizer=0.0).fit(data, duration_col="time", event_col="event")


def run(root: Path) -> dict:
    raw = validate_source((root / "data/raw/lung.csv").read_bytes())
    data = prepare(raw)
    out = root / "results"
    out.mkdir(parents=True, exist_ok=True)
    (root / "data/processed").mkdir(parents=True, exist_ok=True)
    data.to_csv(root / "data/processed/analysis.csv", index=True)
    tables = {}
    tables['missingness'] = pd.DataFrame({"column": REQUIRED, "missing": [int(raw[c].isna().sum()) for c in REQUIRED]})
    groups = [("Overall", data), ("Male", data[data.female == 0]), ("Female", data[data.female == 1])]
    tables['cohort_summary'] = pd.DataFrame([
        {"group": name, "n": len(d), "events": int(d.event.sum()),
         "censored": int((1-d.event).sum()), "median_age_years": float(d.age_decade.median()*10),
         "median_observed_time_days": float(d.time.median())} for name, d in groups
    ])
    fig, ax = plt.subplots(figsize=(9, 7))
    fitters, curves, risk = [], [], []
    medians = {}
    ticks = [0, 200, 400, 600, 800, 1000]
    for name, d in groups[1:]:
        km = KaplanMeierFitter(alpha=0.05).fit(d.time, d.event, label=name)
        fitters.append(km)
        km.plot_survival_function(ax=ax, ci_show=True, show_censors=True)
        curves.append(pd.DataFrame({"group": name, "time": km.survival_function_.index,
                                   "survival": km.survival_function_.iloc[:, 0].values,
                                   "ci_lower": km.confidence_interval_.iloc[:, 0].values,
                                   "ci_upper": km.confidence_interval_.iloc[:, 1].values}))
        medians[name] = float(km.median_survival_time_)
        for t in ticks:
            risk.append({"group": name, "time": t, "at_risk_start": int((d.time >= t).sum())})
    ax.set(title="NCCTG lung cancer: complete-case survival by recorded sex",
           xlabel="Time (days)", ylabel="Survival probability", ylim=(0, 1.02))
    ax.set_xticks(ticks)
    add_at_risk_counts(*fitters, ax=ax, rows_to_show=['At risk'], xticks=ticks,
                       at_risk_count_from_start_of_period=True)
    fig.subplots_adjust(bottom=0.30)
    fig.savefig(out / "kaplan_meier.png", dpi=150)
    plt.close(fig)
    tables['km_curves'] = pd.concat(curves, ignore_index=True)
    tables['at_risk'] = pd.DataFrame(risk)
    male, female = groups[1][1], groups[2][1]
    lr = logrank_test(male.time, female.time, male.event, female.event)
    tables['logrank'] = pd.DataFrame([{"comparison": "Male vs female", "test_statistic": lr.test_statistic,
                                      "degrees_of_freedom": 1, "p": lr.p_value}])
    cph = fit_cox(data)
    tables['cox_model'] = cph.summary.rename_axis('covariate').reset_index()
    ph = proportional_hazard_test(cph, data, time_transform="rank")
    tables['ph_test'] = ph.summary.rename_axis('covariate').reset_index()
    tables['ph_test']['flag_p_lt_0_05'] = tables['ph_test'].p < 0.05
    residuals = cph.compute_residuals(data, kind="scaled_schoenfeld")
    residuals['time'] = data.loc[residuals.index, 'time']
    tables['schoenfeld'] = residuals.rename_axis('source_row').reset_index()

    fig, ax = plt.subplots(figsize=(7, 4))
    s = cph.summary.loc[COVARIATES]
    hr = s['exp(coef)'].to_numpy()
    ax.errorbar(hr, range(3), xerr=np.vstack([hr-s['exp(coef) lower 95%'], s['exp(coef) upper 95%']-hr]), fmt='o', capsize=4)
    ax.set_yticks(range(3), [LABELS[c] for c in COVARIATES])
    ax.axvline(1, color='grey', linestyle='--')
    ax.set(xscale='log', xlabel='Adjusted hazard ratio (95% Wald CI)', title='Prespecified Cox model')
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out / 'hazard_ratios.png', dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    ordered = residuals.sort_values('time', kind='stable')
    for ax, cov in zip(axes, COVARIATES):
        ax.scatter(ordered.time, ordered[cov], alpha=0.4, s=15)
        # Descriptive rolling mean, not a fitted time-varying coefficient or CI.
        smooth = ordered[cov].rolling(25, center=True, min_periods=12).mean()
        ax.plot(ordered.time, smooth, color='darkred', label='25-event centered rolling mean')
        ax.axhline(0, color='grey', linestyle='--')
        ax.set(ylabel='Scaled residual', title=f"{LABELS[cov]} | rank-test p={ph.summary.loc[cov, 'p']:.4g}")
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel('Event time (days); one point per observed death')
    fig.tight_layout()
    fig.savefig(out / 'schoenfeld.png', dpi=150)
    plt.close(fig)
    for name, table in tables.items():
        table.to_csv(out / f'{name}.csv', index=False, float_format='%.12g')
    summary = {
        'source_n': len(raw), 'analysis_n': len(data), 'excluded_n': len(raw)-len(data),
        'events': int(data.event.sum()), 'censored': int((1-data.event).sum()),
        'km_median_days': medians, 'logrank_statistic': float(lr.test_statistic),
        'logrank_p': float(lr.p_value), 'concordance_in_sample': float(cph.concordance_index_),
        'ph_alpha': 0.05, 'ph_transform': 'rank',
        'ph_flagged': tables['ph_test'].loc[tables['ph_test'].flag_p_lt_0_05, 'covariate'].tolist(),
        'covariates': COVARIATES,
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    lines = ["# Computed analysis report", "", "Generated by `python -m survival_analysis analyze`; all estimates use the same complete-case cohort.", "",
             f"Source: {len(raw)} records; analyzed: {len(data)}; excluded: {len(raw)-len(data)}; deaths: {summary['events']}; censored: {summary['censored']}.", "",
             "## Kaplan–Meier and log-rank", "",
             f"Median survival: male {medians['Male']:.0f} days; female {medians['Female']:.0f} days.",
             f"Unadjusted log-rank chi-square(1)={lr.test_statistic:.6f}, p={lr.p_value:.8g}.", "",
             "## Adjusted Cox model", "", "| Covariate | HR | 95% CI | Wald p |", "|---|---:|---|---:|"]
    for cov, row in s.iterrows():
        lines.append(f"| {LABELS[cov]} | {row['exp(coef)']:.4f} | {row['exp(coef) lower 95%']:.4f}–{row['exp(coef) upper 95%']:.4f} | {row['p']:.6g} |")
    lines += ["", f"Apparent (training-data) concordance: {cph.concordance_index_:.6f}; not external predictive validation.", "",
              "## Proportional hazards diagnostics", "", "Scaled Schoenfeld residual rank-time tests, unadjusted alpha=0.05:", "",
              "| Covariate | p | Flag |", "|---|---:|---|"]
    for _, row in tables['ph_test'].iterrows():
        lines.append(f"| {row.covariate} | {row.p:.6g} | {bool(row.flag_p_lt_0_05)} |")
    lines += ["", "Flagged covariates: " + (", ".join(summary['ph_flagged']) or "none") + ".",
              "A non-significant test does not establish proportional hazards. The female residual smoother rises during sparse late follow-up; a constant hazard ratio remains uncertain. Inspect the residual plots; rolling means are descriptive, especially unstable in the tail. Tests are covariate-wise, not a global PH test; no multiplicity adjustment.", "",
              "## Interpretation and limitations", "",
              "Associations in a historical observational advanced-cancer cohort do not establish causal effects or treatment benefit. Independent censoring and linear age/ECOG effects are assumed. Complete-case analysis can be biased. Sparse late follow-up limits precision.", "",
              "External validation: not done. Imputation: not done. Nonlinear functional-form assessment: not done. Time-varying coefficients or PH-remedial refits: not done. Global PH test: not done. Causal inference: not done."]
    (out / 'report.md').write_text('\n'.join(lines) + '\n')
    paths = [out / f'{name}.csv' for name in TABLES] + [out / f'{name}.png' for name in FIGURES]
    paths += [out / 'summary.json', out / 'report.md', root / 'data/processed/analysis.csv']
    manifest = {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest() for path in paths}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    verify(root)
    return summary


def verify(root: Path) -> dict:
    out = root / 'results'
    manifest = json.loads((out / 'manifest.json').read_text())
    expected = {f'results/{name}.csv' for name in TABLES} | {f'results/{name}.png' for name in FIGURES}
    expected |= {'results/summary.json', 'results/report.md', 'data/processed/analysis.csv'}
    if set(manifest) != expected:
        raise ValueError('Unexpected artifact inventory')
    for name, digest in manifest.items():
        path = root / name
        if not path.is_file() or path.stat().st_size == 0 or sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Artifact missing, empty or changed: {name}')
    cox = pd.read_csv(out / 'cox_model.csv')
    ph = pd.read_csv(out / 'ph_test.csv')
    if set(cox.covariate) != set(COVARIATES) or len(cox) != 3 or set(ph.covariate) != set(COVARIATES) or len(ph) != 3:
        raise ValueError('Expected exactly three covariates')
    estimates = cox[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
    if not np.isfinite(estimates.to_numpy()).all() or not (estimates.iloc[:, :3] > 0).all().all():
        raise ValueError('Invalid Cox estimates')
    if not ((estimates.iloc[:, 1] <= estimates.iloc[:, 0]) & (estimates.iloc[:, 0] <= estimates.iloc[:, 2])).all():
        raise ValueError('Invalid confidence interval')
    for table in [cox, ph, pd.read_csv(out / 'logrank.csv')]:
        if not np.isfinite(table.p).all() or not table.p.between(0, 1).all():
            raise ValueError('Invalid p-value')
    if not np.isfinite(ph.test_statistic).all():
        raise ValueError('Invalid PH statistic')
    km = pd.read_csv(out / 'km_curves.csv')
    for _, curve in km.groupby('group'):
        if not curve[['survival', 'ci_lower', 'ci_upper']].apply(lambda x: x.between(0, 1).all()).all():
            raise ValueError('Invalid KM bounds')
        if (curve.survival.diff().dropna() > 1e-10).any():
            raise ValueError('Increasing KM survival')
    return {'hashed_artifacts': len(manifest), 'result_csv': len(TABLES), 'figures': len(FIGURES)}
