"""Scientific invariants, independently computed examples, and pipeline contracts."""
from hashlib import sha256
import json
from pathlib import Path
from unittest.mock import patch

from lifelines import KaplanMeierFitter
import numpy as np
import pandas as pd
from PIL import Image
import pytest
from scipy.optimize import minimize

from survival_analysis.analysis import COVARIATES, FIGURES, TABLES, fit_cox, run, verify
from survival_analysis.data import SOURCE_SHA256, download, prepare, validate_raw, validate_source

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def raw():
    return pd.read_csv(ROOT / 'data/raw/lung.csv')


@pytest.fixture(scope='module')
def generated(tmp_path_factory):
    root = tmp_path_factory.mktemp('pipeline')
    (root / 'data/raw').mkdir(parents=True)
    (root / 'data/raw/lung.csv').write_bytes((ROOT / 'data/raw/lung.csv').read_bytes())
    run(root)
    return root


def test_pinned_source():
    payload = (ROOT / 'data/raw/lung.csv').read_bytes()
    assert sha256(payload).hexdigest() == SOURCE_SHA256
    assert len(validate_source(payload)) == 228
    provenance = json.loads((ROOT / 'data/raw/provenance.json').read_text())
    assert provenance['sha256'] == SOURCE_SHA256
    assert provenance['rows'] == 228


def test_corrupt_source_rejected():
    with pytest.raises(ValueError, match='SHA-256'):
        validate_source(b'bad download')


@pytest.mark.parametrize('column,value', [('status', 2), ('sex', 0), ('time', 0), ('age', -1), ('ph.ecog', 5), ('age', np.inf)])
def test_invalid_data_rejected(raw, column, value):
    raw[column] = raw[column].astype(float)
    raw.loc[0, column] = value
    with pytest.raises(ValueError):
        validate_raw(raw)


def test_missing_column_rejected(raw):
    with pytest.raises(ValueError, match='columns'):
        prepare(raw.drop(columns='time'))


def test_complete_case_counts_and_coding(raw):
    data = prepare(raw)
    assert len(data) == 227
    assert data.event.sum() == 164
    assert (1-data.event).sum() == 63
    assert data.female.sum() == 90
    assert 13 not in data.index  # The only missing required value is ECOG.
    assert data.loc[0, 'event'] == 1
    assert data.loc[2, 'event'] == 0
    assert data.loc[0, 'age_decade'] == 7.4
    # Missing unused measurements do not remove otherwise complete patients.
    assert data.loc[0, 'ecog'] == 1


def test_empty_complete_cases_rejected(raw):
    raw['ph.ecog'] = np.nan
    with pytest.raises(ValueError, match='No complete cases'):
        prepare(raw)


def test_km_hand_computed_censoring_and_ties():
    # Two deaths and a censor at time 2 share the risk set of four.
    # At time 3 the last remaining patient dies: S(2)=1-2/4=0.5, S(3)=0.
    km = KaplanMeierFitter().fit([2, 2, 2, 3], [1, 1, 0, 1])
    np.testing.assert_allclose(km.predict([0, 1, 2, 3]), [1, 1, 0.5, 0])


def test_km_all_censored():
    km = KaplanMeierFitter().fit([1, 2, 3], [0, 0, 0])
    np.testing.assert_allclose(km.survival_function_.values, 1)
    assert np.isinf(km.median_survival_time_)


def test_cox_matches_independent_efron_partial_likelihood(raw):
    data = prepare(raw)
    x = data[COVARIATES].to_numpy()
    x = x - x.mean(axis=0)
    time = data.time.to_numpy()
    event = data.event.to_numpy().astype(bool)
    death_times = np.unique(time[event])

    def objective(beta):
        eta = x @ beta
        weight = np.exp(eta)
        loglik = 0.0
        for t in death_times:
            deaths = (time == t) & event
            d = deaths.sum()
            risk_weight = weight[time >= t].sum()
            tied_weight = weight[deaths].sum()
            loglik += eta[deaths].sum()
            loglik -= np.log(risk_weight - np.arange(d) / d * tied_weight).sum()
        return -loglik

    fitted = minimize(objective, np.zeros(3), method='BFGS', options={'gtol': 1e-5})
    assert fitted.success, fitted.message
    np.testing.assert_allclose(fit_cox(data).params_.loc[COVARIATES], fitted.x, atol=2e-5, rtol=0)


def test_numerical_regression(generated):
    # Golden values recorded from the first pinned pipeline run, not external validation.
    summary = json.loads((generated / 'results/summary.json').read_text())
    assert summary['km_median_days'] == {'Male': 270.0, 'Female': 426.0}
    assert summary['logrank_p'] == pytest.approx(0.00155846981751, abs=1e-10)
    model = pd.read_csv(generated / 'results/cox_model.csv').set_index('covariate')
    np.testing.assert_allclose(model.loc[COVARIATES, 'exp(coef)'], [1.11702317483, 0.575445699784, 1.58999066108], rtol=1e-7)


def test_artifact_inventory_and_pngs(generated):
    assert verify(generated) == {'hashed_artifacts': 14, 'result_csv': 8, 'figures': 3}
    assert len(list((generated / 'results').glob('*.csv'))) == len(TABLES) == 8
    assert len(list((generated / 'results').iterdir())) == 14
    for name in FIGURES:
        with Image.open(generated / f'results/{name}.png') as image:
            image.verify()
        with Image.open(generated / f'results/{name}.png') as image:
            assert min(image.size) >= 600
    assert len(pd.read_csv(generated / 'results/schoenfeld.csv')) == 164


def test_at_risk_definition(generated, raw):
    data = prepare(raw)
    risk = pd.read_csv(generated / 'results/at_risk.csv')
    assert len(risk) == 12
    for row in risk.itertuples():
        group = data[data.female == int(row.group == 'Female')]
        assert row.at_risk_start == (group.time >= row.time).sum()


def test_artifact_tampering_detected(generated, tmp_path):
    import shutil
    root = tmp_path / 'tampered'
    shutil.copytree(generated, root)
    (root / 'results/ph_test.csv').write_text('tampered')
    with pytest.raises(ValueError, match='changed'):
        verify(root)


@pytest.mark.parametrize('failure', ['network', 'checksum'])
def test_failed_download_does_not_mutate_existing_files(tmp_path, failure):
    target = tmp_path / 'data/raw'
    target.mkdir(parents=True)
    (target / 'lung.csv').write_bytes(b'original data')
    (target / 'provenance.json').write_bytes(b'original provenance')
    before = {p.name: p.read_bytes() for p in target.iterdir()}
    if failure == 'network':
        with patch('survival_analysis.data.urlopen', side_effect=OSError('offline')):
            with pytest.raises(OSError):
                download(tmp_path)
    else:
        with patch('survival_analysis.data.urlopen') as request:
            request.return_value.__enter__.return_value.read.return_value = b'bad source'
            with pytest.raises(ValueError, match='SHA-256'):
                download(tmp_path)
    assert {p.name: p.read_bytes() for p in target.iterdir()} == before


def test_flagged_ph_is_disclosed(tmp_path, monkeypatch):
    from survival_analysis import analysis
    original = analysis.proportional_hazard_test

    def flagged(*args, **kwargs):
        result = original(*args, **kwargs)
        result._p_value[0] = 0.001
        return result

    (tmp_path / 'data/raw').mkdir(parents=True)
    (tmp_path / 'data/raw/lung.csv').write_bytes((ROOT / 'data/raw/lung.csv').read_bytes())
    monkeypatch.setattr(analysis, 'proportional_hazard_test', flagged)
    summary = run(tmp_path)
    assert len(summary['ph_flagged']) == 1
    assert 'Flagged covariates: ' + summary['ph_flagged'][0] in (tmp_path / 'results/report.md').read_text()
