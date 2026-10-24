"""Pinned public source and explicit validation of lifelines event coding."""
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import tempfile
from urllib.request import urlopen

import numpy as np
import pandas as pd

SOURCE_URL = "https://raw.githubusercontent.com/CamDavidsonPilon/lifelines/v0.30.0/lifelines/datasets/lung.csv"
SOURCE_SHA256 = "16ebd02c38bbb5708186d600cff80ad84d0764974619436f83ba4a605e77eb44"
REQUIRED = ["time", "status", "age", "sex", "ph.ecog"]
COVARIATES = ["age_decade", "female", "ecog"]


def validate_raw(raw: pd.DataFrame) -> None:
    if not set(REQUIRED).issubset(raw.columns):
        raise ValueError("Missing required columns")
    for column in REQUIRED:
        values = raw[column].dropna()
        if not pd.api.types.is_numeric_dtype(values) or not np.isfinite(values).all():
            raise ValueError(f"Nonfinite/nonnumeric values: {column}")
    if not set(raw.status.dropna()).issubset({0, 1}):
        raise ValueError("Expected lifelines status 0=censored, 1=death; not R 1/2")
    if not set(raw.sex.dropna()).issubset({1, 2}):
        raise ValueError("Expected sex 1=male, 2=female")
    if not set(raw['ph.ecog'].dropna()).issubset({0, 1, 2, 3, 4}):
        raise ValueError("Invalid ECOG score")
    if (raw.time.dropna() <= 0).any() or (raw.age.dropna() <= 0).any():
        raise ValueError("Time and age must be positive")


def validate_source(payload: bytes) -> pd.DataFrame:
    if sha256(payload).hexdigest() != SOURCE_SHA256:
        raise ValueError("Source SHA-256 mismatch; refusing replacement")
    raw = pd.read_csv(BytesIO(payload))
    validate_raw(raw)
    if len(raw) != 228:
        raise ValueError("Expected 228 source records")
    return raw


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as file:
        tmp = Path(file.name)
        try:
            file.write(payload)
            file.flush()
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)


def download(root: Path) -> None:
    # Fetch and validate before touching existing files, including on network failure.
    with urlopen(SOURCE_URL, timeout=60) as response:
        payload = response.read()
    raw = validate_source(payload)
    provenance = {
        "url": SOURCE_URL, "sha256": SOURCE_SHA256, "rows": len(raw),
        "dataset": "NCCTG advanced lung cancer (lifelines distribution)",
        "event_coding": {"0": "right censored", "1": "death"},
        "documentation": "https://stat.ethz.ch/R-manual/R-devel/library/survival/html/lung.html",
        "coding_note": "R survival uses 1/2; this lifelines CSV already uses 0/1.",
    }
    _atomic_write(root / "data/raw/lung.csv", payload)
    _atomic_write(root / "data/raw/provenance.json", (json.dumps(provenance, indent=2) + "\n").encode())


def prepare(raw: pd.DataFrame) -> pd.DataFrame:
    validate_raw(raw)
    complete = raw.dropna(subset=REQUIRED)
    if complete.empty:
        raise ValueError("No complete cases")
    return pd.DataFrame({
        "time": complete.time.astype(float), "event": complete.status.astype(int),
        "age_decade": complete.age / 10,
        "female": (complete.sex == 2).astype(int), "ecog": complete['ph.ecog'],
    }).rename_axis("source_row")
