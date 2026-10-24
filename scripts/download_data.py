"""Run from repository root after installing the project."""
from pathlib import Path
from survival_analysis.data import download

if __name__ == '__main__':
    download(Path.cwd())
    print('Downloaded and verified data/raw/lung.csv and provenance.json')
