import argparse
import json
from pathlib import Path

from .analysis import run, verify
from .data import download


def main():
    parser = argparse.ArgumentParser(description='NCCTG survival pipeline')
    parser.add_argument('command', choices=['download', 'analyze', 'verify'])
    parser.add_argument('--root', type=Path, default=Path.cwd(), help='Repository/output root')
    args = parser.parse_args()
    root = args.root.resolve()
    if args.command == 'download':
        download(root)
        print('Validated download: 228 rows, SHA-256 matched')
    else:
        print(json.dumps(run(root) if args.command == 'analyze' else verify(root), indent=2))


if __name__ == '__main__':
    main()
