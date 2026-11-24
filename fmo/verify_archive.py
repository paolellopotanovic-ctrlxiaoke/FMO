"""Verify archived files and optionally reconstruct split evidence payloads."""
from pathlib import Path
import argparse
import hashlib
import json
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument('--reconstruct', action='store_true')
args = parser.parse_args()
root = Path(__file__).resolve().parent
manifest = json.loads((root / 'SOURCE_MANIFEST.json').read_text())

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while data := f.read(1024 * 1024):
            h.update(data)
    return h.hexdigest()

direct = [entry for entry in manifest['files'] if 'payload' not in entry]
for entry in direct:
    path = root / entry['file']
    assert path.stat().st_size == entry['bytes'], path
    assert sha(path) == entry['sha256'], path
for entry in manifest['payloads']:
    digest = hashlib.sha256()
    for part in entry['parts']:
        path = root / part['file']
        assert path.stat().st_size == part['bytes'], path
        assert sha(path) == part['sha256'], path
        with path.open('rb') as f:
            while data := f.read(1024 * 1024):
                digest.update(data)
    assert digest.hexdigest() == entry['sha256'], entry['file']
    if args.reconstruct:
        output = root / entry['file']
        with output.open('wb') as dest:
            for part in entry['parts']:
                with (root / part['file']).open('rb') as source:
                    while data := source.read(1024 * 1024):
                        dest.write(data)
        assert sha(output) == entry['sha256'], output
        with tarfile.open(output) as tar:
            for member in tar:
                if member.isfile():
                    relative = 'results/openfe/' + member.name
                    source = next((x for x in manifest['files'] if x['file'] == relative), None)
                    assert source is not None, relative
                    assert hashlib.sha256(tar.extractfile(member).read()).hexdigest() == source['sha256'], relative
        print('Reconstructed and verified:', output.name)
print('Verified source entries:', len(manifest['files']), '; direct files:', len(direct))
