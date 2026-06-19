#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f flake.nix ]]; then
  echo "error: run this from the dynsyss repository root" >&2
  exit 1
fi

cp -p flake.nix flake.nix.before-ds-repair

python3 <<'PY'
from pathlib import Path
import re

path = Path("flake.nix")
text = path.read_text()

pattern = re.compile(r'(?ms)^(?P<indent>\s*)tpcas\s*=\s*\{.*?^(?P=indent)\};')
match = pattern.search(text)
if not match:
    raise SystemExit("error: could not find the tpcas input block in flake.nix")

block = match.group(0)
block = re.sub(
    r'(?m)^(\s*)url\s*=\s*"[^"]*";',
    r'\1url = "github:hydrastro/tpcas";',
    block,
    count=1,
)
text = text[:match.start()] + block + text[match.end():]

if not re.search(r'(?m)^\s*ds-src\s*=\s*\{', text):
    match = pattern.search(text)
    assert match is not None
    ds_block = '''

    # DS is a separate source dependency of TPCAS.
    ds-src = {
      url = "github:hydrastro/ds/9e3224b3aef9b6e271d02b76ff671b1db4301601";
      flake = false;
    };'''
    text = text[:match.end()] + ds_block + text[match.end():]

path.write_text(text)
PY

rm -f flake.nix.rej Makefile.rej

echo "Repaired flake.nix. Review with:"
echo "  git diff -- flake.nix"
echo
echo "Then refresh the two source inputs with:"
echo "  nix flake lock --update-input tpcas --update-input ds-src"
