#!/usr/bin/env bash
set -euo pipefail

cd "${1:-.}"

if [[ ! -f Makefile ]]; then
  echo "error: Makefile not found in $(pwd)" >&2
  exit 1
fi

cp -p Makefile Makefile.before-ds-root-repair

python3 - <<'PY'
from pathlib import Path
import re

path = Path("Makefile")
text = path.read_text()
original = text

# 1. Remove stale Makefile assignments that force DS back under TPCAS.
#    An ordinary makefile assignment overrides the exported shell variable.
lines = text.splitlines()
out = []
saw_ds_root = False

assignment_re = re.compile(
    r'^(?P<indent>\s*)(?:override\s+)?DS_ROOT\s*(?P<op>[:?+]?=)\s*(?P<rhs>.*)$'
)

for line in lines:
    m = assignment_re.match(line)
    if m:
        saw_ds_root = True
        rhs = m.group("rhs")
        if "TPCAS_DIR" in rhs and "vendor/ds" in rhs:
            out.append(f'{m.group("indent")}DS_ROOT ?=')
            continue
    out.append(line)

text = "\n".join(out) + ("\n" if original.endswith("\n") else "")

# 2. Rewrite all remaining hardcoded vendored-DS paths.
replacements = {
    "$(TPCAS_DIR)/vendor/ds/lib": "$(DS_ROOT)/lib",
    "${TPCAS_DIR}/vendor/ds/lib": "$(DS_ROOT)/lib",
    "$(TPCAS_DIR)/vendor/ds": "$(DS_ROOT)",
    "${TPCAS_DIR}/vendor/ds": "$(DS_ROOT)",
}
for old, new in replacements.items():
    text = text.replace(old, new)

# 3. Ensure there is a non-overriding declaration for DS_ROOT.
if not re.search(r'(?m)^\s*DS_ROOT\s*\?=', text):
    tpcas_match = re.search(r'(?m)^(\s*TPCAS_DIR\s*\?=.*)$', text)
    if tpcas_match:
        insert_at = tpcas_match.end()
        text = text[:insert_at] + "\nDS_ROOT    ?=" + text[insert_at:]
    else:
        # Safe fallback near the top.
        text = "DS_ROOT ?=\n" + text

# 4. Normalize DS_DIR if present in a stale form.
text = re.sub(
    r'(?m)^\s*DS_DIR\s*[:?+]?=.*$',
    'DS_DIR     := $(DS_ROOT)/lib',
    text,
    count=1,
)

path.write_text(text)

if text == original:
    print("Makefile already had no stale vendored-DS paths.")
else:
    print("Repaired Makefile.")
PY

echo
echo "Relevant Make variables after repair:"
grep -nE '^(TPCAS_DIR|DS_ROOT|DS_DIR|CPPFLAGS)[[:space:]]*[:?+]?=|vendor/ds' Makefile || true

if grep -q 'vendor/ds' Makefile; then
  echo
  echo "warning: Makefile still contains vendor/ds references; inspect the lines above." >&2
  exit 2
fi

echo
echo "Backup: Makefile.before-ds-root-repair"
echo "Next:"
echo '  nix develop --command bash -lc '"'"''
echo '    echo "DS_ROOT=$DS_ROOT"'
echo '    test -f "$DS_ROOT/lib/context.h"'
echo '    make clean'
echo '    make DS_ROOT="$DS_ROOT"'
echo "  '"
