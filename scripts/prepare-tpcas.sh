#!/bin/sh
# Stage the migrated frontend without editing its source checkout.
set -eu
source_dir=$1
dest_dir=$2
patch_dir=$3

test -f "$source_dir/src/arena.h" || {
  echo "TPCAS_DIR must point to a TPCAS source checkout" >&2
  exit 1
}
grep -q 'dsa/core/context.h' "$source_dir/src/arena.h" || {
  echo "TPCAS requires the DSA migration. Apply the TPCAS update or refresh its flake input after pushing it." >&2
  exit 1
}
mkdir -p "$dest_dir/src"
rm -f "$dest_dir/.prepared"
cp "$source_dir"/src/*.c "$source_dir"/src/*.h "$dest_dir/src/"
# Do not ask patch to guess the direction: --batch can silently reverse it.
if ! grep -q 'Decimal/scientific constants used by .dyn models' "$dest_dir/src/lex.c"; then
  patch --batch --forward -d "$dest_dir" -p1 < "$patch_dir/tpcas-numbers.patch"
fi
touch "$dest_dir/.prepared"
