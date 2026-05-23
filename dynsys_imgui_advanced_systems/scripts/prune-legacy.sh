#!/usr/bin/env sh
set -eu

apply=0
case "${1:-}" in
  --apply) apply=1 ;;
  "") ;;
  *)
    echo "usage: $0 [--apply]" >&2
    exit 2
    ;;
esac

remove_paths="
src/dynsys.c
dynsys.c
FSEX302.ttf
result.png
result
build
dynsys
dynsys_imgui.zip
dynsys_dear_imgui.zip
"

if [ "$apply" -eq 0 ]; then
  echo "dry run: legacy files/directories that would be removed:"
else
  echo "removing legacy files/directories:"
fi

found=0
for p in $remove_paths; do
  if [ -e "$p" ]; then
    found=1
    echo "  $p"
    if [ "$apply" -eq 1 ]; then
      rm -rf -- "$p"
    fi
  fi
done

if [ "$found" -eq 0 ]; then
  echo "  none"
fi

if [ "$apply" -eq 0 ]; then
  echo
  echo "run '$0 --apply' to actually remove them."
fi
