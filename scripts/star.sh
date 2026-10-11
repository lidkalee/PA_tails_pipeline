#!/usr/bin/env bash
set -euo pipefail
if [[ "$(uname -s)" == Darwin && "$(uname -m)" == arm64 ]]; then
  root=$(cd "$(dirname "$0")/.." && pwd)
  if [[ ! -x "$root/tools/STAR" ]]; then
    echo "Apple Silicon STAR 2.7.11b may silently read zero records. Run: bash scripts/install_star_macos.sh" >&2
    exit 1
  fi
  exec "$root/tools/STAR" "$@"
fi
exec STAR "$@"
