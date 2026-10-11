#!/usr/bin/env bash
# Build the upstream Apple Silicon read-buffer fix from a pinned commit.
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "$(uname -s)" != Darwin || "$(uname -m)" != arm64 ]]; then
  echo "This installer is for Apple Silicon macOS." >&2
  exit 1
fi
if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "Activate the pa-tails environment first." >&2
  exit 1
fi
compiler=$(find "$CONDA_PREFIX/bin" -maxdepth 1 -name 'arm64-apple-darwin*-clang++' -print -quit)
if [[ -z "$compiler" ]]; then
  echo "Conda Clang++ is missing; install clangxx_osx-arm64 in this environment." >&2
  exit 1
fi
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT
git clone --depth 1 --branch fix/macos-libcxx-0reads \
  https://github.com/BenjaminDEMAILLE/STAR.git "$workdir/STAR"
expected=14b1f235927e75ff927be08134c1fb6a00c14d79
actual=$(git -C "$workdir/STAR" rev-parse HEAD)
if [[ "$actual" != "$expected" ]]; then
  echo "STAR source changed; expected $expected, got $actual" >&2
  exit 1
fi
make -C "$workdir/STAR/source" -j 4 STARforMacStatic \
  CXX="$compiler" CXXFLAGS_SIMD= \
  'LDFLAGS_Mac_static=-pthread -lz htslib/libhts.a' BUILD_DATE="$(date -u +%Y-%m-%d)" \
  > "$workdir/build.log" 2>&1 || { tail -50 "$workdir/build.log" >&2; exit 1; }
mkdir -p tools
install -m 755 "$workdir/STAR/source/STAR" tools/STAR
echo "Installed patched STAR at tools/STAR"
