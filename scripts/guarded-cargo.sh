#!/bin/bash
set -euo pipefail
unset CPPFLAGS LDFLAGS
export CARGO_HOME=/Volumes/KIOXIA/Developments/cargo-home
export CARGO_TARGET_DIR=/Volumes/KIOXIA/Developments/cargo-target/rsomics-popgen
export TMPDIR=/Volumes/KIOXIA/Developments/tmp
export RUSTUP_HOME=/Volumes/KIOXIA/Developments/rustup-home
export XDG_CACHE_HOME=/Volumes/KIOXIA/Developments/tmp/cache
export GH_CACHE_DIR=/Volumes/KIOXIA/Developments/tmp/gh-cache
export PYTHONDONTWRITEBYTECODE=1
toolchain=/Volumes/KIOXIA/Developments/rustup-home/toolchains/1.91.0-aarch64-apple-darwin/bin
export CARGO="$toolchain/cargo" RUSTC="$toolchain/rustc" RUSTDOC="$toolchain/rustdoc" RUSTFMT="$toolchain/rustfmt"
export PATH="$toolchain:$PATH"
export CARGO_BUILD_JOBS=2
/usr/sbin/diskutil info -plist /Volumes/KIOXIA | /usr/bin/python3 -c '
import plistlib,sys
v=plistlib.loads(sys.stdin.buffer.read())
assert v.get("MountPoint")=="/Volumes/KIOXIA" and v.get("Internal") is False
print("verified external volume:",v["DeviceNode"])
'
external_device=$(stat -f %d /Volumes/KIOXIA)
for path in "$CARGO_HOME" "$CARGO_TARGET_DIR" "$TMPDIR" "$RUSTUP_HOME" "$toolchain" "$XDG_CACHE_HOME" "$GH_CACHE_DIR"; do
  resolved=$(cd "$path" && pwd -P)
  case "$resolved" in /Volumes/KIOXIA/Developments/*) ;; *) echo "Forbidden path: $resolved" >&2; exit 1;; esac
  test "$(stat -f %d "$resolved")" = "$external_device"
  echo "verified external path: $resolved"
done
for argument in "$@"; do
  case "$argument" in --target-dir|--target-dir=*|--config|--config=*) echo "Storage overrides are forbidden" >&2; exit 1;; esac
done
df -h /Volumes/KIOXIA
target=$("$CARGO" metadata --no-deps --format-version 1 | /usr/bin/python3 -c 'import json,sys; print(json.load(sys.stdin)["target_directory"])')
test "$target" = "$CARGO_TARGET_DIR"
echo "Cargo resolved target: $target"
if [[ "${1:-}" == "--exec" ]]; then shift; exec "$@"; fi
if [[ "${1:-}" == "clippy" ]]; then shift; exec "$toolchain/cargo-clippy" clippy "$@"; fi
if [[ "${1:-}" == "fmt" ]]; then shift; exec "$toolchain/cargo-fmt" fmt "$@"; fi
exec "$CARGO" "$@"
