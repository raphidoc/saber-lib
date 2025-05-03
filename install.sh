#!/usr/bin/env bash
# install.sh — build & install saber-lib
set -euo pipefail

# Allow override of install prefix
: "${PREFIX:=/usr}"

echo "==> Building saber-lib (prefix=${PREFIX})"
BUILD_DIR=build
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}"

echo "==> Compiling"
make -j"$(nproc)"

echo "==> Running tests"
ctest --output-on-failure

echo "==> Installing to ${PREFIX}"
# If you need sudo, user can re-invoke with sudo ./install.sh
make install

echo "==> Verifying pkg-config"
if ! pkg-config --exists saber-lib; then
  echo "ERROR: pkg-config cannot see saber-lib.pc"
  exit 1
fi

echo "Installed saber-lib $(pkg-config --modversion saber-lib) successfully!"
