#!/bin/sh
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
#
# Runs the ASan/UBSan build command from native/sage1/README.md verbatim
# (compiler and output path are the only substitutions).
set -eu
CC_BIN=${1:-cc}
OUT=${2:-/tmp/uor-sage1-test}
"$CC_BIN" -std=c11 -O2 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  native/sage1/sage1.c native/sage1/test_sage1.c -lm -o "$OUT"
"$OUT"
