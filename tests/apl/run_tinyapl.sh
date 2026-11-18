#!/bin/sh
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
#
# run_tinyapl.sh <tinyapl> <engine.tinyapl> <checker> <workdir>
set -eu
OUT="$4/tinyapl_engine_out.txt"
"$1" "$2" > "$OUT" 2> "$4/tinyapl_engine_err.txt"
if [ -s "$4/tinyapl_engine_err.txt" ]; then cat "$4/tinyapl_engine_err.txt" >&2; exit 1; fi
cat "$OUT"
"$3" "$OUT"
