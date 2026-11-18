# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
#
# Locate swipl and its C headers via `swipl --dump-runtime-variables`.
find_program(SWIPL_EXECUTABLE swipl REQUIRED)
execute_process(COMMAND ${SWIPL_EXECUTABLE} --dump-runtime-variables
                OUTPUT_VARIABLE _swipl_vars RESULT_VARIABLE _swipl_rc)
if(NOT _swipl_rc EQUAL 0)
  message(FATAL_ERROR "swipl --dump-runtime-variables failed")
endif()
string(REGEX MATCH "PLBASE=\"([^\"]*)\"" _m "${_swipl_vars}")
set(SWIPL_HOME "${CMAKE_MATCH_1}")
find_path(SWIPL_INCLUDE_DIR SWI-Prolog.h
          HINTS "${SWIPL_HOME}/include" REQUIRED)
message(STATUS "SWI-Prolog: ${SWIPL_EXECUTABLE} (headers in ${SWIPL_INCLUDE_DIR})")
