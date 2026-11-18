# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
cmake_minimum_required(VERSION 3.16)
project(sage1 VERSION 1.0.0 LANGUAGES C)

option(SAGE1_SANITIZE "Build everything with AddressSanitizer + UBSan" OFF)
option(SAGE1_PROLOG   "Build the SWI-Prolog foreign library and grammar tests" OFF)
option(SAGE1_FUTHARK  "Build the Futhark array backend and its cross-check" OFF)
set(SAGE1_FUTHARK_BACKENDS "c;multicore" CACHE STRING
    "Futhark backends to build and cross-check (c, multicore, opencl, cuda)")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SAGE1_WARNINGS -Wall -Wextra -Werror -pedantic)
if(SAGE1_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer
                      -fno-sanitize-recover=undefined)
  add_link_options(-fsanitize=address,undefined)
endif()

find_library(MATH_LIB m)

# ---- C11 kernel (original sources, compiled as-is) ----
add_library(sage1 STATIC native/sage1/sage1.c)
target_include_directories(sage1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/native/sage1)
target_compile_options(sage1 PRIVATE ${SAGE1_WARNINGS})
if(MATH_LIB)
  target_link_libraries(sage1 PUBLIC ${MATH_LIB})
endif()

enable_testing()

add_executable(test_sage1 native/sage1/test_sage1.c)
target_compile_options(test_sage1 PRIVATE ${SAGE1_WARNINGS})
target_link_libraries(test_sage1 PRIVATE sage1)
# Its checks are assert()s: keep them live in Release builds too.
target_compile_options(test_sage1 PRIVATE -UNDEBUG)
add_test(NAME c.original COMMAND test_sage1)

add_executable(test_sage1_extended tests/c/test_sage1_extended.c)
target_compile_options(test_sage1_extended PRIVATE ${SAGE1_WARNINGS})
target_link_libraries(test_sage1_extended PRIVATE sage1)
add_test(NAME c.extended COMMAND test_sage1_extended)

# Exact build command from native/sage1/README.md, run from the repo root.
add_test(NAME c.readme_sanitizer
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/tests/c/readme_sanitizer.sh
                 ${CMAKE_C_COMPILER} ${CMAKE_CURRENT_BINARY_DIR}/uor-sage1-test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(SAGE1_PROLOG)
  include(cmake/SWIProlog.cmake)
  add_subdirectory(prolog)
endif()

if(SAGE1_FUTHARK)
  add_subdirectory(futhark)
endif()
