# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
#
# Futhark array backend: `futhark <backend> --library` for every backend in
# SAGE1_FUTHARK_BACKENDS, each linked into a cross-check against the C kernel.
find_program(FUTHARK_EXECUTABLE futhark REQUIRED)
find_package(Threads)

set(SAGE1_FUT_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sage1.fut)

foreach(backend IN LISTS SAGE1_FUTHARK_BACKENDS)
  set(out_dir ${CMAKE_CURRENT_BINARY_DIR}/${backend})
  set(stem ${out_dir}/sage1_${backend})
  add_custom_command(
    OUTPUT ${stem}.c ${stem}.h
    COMMAND ${CMAKE_COMMAND} -E make_directory ${out_dir}
    COMMAND ${FUTHARK_EXECUTABLE} ${backend} --library -o ${stem} ${SAGE1_FUT_SOURCE}
    DEPENDS ${SAGE1_FUT_SOURCE}
    COMMENT "futhark ${backend} --library sage1.fut"
    VERBATIM)

  add_library(sage1_fut_${backend} STATIC ${stem}.c)
  target_include_directories(sage1_fut_${backend} PUBLIC ${out_dir})
  target_compile_options(sage1_fut_${backend} PRIVATE -O3)
  if(MATH_LIB)
    target_link_libraries(sage1_fut_${backend} PUBLIC ${MATH_LIB})
  endif()
  if(backend STREQUAL "multicore")
    target_link_libraries(sage1_fut_${backend} PUBLIC Threads::Threads)
  elseif(backend STREQUAL "opencl")
    find_package(OpenCL REQUIRED)
    target_link_libraries(sage1_fut_${backend} PUBLIC OpenCL::OpenCL)
  elseif(backend STREQUAL "cuda")
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(sage1_fut_${backend} PUBLIC CUDA::cuda_driver CUDA::nvrtc)
  endif()

  add_executable(test_futhark_${backend} test_futhark_vs_c.c)
  target_compile_definitions(test_futhark_${backend} PRIVATE
    FUT_BACKEND="${backend}" FUT_HEADER="sage1_${backend}.h")
  target_compile_options(test_futhark_${backend} PRIVATE -Wall -Wextra)
  target_link_libraries(test_futhark_${backend} PRIVATE sage1 sage1_fut_${backend})
  add_test(NAME futhark.${backend} COMMAND test_futhark_${backend})
endforeach()
