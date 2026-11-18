# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab

# Foreign library loaded by swipl; libswipl symbols resolve from the host.
add_library(sage1_pl MODULE sage1_pl.c)
set_target_properties(sage1_pl PROPERTIES PREFIX "")
target_include_directories(sage1_pl PRIVATE ${SWIPL_INCLUDE_DIR})
target_compile_options(sage1_pl PRIVATE -Wall -Wextra -Werror)
target_link_libraries(sage1_pl PRIVATE sage1)
if(SAGE1_SANITIZE)
  # An ASan runtime cannot be injected into an unsanitized swipl process.
  target_compile_options(sage1_pl PRIVATE -fno-sanitize=all)
  target_link_options(sage1_pl PRIVATE -fno-sanitize=all)
endif()

add_test(NAME prolog.grammar
         COMMAND ${SWIPL_EXECUTABLE} -q -g run_tests -t halt
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_sage1_grammar.pl)
set_tests_properties(prolog.grammar PROPERTIES
  ENVIRONMENT "SAGE1_PL_LIB=$<TARGET_FILE:sage1_pl>")
