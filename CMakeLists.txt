cmake_minimum_required(VERSION 3.16)
project(MyMD
        VERSION 0.1.0
        DESCRIPTION "Molecular dynamics of ellipsoidal packings"
        LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Language and defaults
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# The source is C++17 but leans on a few GNU-isms inherited from its C++98
# days, so gnu++17 rather than c++17.
set(CMAKE_CXX_EXTENSIONS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---------------------------------------------------------------------------
# Dependencies
#
# GSL is used for exactly one thing: the 4x4 non-symmetric eigensolver in the
# ellipsoid-ellipsoid contact test.  'make gsl' builds a copy into .deps/gsl;
# point CMake at it if it is there, otherwise let find_package look in the
# usual places (Homebrew, /usr, ...).
# ---------------------------------------------------------------------------
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.deps/gsl/include/gsl")
  set(GSL_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.deps/gsl" CACHE PATH
      "Bundled GSL built by 'make gsl'")
endif()

find_package(GSL REQUIRED)

# ---------------------------------------------------------------------------
# The solver.  It is one translation unit by design: main.cc includes the
# whole thing (see docs/ARCHITECTURE.md).
# ---------------------------------------------------------------------------
add_executable(ellipmd main.cc)
target_include_directories(ellipmd PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  ${GSL_INCLUDE_DIRS})
target_link_libraries(ellipmd PRIVATE ${GSL_LIBRARIES})
target_compile_options(ellipmd PRIVATE -Wall)
# single source of truth for the version: CMakeLists.txt's project(VERSION)
target_compile_definitions(ellipmd PRIVATE
  MYMD_VERSION="${PROJECT_VERSION}")

# ---------------------------------------------------------------------------
# Analysis utilities.  See tools/Makefile for the ones that were dropped
# because they cannot be built outside the author's original setup.
# ---------------------------------------------------------------------------
set(MYMD_TOOLS coord_convert packing_density fabric compact periodic2full)

foreach(tool IN LISTS MYMD_TOOLS)
  add_executable(${tool} "tools/${tool}.cc")
  target_include_directories(${tool} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${GSL_INCLUDE_DIRS})
  target_link_libraries(${tool} PRIVATE ${GSL_LIBRARIES})
  # the tools include the solver headers by relative path
  set_target_properties(${tool} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endforeach()

set_target_properties(ellipmd PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# ---------------------------------------------------------------------------
# Tests
#
# check_physics.py compares against committed reference snapshots and checks
# invariants; --elf tells it which binary to run, so it works from any build
# directory.  See ROADMAP.md Phase 0.
# ---------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter)

if(Python3_Interpreter_FOUND)
  enable_testing()
  add_test(NAME physics
           COMMAND "${Python3_EXECUTABLE}"
                   "${CMAKE_CURRENT_SOURCE_DIR}/bench/check_physics.py"
                   --elf "$<TARGET_FILE:ellipmd>")
  set_tests_properties(physics PROPERTIES TIMEOUT 600)
endif()

# ---------------------------------------------------------------------------
# Optional sanitizer build
#
#   cmake -B build-asan -DMYMD_SANITIZE=address,undefined
#
# Clang is needed for this on macOS because Homebrew GCC ships no sanitizer
# runtime.  It is what found the DisBetaDistribution heap overflow; see
# CHANGELOG.md.
# ---------------------------------------------------------------------------
set(MYMD_SANITIZE "" CACHE STRING
    "Comma-separated sanitizers, e.g. address,undefined")

if(MYMD_SANITIZE)
  message(STATUS "Sanitizers enabled: ${MYMD_SANITIZE}")
  add_compile_options(-fsanitize=${MYMD_SANITIZE} -fno-omit-frame-pointer -g -O1)
  add_link_options(-fsanitize=${MYMD_SANITIZE})
endif()
