# CMake build for the C++ tests.  The pybind11 module is best built with
# build.sh (which handles the MinGW <-> MSVC-Python details); this file targets
# the header-only library and its test executables.
#
#   cmake -S . -B build -G Ninja
#   cmake --build build
#   ctest --test-dir build --output-on-failure
cmake_minimum_required(VERSION 3.20)
project(alcp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Header-only interface library.
add_library(alcp INTERFACE)
target_include_directories(alcp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# GMP / gmpxx (no standard CMake module; locate by name).  Derive the toolchain
# lib/ directory from the compiler location so MSYS2's ucrt64 libs are found.
get_filename_component(_cxx_bin ${CMAKE_CXX_COMPILER} DIRECTORY)
get_filename_component(_toolchain_root ${_cxx_bin} DIRECTORY)
find_library(GMP_LIB gmp REQUIRED HINTS ${_toolchain_root}/lib)
find_library(GMPXX_LIB gmpxx REQUIRED HINTS ${_toolchain_root}/lib)
target_link_libraries(alcp INTERFACE ${GMPXX_LIB} ${GMP_LIB})

enable_testing()
# crosscheck is built as "xcheck" — Windows Defender blocks a file named
# literally crosscheck.exe by name heuristic.
foreach(pair smoke:smoke smoke_fq:smoke_fq test_all:test_all xcheck:crosscheck)
  string(REPLACE ":" ";" parts ${pair})
  list(GET parts 0 exe)
  list(GET parts 1 src)
  add_executable(${exe} tests/${src}.cpp)
  target_link_libraries(${exe} PRIVATE alcp)
endforeach()

add_test(NAME smoke      COMMAND smoke)
add_test(NAME smoke_fq   COMMAND smoke_fq)
add_test(NAME test_all   COMMAND test_all)
# (xcheck needs python/oracle_vectors.txt; generate it via build.sh or
#  `python python/oracle.py python/oracle_vectors.txt` before running.)
