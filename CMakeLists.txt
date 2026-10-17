# CMake para la biblioteca header-only y sus tests C++ (CTest).
#
# El módulo pybind11 NO se compila aquí: usa «python build.py module», que
# resuelve los detalles de g++ (MinGW) frente a CPython (MSVC) en Windows.
# build.py es también la vía principal para los tests («python build.py tests»);
# este archivo es la alternativa para IDE y para quien prefiera CMake/CTest.
#
#   cmake -S . -B build/cmake -G Ninja
#   cmake --build build/cmake
#   ctest --test-dir build/cmake --output-on-failure
#
# Se usa build/cmake (y no build/) para no mezclarse con lo que genera build.py.
cmake_minimum_required(VERSION 3.20)
project(alcp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Biblioteca de interfaz (solo cabeceras).
add_library(alcp INTERFACE)
target_include_directories(alcp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# GMP / gmpxx (no hay módulo estándar de CMake; se buscan por nombre).  El
# directorio lib/ de la toolchain se deduce de la ruta del compilador para que
# se encuentren las bibliotecas de MSYS2 ucrt64.
get_filename_component(_cxx_bin ${CMAKE_CXX_COMPILER} DIRECTORY)
get_filename_component(_toolchain_root ${_cxx_bin} DIRECTORY)
if(MINGW)
  # En MinGW los tests se enlazan estáticamente (ver más abajo): se buscan las
  # bibliotecas .a y no las de importación .dll.a.
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
endif()
find_library(GMP_LIB gmp REQUIRED HINTS ${_toolchain_root}/lib)
find_library(GMPXX_LIB gmpxx REQUIRED HINTS ${_toolchain_root}/lib)
target_link_libraries(alcp INTERFACE ${GMPXX_LIB} ${GMP_LIB})

enable_testing()

# Un ejecutable y un test de CTest por cada tests/*.cpp (check.hpp es el arnés
# CHECK, no un test).  CONFIGURE_DEPENDS: un test nuevo se detecta al compilar
# sin tener que reconfigurar a mano.
file(GLOB ALCP_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
foreach(src IN LISTS ALCP_TEST_SOURCES)
  get_filename_component(name ${src} NAME_WE)
  # crosscheck se compila como «xcheck»: Windows Defender bloquea por su
  # nombre un ejecutable llamado crosscheck.exe.
  if(name STREQUAL "crosscheck")
    set(name xcheck)
  endif()
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE alcp)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # -UNDEBUG: los tests comprueban siempre, también en Release.
    target_compile_options(${name} PRIVATE -Wall -Wextra -UNDEBUG)
  endif()
  if(MINGW)
    # Enlace estático en MinGW: el ejecutable no depende del PATH (en Git Bash,
    # la libstdc++-6.dll de /mingw64/bin es incompatible con la de ucrt64).
    target_link_options(${name} PRIVATE -static)
  endif()
  if(NOT name STREQUAL "xcheck")
    add_test(NAME ${name} COMMAND ${name})
  endif()
endforeach()

# xcheck compara el C++ con las prácticas Python originales usando los vectores
# que genera python-bindings/oracle.py.  Si hay intérprete de Python, un test
# «oracle» (fixture) los regenera antes de xcheck; si no lo hay, xcheck solo se
# registra cuando el archivo ya existe (p. ej., tras «python build.py oracle»).
set(ALCP_ORACLE_VECTORS "${CMAKE_CURRENT_SOURCE_DIR}/build/oracle_vectors.txt"
    CACHE FILEPATH "Vectores del oráculo para xcheck")
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  add_test(NAME oracle
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python-bindings/oracle.py
                   ${ALCP_ORACLE_VECTORS})
  set_tests_properties(oracle PROPERTIES FIXTURES_SETUP oracle_vectors
                       ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
  add_test(NAME xcheck COMMAND xcheck ${ALCP_ORACLE_VECTORS})
  set_tests_properties(xcheck PROPERTIES FIXTURES_REQUIRED oracle_vectors)
elseif(EXISTS ${ALCP_ORACLE_VECTORS})
  add_test(NAME xcheck COMMAND xcheck ${ALCP_ORACLE_VECTORS})
else()
  message(STATUS "xcheck no se registra en CTest: no hay Python ni existe "
                 "${ALCP_ORACLE_VECTORS} (genéralo con «python build.py oracle» y reconfigura)")
endif()
