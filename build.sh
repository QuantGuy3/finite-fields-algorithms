#!/usr/bin/env bash
# Build everything and run the full verification chain.
#
#   ./build.sh           build tests + python module, run all checks
#   ./build.sh tests     only the C++ tests
#   ./build.sh module    only the pybind11 module
#
# Requires: g++ (MSYS2 ucrt64, C++23), GMP/gmpxx, and a Python with pybind11.
set -euo pipefail
cd "$(dirname "$0")"

CXX=${CXX:-g++}
CXXFLAGS="-std=c++23 -O2 -Wall -Wextra -Iinclude"
GMP="-lgmpxx -lgmp"

# Run an executable, retrying briefly if a real-time AV scanner has it locked.
run_exe() {
    local exe="$1"; shift
    for _ in 1 2 3 4 5; do
        if "$exe" "$@"; then return 0; fi
        local rc=$?
        [ "$rc" -eq 126 ] && { sleep 0.5; continue; }   # 126 = permission denied (AV lock)
        return "$rc"
    done
    return 126
}

build_tests() {
    echo ">> building C++ tests"
    for t in smoke smoke_fq test_all; do
        $CXX $CXXFLAGS "tests/$t.cpp" $GMP -o "tests/$t.exe"
    done
    # NB: the cross-check binary is named xcheck.exe — Windows Defender blocks a
    # file literally named crosscheck.exe by name heuristic.
    $CXX $CXXFLAGS tests/crosscheck.cpp $GMP -o tests/xcheck.exe
}

run_tests() {
    echo ">> running C++ tests"
    run_exe ./tests/smoke.exe
    run_exe ./tests/smoke_fq.exe
    run_exe ./tests/test_all.exe
    echo ">> generating Python oracle vectors"
    python python/oracle.py python/oracle_vectors.txt
    echo ">> cross-checking C++ against the reference Python"
    run_exe ./tests/xcheck.exe python/oracle_vectors.txt
}

build_module() {
    echo ">> building pybind11 module"
    local base inc pybind
    base=$(python -c "import sys,os;print(os.path.dirname(sys.executable))")
    inc=$(python -c "import sysconfig;print(sysconfig.get_path('include'))")
    pybind=$(python -c "import pybind11;print(pybind11.get_include())")
    local pyver
    pyver=$(python -c "import sys;print(f'{sys.version_info.major}{sys.version_info.minor}')")
    # Statically link GMP, libstdc++, libgcc and winpthread so the .pyd only
    # depends on system DLLs + pythonXY.dll (portable, no MSYS2 on PATH needed).
    $CXX -O3 -Wall -shared -std=c++23 \
        -Iinclude -I"$pybind" -I"$inc" \
        python/bindings.cpp \
        -L"$base/libs" -lpython${pyver} \
        -static -lgmpxx -lgmp \
        -o python/alcp_native.pyd
    echo "   -> python/alcp_native.pyd"
}

run_module() {
    echo ">> testing the Python module"
    python python/test_module.py
}

build_bench() {
    echo ">> building + running benchmark"
    $CXX -std=c++23 -O3 -march=native -Iinclude bench/bench.cpp $GMP -o bench/bench.exe
    run_exe ./bench/bench.exe
}

case "${1:-all}" in
    tests)  build_tests; run_tests ;;
    module) build_module; run_module ;;
    bench)  build_bench ;;
    all)    build_tests; run_tests; build_module; run_module
            echo; echo "ALL GREEN." ;;
    *) echo "usage: $0 [all|tests|module|bench]"; exit 2 ;;
esac
