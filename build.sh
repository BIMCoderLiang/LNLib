#!/usr/bin/env bash
#
# LNLib one-click build script for Linux.
#
set -euo pipefail

usage() {
    cat <<'EOF'
LNLib Linux build script

Usage: ./build.sh [options]

Options:
    --static    Build the static library only (libLNLib.a)
    --debug     Build in Debug mode instead of Release
    --help      Show this help message

Default (no options): build BOTH the shared library (libLNLib.so)
and the static library (libLNLib.a) in Release mode.

Output:
    Libraries are placed in build/lib/ and installed to the project's
    install/ directory:
        install/lib/libLNLib.so | install/lib/libLNLib.a
        install/include/               (public headers)
        install/lib/cmake/LNLib/       (CMake package: find_package(LNLib),
                                        link target LNLib::LNLib)

Note:
    Unit test support (--test) is planned as a separate task.
EOF
}

BUILD_SHARED=ON
BUILD_STATIC=ON
BUILD_TYPE=Release

while [ $# -gt 0 ]; do
    case "$1" in
        --static) BUILD_SHARED=OFF; BUILD_STATIC=ON ;;
        --debug) BUILD_TYPE=Debug ;;
        --help) usage; exit 0 ;;
        *)
            echo "Error: unknown option '$1' (use --help for usage)" >&2
            exit 1
            ;;
    esac
    shift
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INSTALL_DIR="${SCRIPT_DIR}/install"

build_lib() {
    local shared=$1
    # Reconfigure only when the cached options differ; cmake --build itself
    # regenerates when CMakeLists change, so a skipped configure is safe.
    local cached_shared cached_type
    cached_shared=$(sed -n 's/^LNLib_BUILD_SHARED:BOOL=//p' build/CMakeCache.txt 2>/dev/null | head -1)
    cached_type=$(sed -n 's/^CMAKE_BUILD_TYPE:STRING=//p' build/CMakeCache.txt 2>/dev/null | head -1)
    if [ "${cached_shared:-}" != "${shared}" ] || [ "${cached_type:-}" != "${BUILD_TYPE}" ]; then
        cmake -S . -B build \
            -DLNLib_BUILD_SHARED=${shared} \
            -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
            -DENABLE_UNIT_TESTS=OFF
    fi
    cmake --build build --target LNLib -j"$(nproc)"
    # Install right after each form build: the CMake cache only keeps the last
    # LNLib_BUILD_SHARED configuration, so a single install at the end would
    # silently miss the library form that was built first.
    cmake --install build --prefix "${INSTALL_DIR}"
}

if [ "${BUILD_SHARED}" = "ON" ]; then
    build_lib ON
fi
if [ "${BUILD_STATIC}" = "ON" ]; then
    build_lib OFF
fi

# Remove artifacts of library forms that were not built this run,
# so build/lib/ and install/lib/ never mix stale output of the other form.
if [ "${BUILD_SHARED}" = "OFF" ]; then
    rm -f build/lib/libLNLib.so install/lib/libLNLib.so
fi
if [ "${BUILD_STATIC}" = "OFF" ]; then
    rm -f build/lib/libLNLib.a install/lib/libLNLib.a
fi

echo
echo "Build finished."
if [ "${BUILD_SHARED}" = "ON" ]; then
    echo "  Shared: build/lib/libLNLib.so -> install/lib/libLNLib.so (${BUILD_TYPE})"
fi
if [ "${BUILD_STATIC}" = "ON" ]; then
    echo "  Static: build/lib/libLNLib.a -> install/lib/libLNLib.a (${BUILD_TYPE})"
fi
echo "  Headers: install/include/ | CMake package: install/lib/cmake/LNLib/ (find_package(LNLib), target LNLib::LNLib)"
