#!/usr/bin/env bash

if [[ ! -f CMakeLists.txt || ! -d src ]] ; then
    echo "Are you in the correct directory? Only execute this script as ./clean.sh!"
    exit 1
fi

rm -rf CMakeCache.txt \
       CTestTestfile.cmake \
       CMakeFiles \
       Makefile \
       libdocument-lib-prefix \
       cmake_install.cmake \
       doc/html \
       doc/Doxyfile \
       doxygen-prefix \
       biointerchange \
       biointerchange.dylib \
       biointerchange.build \
       biointerchange.xcodeproj \
       libdocument-lib-prefix \
       googletest-lib-prefix \
       CMakeScripts \
       src/doxygen \
       src/googletest \
       src/libdocument \
       Debug

