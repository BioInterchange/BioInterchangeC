#!/usr/bin/env bash

set -e

version=`./biointerchange -v`

tar czf biointerchange-${version}.tar.gz \
    xcode.sh \
    clean.sh \
    build-release.sh \
    build-deb.sh \
    build-brew.sh \
    build-tarball.sh \
    README.md \
    LICENSE.txt \
    CMakeLists.txt \
    test \
    test-data/simplepy \
    templates \
    examples \
    context \
    scripts \
    lib \
    include \
    src/biointerchange

