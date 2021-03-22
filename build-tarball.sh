#!/usr/bin/env bash

set -e

version=`./biointerchange -v | sed -E 's/\+.*//'`

# Note: Does not copy GFF, GVF, VCF, etc., examples.
tar czf biointerchange-${version}.tar.gz \
    xcode.sh \
    clean.sh \
    build-release.sh \
    build-deb.sh \
    build-brew-binary.sh \
    build-tarball.sh \
    README.md \
    LICENSE.txt \
    CMakeLists.txt \
    test \
    test-data/simplepy \
    templates \
    context \
    scripts \
    lib \
    include \
    src/biointerchange

