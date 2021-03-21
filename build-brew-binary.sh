#!/usr/bin/env bash

set -e

version="`./biointerchange -v`"

mkdir -p brew/biointerchange

cp biointerchange brew/biointerchange
cp LICENSE.txt brew/biointerchange/license.txt

cd brew

if [[ -f "biointerchange-$version.tar.gz" ]] ; then
    rm -f "biointerchange-$version.tar.gz"
fi

tar cfz "biointerchange-$version.tar.gz" biointerchange

cd ..

shasum -a 256 "brew/biointerchange-$version.tar.gz"

