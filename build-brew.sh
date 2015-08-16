#!/usr/bin/env bash

set -e

license_version=2
version="`./biointerchange -v`"

mkdir -p brew/biointerchange

cp biointerchange brew/biointerchange
cp biointerchange-l$license_version.txt brew/biointerchange/license.txt

cd brew

if [[ -f "biointerchange-$version.tar.gz" ]] ; then
    rm -f "biointerchange-$version.tar.gz"
fi

tar cfz "biointerchange-$version.tar.gz" biointerchange

cd ..

shasum "brew/biointerchange-$version.tar.gz"

