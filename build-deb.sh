#!/usr/bin/env bash

set -e

mkdir -p deb/DEBIAN
mkdir -p deb/usr/bin

cp biointerchange deb/usr/bin

version=`./biointerchange -v`
today=`date`

<templates/control sed 's/VERSION/'$version'/' > deb/DEBIAN/control
<templates/changelog sed 's/VERSION/'$version'/' | \
    sed 's/DATE/'"$today"'/' > deb/DEBIAN/changelog

vim deb/DEBIAN/changelog

cat deb/DEBIAN/changelog >> templates/changelog

dpkg-deb -b deb biointerchange_${version}_amd64.deb

mkdir -p debs/amd64

mv biointerchange_${version}_amd64.deb debs/amd64

cd debs

dpkg-scanpackages amd64 | gzip -9c > amd64/Packages.gz

cd ..


