#!/usr/bin/env bash

set -e

cmake -DCMAKE_BUILD_TYPE=Release -G 'Unix Makefiles'

make python

if [[ "`uname`" == "Darwin" ]] ; then
    make openssl-lib
fi

make curl-lib
#make googletest-lib
make libdocument-lib

make biointerchange

strip biointerchange

