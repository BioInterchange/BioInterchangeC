#!/usr/bin/env bash

set -e

cmake -DCMAKE_BUILD_TYPE=Release -G 'Unix Makefiles'

if [[ "`uname`" == "Darwin" ]] ; then
    make openssl-lib
fi

make curl-lib
make libdocument-lib

make biointerchange

strip biointerchange

