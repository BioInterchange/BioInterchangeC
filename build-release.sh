#!/usr/bin/env bash

set -e

cmake -DCMAKE_BUILD_TYPE=Release -G 'Unix Makefiles'

make python
make openssl-lib
make curl-lib
#make googletest-lib
make libdocument-lib

make biointerchange

strip biointerchange

