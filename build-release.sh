#!/usr/bin/env bash

set -e

cmake -DCMAKE_BUILD_TYPE=Release -G 'Unix Makefiles'

make libdocument-lib

make biointerchange

strip biointerchange

