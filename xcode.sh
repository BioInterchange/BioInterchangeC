#!/usr/bin/env bash

if [[ ! -d biointerchange.xcodeproj ]] ; then
    cmake -G Xcode
fi

/Applications/Xcode.app/Contents/MacOS/Xcode

