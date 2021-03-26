#!/usr/bin/env bash

set -e

git rev-list --all \
    | wc -l \
    | tr -d ' '

