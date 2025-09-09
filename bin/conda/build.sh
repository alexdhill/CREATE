#!/bin/bash

for env in `find modules -mindepth 1 -type f -name "*.yaml"`; do
    build="$(basename $env | sed 's/.yaml$//')"
    mamba env create -f $env \
        2>&1 > logs/$build.log
done