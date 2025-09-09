#!/bin/bash

for build in `find modules -mindepth 1 -type d`; do
    dockerfile="${build}/Dockerfile"
    name="$(basename ${build})"
    version="$(grep -m1 '_VER' ${dockerfile} | cut -d' ' -f3)"
    echo "Building alexdhill/create:${name}-${version}..."
    docker build -t alexdhill/create:${name}-${version} ${build} \
        > logs/${name}-${version}.log 2>&1 &
done