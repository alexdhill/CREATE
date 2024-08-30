#!/bin/bash

docker run build \
    --progress plain \
    -t alexdhill/create:${1} \
    modules/${1} \
&& docker push alexdhill/create:${1}