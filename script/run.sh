#!/bin/bash

set -x

#touch ../build && rm ../build -rf
mkdir -p ../build
cd ../build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

./src/virgo_plus_run ../data/SHA256_64.pws

cd ../script