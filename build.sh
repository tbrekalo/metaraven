#!/bin/bash
# build procedure for metaraven

mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make

./bin/metaraven