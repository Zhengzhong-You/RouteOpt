#!/bin/bash

rm -rf build bin

mkdir build

cd build

cmake -DCMAKE_BUILD_TYPE=Release ..

make -j$(nproc)

cd ..

./bin/cvrp instance/P-n20-k2.vrp -u 217
