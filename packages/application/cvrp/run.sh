#!/bin/bash

rm -rf build bin
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..

make -j$(nproc)

cd ..

./bin/cvrp instances/A-n60-k9.vrp -u 1355
#./bin/cvrp instances/P-n101-k4.vrp -u 682