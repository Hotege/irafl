#!/bin/bash

rm build -rf

mkdir build

pushd build
cmake .. \
-DCMAKE_INSTALL_PREFIX=install
popd
