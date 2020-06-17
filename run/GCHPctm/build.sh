#!/bin/bash

source gchp.env
mkdir -P build
# Add option to build from scratch or rebuild
cd build
cmake ../CodeDir
make -j
make install
cd ..
