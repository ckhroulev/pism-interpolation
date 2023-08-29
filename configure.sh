#!/bin/bash

set -x
set -u
set -e

YAC=$HOME/local/dkrz
PISM=$HOME/local/pism
NETCDF=$HOME/local/netcdf

prefix=$HOME/local/pism-interp

cmake -B build -S . \
      -DCMAKE_COLOR_MAKEFILE=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -DCMAKE_INSTALL_PREFIX=${prefix} \
      -DCMAKE_PREFIX_PATH="${YAC};${PISM};${NETCDF}"

rm -f compile_commands.json
ln -s build/compile_commands.json .
