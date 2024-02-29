#!/bin/bash

set -u
set -e
set -x

prefix=$HOME/local/dkrz

build_yaxt () {

rm -rf yaxt
git clone -b release-0.10.0 https://gitlab.dkrz.de/dkrz-sw/yaxt.git
pushd yaxt

autoreconf -i

./configure --prefix=${prefix}

make all
make install

popd
}

build_yac () {

rm -rf yac
git clone -b release-3.1.1 https://gitlab.dkrz.de/dkrz-sw/yac.git
pushd yac

autoreconf -i

netcdf=$HOME/local/netcdf

export CC=mpicc FC=mpifort CFLAGS="-O3 -g -march=native"
./configure --prefix=${prefix} \
            --with-yaxt-root=${prefix} \
            --with-netcdf-root=${netcdf}

make all
make install

popd
}

build_yaxt
build_yac
