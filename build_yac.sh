#!/bin/bash

set -u
set -e
set -x

prefix=$HOME/local/dkrz

build_yaxt () {

git clone -b release-0.10.0 https://gitlab.dkrz.de/dkrz-sw/yaxt.git || (cd yaxt && git pull)

pushd yaxt

autoreconf -i

./configure --prefix=${prefix}

make all
make install

popd
}

build_yac () {

git clone -b release-3.0.2_p1 https://gitlab.dkrz.de/dkrz-sw/yac.git || (cd yac && git pull)

pushd yac

autoreconf -i

netcdf=$HOME/local/netcdf

export CC=mpicc FC=mpifort CFLAGS="-O3 -g -march=native"
./configure --prefix=${prefix} \
            --with-yaxt-root=${prefix} \
            --disable-netcdf

make all
make install

popd
}

build_yaxt
build_yac
