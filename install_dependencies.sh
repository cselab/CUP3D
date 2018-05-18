#!/bin/bash

set -e

# https://stackoverflow.com/a/23378780/2203044
# LOGICAL_CPU_COUNT=$([[ $(uname) = 'Darwin' ]] && sysctl -n hw.logicalcpu_max || lscpu -p | egrep -v '^#' | wc -l)
PHYSICAL_CPU_COUNT=$([[ $(uname) = 'Darwin' ]] && sysctl -n hw.physicalcpu_max || lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)

JOBS=${JOBS:-$PHYSICAL_CPU_COUNT}

SOURCES=${SOURCES:-$PWD/dependencies/}
INSTALL_PATH=${INSTALL_PATH:-$PWD/dependencies/build/}

# Shorthands for version. Note that changing this number may not be enough for
# the script to work properly!
CMAKE_VERSION=3.11.1
CMAKE_SHA_256='57bebc6ca4d1d42c6385249d148d9216087e0fda57a47dc5c858790a70217d0c  cmake-3.11.1.tar.gz'
CMAKE_VERSION_SHORT=3.11

FFTW_VERSION=3.3.7
FFTW_MD5='0d5915d7d39b3253c1cc05030d79ac47  fftw-3.3.7.tar.gz'

HDF5_VERSION=1.10.1
HDF5_URL='https://www.hdfgroup.org/package/source-cmake-unix/?wpdmdl=4411&refresh=5afda694193441526572692'
HDF5_MD5='9121cf2f4d155867f61ea2f75a5f9b6d  CMake-hdf5-1.10.1.tar.gz'

GSL_VERSION=2.1

# These echos are kind of useles...
echo "This script installs locally the following CubismUP_3D dependencies:"
echo "  - CMake ${CMAKE_VERSION} (required at least 3.2)"
echo "  - FFTW ${FFTW_VERSION}"
echo "  - HDF5 ${HDF5_VERSION}"
echo "  - GSL ${GSL_VERSION}"
echo

BASEPWD=$PWD

echo "Installing CMake ${CMAKE_VERSION}..."
wget -nc https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}.tar.gz -P $SOURCES
cd $SOURCES
sha256sum --quiet -c - <<< $CMAKE_SHA_256
tar --skip-old-files -xzvf cmake-${CMAKE_VERSION}.tar.gz
cd cmake-${CMAKE_VERSION}
./bootstrap --parallel=${JOBS} --prefix=$INSTALL_PATH/cmake-${CMAKE_VERSION}/
make -j${JOBS}
make install -j${JOBS}
cd $BASEPWD


echo "Installing FFTW ${FFTW_VERSION}..."
wget -nc http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz -P $SOURCES
cd $SOURCES
md5sum --quiet -c - <<< $FFTW_MD5
tar --skip-old-files -xzvf fftw-${FFTW_VERSION}.tar.gz
cd fftw-${FFTW_VERSION}
./configure --prefix=$INSTALL_PATH/fftw-${FFTW_VERSION}/ --enable-mpi --enable-threads --enable-shared
make -j${JOBS}
make install -j${JOBS}
cd $BASEPWD


echo "Installing HDF5 ${HFD5_VERSION}..."
wget ${HDF5_URL} -O $SOURCES/CMake-hdf5-${HDF5_VERSION}.tar.gz
cd $SOURCES
md5sum --quiet -c <<< $HDF5_MD5
tar --skip-old-files -xzvf CMake-hdf5-${HDF5_VERSION}.tar.gz
cd CMake-hdf5-${HDF5_VERSION}/hdf5-${HDF5_VERSION}
./configure --prefix=$INSTALL_PATH/hdf5-${HDF5_VERSION}-parallel/ --enable-parallel
make -j $JOBS
make install -j $JOBS
cd $BASEPWD


echo "Installing GSL ${GSL_VERSION}..."
wget -nc http://mirror.switch.ch/ftp/mirror/gnu/gsl/gsl-${GSL_VERSION}.tar.gz -P $SOURCES
cd $SOURCES
tar --skip-old-files -xzvf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}
./configure --prefix=$INSTALL_PATH/gsl-${HDF5_VERSION} --enable-parallel
make -j $JOBS
make install -j $JOBS


echo "======================================================================"
echo "Done! Run or add to ~/.bashrc the following commands"
echo "to capture CMake, FFTW and HDF5 in your environment:"
echo
echo "export PATH=$INSTALL_PATH/cmake-${CMAKE_VERSION}/bin:\$PATH"
echo "export FFTWDIR=$INSTALL_PATH/fftw-${FFTW_VERSION}/"
echo "export HDF5_ROOT=$INSTALL_PATH/hdf5-${HDF5_VERSION}-parallel/"
echo "export GSL_ROOT_DIR=$INSTALL_PATH/gsl-${GSL_VERSION}/"
