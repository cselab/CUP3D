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
CMAKE_VERSION_MAJOR=3.11  # Naming not correct really, but ok.

FFTW_VERSION=3.3.7
FFTW_MD5='0d5915d7d39b3253c1cc05030d79ac47  fftw-3.3.7.tar.gz'

HDF5_VERSION=1.10.1
HDF5_URL='https://www.hdfgroup.org/package/source-cmake-unix/?wpdmdl=4411&refresh=5afda694193441526572692'
HDF5_MD5='9121cf2f4d155867f61ea2f75a5f9b6d  CMake-hdf5-1.10.1.tar.gz'

echo "This script installs the following CubismUP_3D dependencies:"
echo "  - CMake ${CMAKE_VERSION} (required at least 3.2)"
echo "  - FFTW ${FFTW_VERSION}"
echo "  - HDF5 ${HDF5_VERSION}"
echo
echo "Target folder for sources:
echo
echo "

BASEPWD=$PWD

echo "Installing CMake ${CMAKE_VERSION}..."
wget -nc https://cmake.org/files/v${CMAKE_VERSION_MAJOR}/cmake-${CMAKE_VERSION}.tar.gz -P $SOURCES
cd $SOURCES
sha256sum --quiet -c - <<< $CMAKE_SHA_256
tar --skip-old-files -xzvf cmake-${CMAKE_VERSION}.tar.gz
cd cmake-${CMAKE_VERSION}
./bootstrap --parallel=${JOBS} --prefix=$INSTALL_PATH/cmake-${CMAKE_VERSION}/
make -j${JOBS}
make install -j${JOBS}
cd $BASEPWD
export PATH=$INSTALL_PATH/cmake-${CMAKE_VERSION}/bin:$PATH  # Required by HDF5!

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
cd CMake-hdf5-${HDF5_VERSION}
# FIXME: Parallel flags are not correct!!
ctest -DHDF5_ENABLE_PARALLEL:BOOL=ON -S HDF5config.cmake,BUILD_GENERATOR=Unix -C Release -V -O hdf5.log


echo "Done! Run or add to .bashrc the following commands to capture CMake, FFTW and HDF5 in your environment:"
echo "export PATH=$INSTALL_PATH/cmake-${CMAKE_VERSION}/bin:\$PATH"
echo "export FFTWDIR=$INSTALL_PATH/fftw-${FFTW_VERSION}/"
echo "export HDF5_ROOT=$SOURCES/CMake-hdf5-${HDF5_VERSION}/HDF_Group/HDF5/${HDF5_VERSION}"
