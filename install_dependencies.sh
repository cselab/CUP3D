#!/bin/bash

# TODO: sha256 and md5 check for Mac.

set -e

# https://stackoverflow.com/a/23378780/2203044
# LOGICAL_CORE_COUNT=$([[ $(uname) = 'Darwin' ]] && sysctl -n hw.logicalcpu_max || lscpu -p | egrep -v '^#' | wc -l)
PHYSICAL_CORE_COUNT=$([[ $(uname) = 'Darwin' ]] && sysctl -n hw.physicalcpu_max || lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)

# Parameters modifiable from environment.
JOBS=${JOBS:-$PHYSICAL_CORE_COUNT}
SOURCES=${SOURCES:-$PWD/dependencies}
INSTALL_PATH=${INSTALL_PATH:-$PWD/dependencies/build}

# Shorthands for version. Note that changing this number may not be enough for the script to work properly!
CMAKE_VERSION=3.11.1
CMAKE_VERSION_SHORT=3.11
CMAKE_SHA_256='57bebc6ca4d1d42c6385249d148d9216087e0fda57a47dc5c858790a70217d0c  cmake-3.11.1.tar.gz'

FFTW_VERSION=3.3.7
FFTW_MD5='0d5915d7d39b3253c1cc05030d79ac47  fftw-3.3.7.tar.gz'

HDF5_VERSION=1.10.1
HDF5_MD5='43a2f9466702fb1db31df98ae6677f15  hdf5-1.10.1.tar.gz'
HDF5_URL='https://www.hdfgroup.org/package/source-gzip/?wpdmdl=4301&refresh=5afee8d8a45151526655192'

GSL_VERSION=2.1

# Other shorthands.
TAR="tar --keep-newer-files"

# Flags. By default all are disabled.
INSTALL_CMAKE=
INSTALL_FFTW=
INSTALL_HDF5=
INSTALL_GSL=
UNKNOWN_ARGUMENT=
PRINT_EXPORT=

# Determine what users wants.
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -a|--all) INSTALL_CMAKE=1; INSTALL_FFTW=1; INSTALL_HDF5=1; INSTALL_GSL=1; shift ;;
        -e|--export) PRINT_EXPORT=1; shift ;;
        --cmake) INSTALL_CMAKE=1; shift ;;
        --fftw) INSTALL_FFTW=1; shift ;;
        --hdf5) INSTALL_HDF5=1; shift ;;
        --gsl) INSTALL_GSL=1; shift ;;
        *) UNKNOWN_ARGUMENT=1; shift ;;
    esac
done


if [ -z "$INSTALL_CMAKE" -a -z "$INSTALL_FFTW" -a -z "$INSTALL_HDF5" -a -z "$INSTALL_GSL" -a -z "$PRINT_EXPORT" -o -n "$UNKNOWN_ARGUMENT" ]; then
    echo "Usage:
    ./install_dependencies [-a | --all | [[--cmake] [--fftw] [--hdf5] [--gsl]]]

Arguments:
  -a,  --all    - Install all available libraries and tools
  -e,  --export - Print export commands for all libraries and tools
                  (assuming they are installed)
  --cmake       - Install CMake ${CMAKE_VERSION} (required at least 3.2)
  --fftw        - Install FFTW ${FFTW_VERSION}
  --hdf5        - Install HDF5 ${HDF5_VERSION}
  --gsl         - Install GSL ${GSL_VERSION}

All libraries and tools are installed locally in the dependencies/ folder.
Note that this script tries to not redo everything from scratch if run multiple
times. Therefore, in case of errors, try erasing dependencies/ folder."
    exit
fi


BASEPWD=$PWD

if [ -n "$INSTALL_CMAKE" ]; then
    echo "Installing CMake ${CMAKE_VERSION}..."
    wget -nc https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}.tar.gz -P $SOURCES
    cd $SOURCES
    [ -x "$(command -v sha256sum)" ] && sha256sum --quiet -c - <<< $CMAKE_SHA_256
    $TAR -xzvf cmake-${CMAKE_VERSION}.tar.gz
    cd cmake-${CMAKE_VERSION}
    ./bootstrap --parallel=${JOBS} --prefix=$INSTALL_PATH/cmake-${CMAKE_VERSION}/
    make -j${JOBS}
    make install -j${JOBS}
    cd $BASEPWD
fi


if [ -n "$INSTALL_FFTW" ]; then
    echo "Installing FFTW ${FFTW_VERSION}..."
    wget -nc http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz -P $SOURCES
    cd $SOURCES
    [ -x "$(command -v md5sum)" ] && md5sum --quiet -c - <<< $FFTW_MD5
    $TAR -xzvf fftw-${FFTW_VERSION}.tar.gz
    cd fftw-${FFTW_VERSION}
    ./configure --prefix=$INSTALL_PATH/fftw-${FFTW_VERSION}/ --enable-mpi --enable-threads --enable-shared
    make -j${JOBS}
    make install -j${JOBS}
    cd $BASEPWD
fi


if [ -n "$INSTALL_HDF5" ]; then
    echo "Installing HDF5 ${HDF5_VERSION}..."
    mkdir -p $SOURCES
    wget ${HDF5_URL} -O $SOURCES/hdf5-${HDF5_VERSION}.tar.gz
    cd $SOURCES
    [ -x "$(command -v md5sum)" ] && md5sum --quiet -c <<< $HDF5_MD5
    $TAR -xzvf hdf5-${HDF5_VERSION}.tar.gz
    cd hdf5-${HDF5_VERSION}
    CC=mpicc ./configure --prefix=$INSTALL_PATH/hdf5-${HDF5_VERSION}-parallel/ --enable-parallel
    make -j $JOBS
    make install -j $JOBS
    cd $BASEPWD
fi


if [ -n "$INSTALL_GSL" ]; then
    echo "Installing GSL ${GSL_VERSION}..."
    wget -nc http://mirror.switch.ch/ftp/mirror/gnu/gsl/gsl-${GSL_VERSION}.tar.gz -P $SOURCES
    cd $SOURCES
    $TAR -xzvf gsl-${GSL_VERSION}.tar.gz
    cd gsl-${GSL_VERSION}
    ./configure --prefix=$INSTALL_PATH/gsl-${GSL_VERSION} --enable-parallel
    make -j $JOBS
    make install -j $JOBS
    cd $BASEPWD
fi


[[ -z "$PRINT_EXPORT" ]] && echo "======================================================================"
[[ -z "$PRINT_EXPORT" ]] && echo "Done! Run or add to ~/.bashrc the following command(s)"
[[ -z "$PRINT_EXPORT" ]] && echo "to capture installed dependencies in your environment:"
echo
[ -n "$PRINT_EXPORT" -o -n "$INSTALL_CMAKE" ] && echo "export PATH=$INSTALL_PATH/cmake-${CMAKE_VERSION}/bin:\$PATH"
[ -n "$PRINT_EXPORT" -o -n "$INSTALL_FFTW" ] && echo "export FFTWDIR=$INSTALL_PATH/fftw-${FFTW_VERSION}/"
[ -n "$PRINT_EXPORT" -o -n "$INSTALL_HDF5" ] && echo "export HDF5_ROOT=$INSTALL_PATH/hdf5-${HDF5_VERSION}-parallel/"
[ -n "$PRINT_EXPORT" -o -n "$INSTALL_GSL" ] && echo "export GSL_ROOT_DIR=$INSTALL_PATH/gsl-${GSL_VERSION}/"
