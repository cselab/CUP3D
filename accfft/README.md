## AccFFT

AccFFT is a new massively parallel FFT library for CPU/GPU architectures.
It is specifically designed with the goal of achieving maximum performance.

## Installation

Please see http://accfft.org/articles/install/ for installation guide.

*ON DAINT*
Modules to load:
```
module load daint-gpu
module swap PrgEnv-cray PrgEnv-gnu
module load cray-python/3.6.1.1 cray-hdf5-parallel cray-fftw
module load cudatoolkit/9.2.148_3.19-6.0.7.1_2.1__g3d9acc8 CrayGNU/.18.08
module load GSL/2.5-CrayGNU-18.08 CMake/3.12.0
```
In this folder: `mkdir build && cd build`
Ask for a node: `salloc -N 1 -C gpu -p debug`
Install:
```
export FFTWDIR=$FFTW_DIR/../
srun cmake -DBUILD_GPU=true -DBUILD_SHARED=false  ..
srun make -j
srun make install
```

## Documentation
A detailed step by step documentation on how to use the library
can be found at http://accfft.org/documentation/.

Original research paper introducing AccFFT library:
http://arxiv.org/abs/1506.07933

## License
AccFFT is distributed under GNU GENERAL PUBLIC LICENSE Version 2.
Please see LICENSE file.

