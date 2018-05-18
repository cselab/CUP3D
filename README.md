# Compilation:

Daint:
```shell
module load gcc
module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel
module load fftw
module load daint-gpu
module load cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1
module load GSL/2.1-CrayGNU-2016.11
```

Euler:
```shell
module load new modules gcc/6.3.0 open_mpi/2.1.1 fftw/3.3.4 binutils/2.25 gsl/1.16 hwloc/1.11.0 fftw_sp/3.3.4
```

## Manual installation of dependencies

If `module load` fails, try installing the dependencies manually by running `./install_dependencies.sh`. The provided script installs CMake, FFTW and HDF5, by default in the folder `dependencies/`. After it's finished, copy paste the printed `export` commands into `.bashrc` (or make sure to have those environment variables set before compiling CubismUP).

Installation takes 5-10 mins on Falcon/Panda. To specify number of parallel jobs in `make`, write
```
JOBS=10 ./install_dependencies.sh
```
The default number of jobs is equal to the number of physical cores.

Script works so far only on Linux.

Requirements:

- Cubism requires C++11-compatible compiler,
- Make sure `mpicc` points to a MPI-compatible `C` compiler, and that `mpic++` points to a MPI-compatible `C++` compiler.
