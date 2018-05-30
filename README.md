<!--
TO VIEW THIS FILE, RUN THE FOLLOWING:
    python3 -m pip install --user grip
    python3 -m grip README.md --export README.html

AND OPEN
    README.html
-->

# Cluster-specific modules

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

# Compilation

CubismUP depends on the following 3rd party tools and libraries:

  - C++11 compiler supporting OpenMP
  - MPI (*)
  - CMake (3.2 or higher)
  - FFTW (3.3.7) (**)
  - HDF5 (1.10.1) (**)
  - GSL (2.1) (**)

(\*) If manually installing dependencies (see below), make sure `mpicc` points to a MPI-compatible `C` compiler, and that `mpic++` point to a MPI-compatible `C++` compiler.<br>
(\*\*) Possibly higher versions work too.

We suggest first trying to compile the code with the libraries already installed on the target machine or cluster
(if available, dependencies may be loaded with `module load ...` or `module load new ...`).
Otherwise, an installation script is provided for all dependencies (except MPI and C++ compiler, for which we assume are already available).

To compile, run:
```bash
cmake .
make
```

If `cmake .` fails, use the provided script for installing the missing dependencies:
```bash
# Step 1: Install dependencies
./install_dependencies.sh --all

# Step 2: Append the export commands to ~/.bashrc or ~/.bash_profile (Mac):
./install_dependencies.sh --export >> ~/.bashrc
# or
./install_dependencies.sh --export >> ~/.bash_profile

# Step 3:
source ~/.bashrc
# or
source ~/.bash_profile

# Step 4: Try again
cmake .
```


## Other options and further info

The `--all` flag installs all available dependencies. To install only some of them, run `./install_dependencies.sh` to see the full list of flags.

All dependencies are installed in the folder `./dependencies/`.
Full installation takes 5-15 minutes, depending on the machine.
To specify number of parallel jobs in the internal `make`, write
```
JOBS=10 ./install_dependencies.sh [flags]
```
The default number of jobs is equal to the number of physical cores.


## Compiling on Mac

The default compiler `clang` on Mac does not support OpenMP. It is therefore necessary either to install `clang` OpenMP extension or to install e.g. `g++` compiler. The following snippet shows how to compile with `g++-7`:
```bash
CXX=g++-7 cmake .
MPICH_CXX=g++-7 make
# or
OMPI_CXX=g++-7 make
```


## Troubleshooting

If `cmake .` keeps failing, delete the file `CMakeCache.txt` and try again.
