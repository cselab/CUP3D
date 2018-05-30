<!--
TO VIEW THIS FILE, RUN THE FOLLOWING:
    python3 -m pip install --user grip
    python3 -m grip README.md --export README.html

AND OPEN
    README.html

OR USE WEB SERVER VARIANT (NOTE: 60 UPDATES/HOUR LIMIT!!)
    python3 -m grip README.md
-->

# Compilation

## Quick compilation

To compile, run:
```bash
cd makefiles/
cmake .
make
```

If that doesn't work, read the following section.

## Detailed instructions

CubismUP uses CMake to automatically detect required dependencies and to compile the code.
If the dependencies are missing, they can be easily downloaded, compiled and locally installed using the provided script `install_dependencies.sh` (details below).
If the dependencies are already available, but CMake does not detect them, appropriate environment variables specifying their path have to be defined.

CubismUP requries the following 3rd party libraries:

| Dependency            | Environment variable pointing to the existing installation |
|-----------------------|----------------------------------|
| FFTW (3.3.7) (\*)     | $FFTWDIR                         |
| HDF5 (1.10.1) (\*)    | $HDF5_ROOT                       |
| GSL (2.1) (\*)        | $GSL_ROOT_DIR                    |
| MPI (\*\*)            | [See instruction][mpi-path]      |

(\*) Possibly higher versions work too.<br>
(\*\*) Especially if installing the dependencies, make sure that `mpicc` points to a MPI-compatible `C` compiler, and `mpic++` to a MPI-compatible `C++` compiler.

We suggest first trying to compile the code with the libraries already installed on the target machine or cluster.
If available, dependencies may be loaded with `module load ...` or `module load new ...`.
If `module load` is not available, but libraries are installed, set the above mentioned environment variables.

## Manually installing dependencies

To install the missing dependencies, run the following code (from the repository root folder):
```bash
# Step 1: Install dependencies
./install_dependencies.sh --all

# Step 2: Append the export commands to ~/.bashrc or ~/.bash_profile
./install_dependencies.sh --export >> ~/.bashrc
# or
./install_dependencies.sh --export >> ~/.bash_profile  # (Mac)

# Step 3:
source ~/.bashrc
# or
source ~/.bash_profile  # (Mac)

# Step 4: Try again
cmake .
```

## Other options and further info

The `--all` flag installs all dependencies known to the script (FFTW, HDF5, GSL, as well as CMake itself).
If only some dependencies are missing, pass instead flags like `--cmake`, `--fftw` and othes.
Run `./install_dependencies.sh` to get the full list of available flags.

All dependencies are installed in the folder `dependencies/`.
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



# Cluster-specific modules

Piz Daint:
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


[mpi-path]: https://stackoverflow.com/questions/43054602/custom-mpi-path-in-cmake-project
