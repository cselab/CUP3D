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

Also dependent on smarties directory (location for smarties must be your home directory):
`git@gitlab.ethz.ch:mavt-cse/smarties.git`

Note: If you are putting these `module load`s into a script,
use `source script_name`, not `./script_name`.
