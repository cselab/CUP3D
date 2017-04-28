Daint compilation:

```shell
module load gcc
module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel
module load fftw
module load daint-gpu
module load cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1
module load GSL/2.1-CrayGNU-2016.11
```

Also dependent on smarties directory (location for smarties must be your home directory):
`git@gitlab.ethz.ch:mavt-cse/smarties.git`

Note: If you are putting these `module load`s into a script,
use `source script_name`, not `./script_name`.
