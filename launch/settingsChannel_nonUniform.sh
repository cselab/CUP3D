#!/bin/bash
NNODEX=${NNODEX:-4}
NNODEY=1
NNODE=$(($NNODEX * $NNODEY))

BPDX=${BPDX:-8}
BPDY=${BPDY:-${BPDX}} #${BPDY:-32}
BPDZ=${BPDZ:-$((${BPDX}/2))} #${BPDZ:-32}

NU=${NU:-0.0004} # RE = halfHeight * U / nu = 2500

FACTORY=''

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -extentx 6.2831853072 -extenty 2 -extentz 4.7123889804"
OPTIONS+=" -mesh_density_y SinusoidalDensity -eta_y 0.9" # first mesh point at 0.012523878 (Kim and Moin 1987)

# stuff to test with iterative penalization (almost uniform grid):
#OPTIONS+=" -extentx 2 -extenty 2 -extentz 1"
#OPTIONS+=" -mesh_density_y SinusoidalDensity -eta_y 0.5" # first mesh point at y+ = 0.018313569

OPTIONS+=" -useSolver petsc"
OPTIONS+=" -iterativePenalization 0"
OPTIONS+=" -dump2D 1 -dump3D 1 -tdump 1"
OPTIONS+=" -BC_x periodic -BC_y wall -BC_z periodic -fadeLen 0 -initCond channelRandom"
OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1 -tend 10 -uMax_forced 1 -compute-dissipation 1"
OPTIONS+=" -nu ${NU}"
