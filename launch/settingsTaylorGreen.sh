#!/bin/bash
NNODEX=1
NNODEY=1
NNODE=$(($NNODEX * $NNODEY))

FACTORY=''

OPTIONS=
OPTIONS+=" -bpdx 16 -bpdy 16 -bpdz 16"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1"
OPTIONS+=" -BC_x periodic -BC_y periodic -BC_z periodic -initCond taylorGreen"
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1 -tend 10 -uMax_forced 1 -compute-dissipation 1"
# RE = 8/12 Ly Umax / \nu = 66.666
OPTIONS+=" -nu 0.0002"
