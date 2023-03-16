#!/bin/bash

NNODE=4
BPDX=${BPDX:-16}
BPDY=${BPDY:-${BPDX}}
BPDZ=${BPDZ:-$((${BPDX}/8))}

NU=${NU:-0.1} # RE = 8/12 Ly Umax / \nu = 66.666

FACTORY=''

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0.1 -tend 100.00 "
OPTIONS+=" -CFL 0.3 -nu ${NU} -implicitDiffusion 1"
OPTIONS+=" -BC_x periodic -BC_y wall -BC_z periodic -initCond channel -fixedMassFlux 0 -uMax_forced 1"
OPTIONS+=" -nu ${NU}"
OPTIONS+=" -dumpP 1"
OPTIONS+=" -dumpOmega 1"
OPTIONS+=" -dumpOmegaX 1"
OPTIONS+=" -dumpOmegaY 1"
OPTIONS+=" -dumpOmegaZ 1"
OPTIONS+=" -dumpVelocity 1"
OPTIONS+=" -dumpVelocityX 1"
OPTIONS+=" -dumpVelocityY 1"
OPTIONS+=" -dumpVelocityZ 1"
