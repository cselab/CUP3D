#!/bin/bash
NNODE=1

BPDX=${BPDX:-2}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}
LEVELS=${LEVELS:-5}
NU=${NU:-0.0625}

FACTORY='Pipe L=2.5 xpos=1.3 ypos=1.3 halflength=2.6 zvel=0.0 bForcedInSimFrame=1 bFixFrameOfRef=1 
'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -extentx 2.6"
OPTIONS+=" -tdump 0.1 -tend 100.00 "
OPTIONS+=" -CFL 0.3 -nu ${NU} -implicitDiffusion 1"
OPTIONS+=" -levelMax ${LEVELS} -levelStart 4 -Rtol 0.4 -Ctol 0.1"
OPTIONS+=" -lambda 1e12"
OPTIONS+=" -dumpP 1"
OPTIONS+=" -dumpChi 1"
OPTIONS+=" -dumpOmega 1"
OPTIONS+=" -dumpOmegaX 1"
OPTIONS+=" -dumpOmegaY 1"
OPTIONS+=" -dumpOmegaZ 1"
OPTIONS+=" -dumpVelocity 1"
OPTIONS+=" -dumpVelocityX 1"
OPTIONS+=" -dumpVelocityY 1"
OPTIONS+=" -dumpVelocityZ 1"
OPTIONS+=" -BC_x freespace"
OPTIONS+=" -BC_y freespace"
OPTIONS+=" -BC_z freespace"
