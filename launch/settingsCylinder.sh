#!/bin/bash
NNODE=64
PSOLVER="cuda_iterative"
# L is the diameter.
FACTORY='Cylinder L=0.2 xpos=0.5 xvel=0.2 yvel=0.0 zvel=0.0 bFixFrameOfRef=1 bForcedInSimFrame=1 halflength=0.5'
OPTIONS=
OPTIONS+=" -bpdx 8 -bpdy 4 -bpdz 4"
OPTIONS+=" -extentx 2.0"
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -tdump 0.1 -tend 200 "
OPTIONS+=" -CFL 0.5 -nu 0.00001 "
OPTIONS+=" -levelMax 6 -levelStart 4 -Rtol 5.0 -Ctol 0.5"
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
