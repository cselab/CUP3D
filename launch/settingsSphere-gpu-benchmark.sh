#!/bin/bash
NNODE=${NNODE:-10}
LEVELS=${LEVELS:-6}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}
TEND=${TEND:-1.0}


PSOLVER="iterative"
#PSOLVER="cuda_iterative"

FACTORY='Sphere L=0.125 xpos=0.6 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1 bBreakSymmetry=1'
OPTIONS=
OPTIONS+=" -levelMax ${LEVELS} -tend ${TEND} -poissonTol ${PT} -poissonTolRel ${PTR} -poissonSolver ${PSOLVER}"
OPTIONS+=" -nu 0.00003125" 
OPTIONS+=" -bpdx 8 -bpdy 2 -bpdz 2 -extentx 4.0"
OPTIONS+=" -CFL 0.5 -lambda 1e8 -levelStart 2 -Rtol 5.00 -Ctol 0.10"
OPTIONS+=" -dumpOmegaX 1 -dumpOmegaY 1 dumpOmegaZ 1 -tdump 0.25 "
