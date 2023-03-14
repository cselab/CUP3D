#!/bin/bash
NNODE=${NNODE:-16}
LEVELS=${LEVELS:-4}
# PSOLVER="iterative" #CPU Poisson solver
PSOLVER="cuda_iterative" #GPU Poisson solver

FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.5 bFixFrameOfRef=1 CorrectPosition=0 CorrectPositionZ=0 CorrectRoll=0 heightProfile=larval widthProfile=larval
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER} -implicitDiffusion 1"
OPTIONS+=" -extent 2.0 -bpdx 4 -bpdy 2 -bpdz 2"
OPTIONS+=" -tdump 0.1 -tend 10.0" 
OPTIONS+=" -CFL 0.3 -nu 0.001 -lambda 1e12 " # Re=40
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -levelMax 9 -levelMaxVorticity ${LEVELS} -levelStart 3 -Rtol 1.0 -Ctol 0.1"
# OPTIONS+=" -restart 0 -checkpointsteps 1000 "

