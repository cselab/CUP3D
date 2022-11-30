#!/bin/bash
NNODE=${NNODE:-8}
PSOLVER="iterative"

CFL=$3
RTOL=$4
CTOL=$5
LEVELMAX=$6
LEVELMAXVORTICITY=$((LEVELMAX-1))
BX=$7
BY=$((BX/2))
BZ=$((BX/2))

FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.5 bFixFrameOfRef=1 CorrectRoll=0 CorrectPosition=0 CorrectPositionZ=0
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -extent 2.0 -bpdx ${BX} -bpdy ${BY} -bpdz ${BZ}"
OPTIONS+=" -tdump 0 -tend 20.0 -freqDiagnostics 10 "
OPTIONS+=" -restart 0 -checkpointsteps 1000000 "
OPTIONS+=" -nu 0.000008 -lambda 1e10 -levelStart 3"
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -CFL ${CFL} "
OPTIONS+=" -levelMax ${LEVELMAX} "
OPTIONS+=" -levelMaxVorticity ${LEVELMAXVORTICITY} "
OPTIONS+=" -Rtol ${RTOL} -Ctol ${CTOL}"
