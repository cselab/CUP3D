#!/bin/bash
NNODE=${NNODE:-16}
PSOLVER="iterative"

# Re = 5000 = L^2/T/NU => NU=0.000008, for L=0.2 and T=1.0

FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.5 ypos=0.50 bFixToPlanar=0 bFixFrameOfRef=1 CorrectPosition=1 CorrectZ=1 CorrectRoll=1
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.8 ypos=0.55 bFixToPlanar=0 bFixFrameOfRef=0 CorrectPosition=1 CorrectZ=1 CorrectRoll=1
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.8 ypos=0.45 bFixToPlanar=0 bFixFrameOfRef=0 CorrectPosition=1 CorrectZ=1 CorrectRoll=1
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -extent 2.0 -bpdx 4 -bpdy 2 -bpdz 2"
OPTIONS+=" -tdump 1.0 -tend 100.0"
OPTIONS+=" -CFL 0.8 -nu 0.000008 -lambda 1e8 "
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -levelMax 6 -levelStart 3 -Rtol 2.0 -Ctol 0.1"
OPTIONS+=" -restart 0 -checkpointsteps 1000 "
