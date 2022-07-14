#!/bin/bash
NNODE=${NNODE:-32}
PT=${PT:-1e-5}
PTR=${PTR:-1e-3}

PSOLVER="iterative"
#PSOLVER="cuda_iterative"

# L=0.2 stefanfish Re=1'000 <-> NU=0.00004
NU=${NU:-0.00001}

FACTORY='StefanFish L=0.2 T=1 xpos=0.3 zpos=0.55 bFixToPlanar=0 bFixFrameOfRef=1 Correct=1 bCorrectZ=1 heightProfile=danio widthProfile=stefan'

OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -bpdx 4 -bpdy 4 -bpdz 4"
OPTIONS+=" -tdump 0.0 -tend 100.0"
OPTIONS+=" -CFL 0.4 -nu ${NU} -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -levelMax 5 -levelStart 4 -Rtol 40.00 -Ctol 1.00"
OPTIONS+=" -restart 0 -checkpointsteps 250 "
