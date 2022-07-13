#!/bin/bash
NNODE=${NNODE:-32}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}

PSOLVER="iterative"
#PSOLVER="cuda_iterative"

# L=0.2 stefanfish Re=1'000 <-> NU=0.00004
NU=${NU:-0.00004}

FACTORY='StefanFish L=0.2 T=1 xpos=0.3 bFixToPlanar=0 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan'

OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -bpdx 16 -bpdy 8 -bpdz 8"
OPTIONS+=" -tdump 0.0 -tend 10.0"
OPTIONS+=" -CFL 0.4 -nu ${NU} -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -levelMax 5 -levelStart 1 -Rtol 4.00 -Ctol 1.00"
OPTIONS+=" -restart 0 -checkpointsteps 250 "
