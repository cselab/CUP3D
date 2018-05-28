#!/bin/bash
BASENAME=Stefan_Dcyl_08
NNODE=4
NNODEX=4
NNODEY=1
#WCLOCK=12:00

FACTORY='IF3D_DCylinder L=0.1 xpos=0.2 xvel=0.3 computeForces=0 bFixFrameOfRef=1 bForcedInSimFrame=1
IF3D_StefanFish L=0.3 T=1 xpos=0.45 bFixToPlanar=1 xvel=0.3 bForcedInSimFrame=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 8"
OPTIONS+=" -dump2D 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.3"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00001"
OPTIONS+=" -tend 20 -tdump 0.05"


