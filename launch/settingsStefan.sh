#!/bin/bash
BASENAME=Stefan_CFL01_BPD64
NNODE=4
NNODEX=4
NNODEY=1

FACTORY='IF3D_StefanFish L=0.3 T=1 xpos=0.25 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -nActions 0"
OPTIONS+=" -bpdx 64 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.3"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.000018"
OPTIONS+=" -tend 8 -tdump 0"
