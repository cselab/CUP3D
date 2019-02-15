#!/bin/bash
NNODE=4
NNODEX=4
NNODEY=1

FACTORY='IF3D_StefanFish L=0.3 T=1 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1
'

OPTIONS=
OPTIONS+=" -nActions 0"
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 16"
OPTIONS+=" -dump2D 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.3"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.0000225"
OPTIONS+=" -tend 8 -tdump 0.05"
