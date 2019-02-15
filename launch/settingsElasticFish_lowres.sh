#!/bin/bash
NNODE=4
NNODEX=4
NNODEY=1
#WCLOCK=12:00

FACTORY='IF3D_ElasticFishOperator L=0.5 T=1 xvel=0.2 xpos=0.65 bFixToPlanar=1 bForcedInSimFrame=1
'

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 8"
OPTIONS+=" -dump2D 1  -dump3D 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.5"
OPTIONS+=" -use-dlm 10"
OPTIONS+=" -nu 0.00003"
OPTIONS+=" -tend 10 -tdump 0.05"
