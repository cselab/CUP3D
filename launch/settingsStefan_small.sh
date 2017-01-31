#!/bin/bash
LAMBDA=1e5
BASENAME=StefanSmall_BPD32_trickier_04
NNODE=2
NNODEX=2
NNODEY=1
#WCLOCK=48:00
#WSECS=172000
WCLOCK=12:00
WSECS=43000
FFACTORY=factoryStefanSmall

OPTIONS=
OPTIONS+=" -nActions 0"
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 4"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.5"
OPTIONS+=" -lambda ${LAMBDA}"
OPTIONS+=" -nu 0.0000625"
OPTIONS+=" -tend 20 -tdump 1"
