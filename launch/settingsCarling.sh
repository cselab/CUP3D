#!/bin/bash

BASENAME=FlowPastCarlingFishRe0550_CFL01_finned_fixed
NNODE=8
NNODEX=8
NNODEY=1
#WCLOCK=48:00
#WSECS=172000
WCLOCK=48:00
WSECS=43000
FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0001136363636"
OPTIONS+=" -tend 8"
