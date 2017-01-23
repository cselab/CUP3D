#!/bin/bash
LAMBDA=1e5
BASENAME=Stefan_CFL01_lambda${LAMBDA}_BPD64
NNODE=4
NNODEX=4
NNODEY=1
#WCLOCK=48:00
#WSECS=172000
WCLOCK=12:00
WSECS=43000
FFACTORY=factoryStefan

OPTIONS=
OPTIONS+=" -nActions 0"
OPTIONS+=" -bpdx 64 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.3"
OPTIONS+=" -lambda ${LAMBDA}"
OPTIONS+=" -nu 0.000018"
OPTIONS+=" -tend 8 -tdump 0"
