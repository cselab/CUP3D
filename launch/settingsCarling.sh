#!/bin/bash

BASENAME=FlowPastCarlingFishRe0400_CFL01_lambda1e5
NNODE=8
NNODEX=8
NNODEY=1
#WCLOCK=48:00
#WSECS=172000
WCLOCK=24:00
WSECS=43000
FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1 -tdump 0.0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e5"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 0.00015625"
OPTIONS+=" -tend 8 -tdump 0"
