#!/bin/bash

BASENAME=FlowPastCarlingFishRe0550_CFL01_testSave_03
NNODE=4
NNODEX=4
NNODEY=1
#WCLOCK=48:00
#WSECS=172000
WCLOCK=24:00
WSECS=43000
FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 1 -restart 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1 -tdump 0.5"
OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0001136363636"
OPTIONS+=" -tend 6"
