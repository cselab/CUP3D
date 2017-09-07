#!/bin/bash

BASENAME=carling_parahinge_00
NNODE=8
NNODEX=4
NNODEY=2
#WCLOCK=48:00
#WSECS=172000
WCLOCK=24:00:00
WSECS=43000
FFACTORY=factoryCarling
#FFACTORY=factoryStefan

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 32 -bpdz 8"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1" # -tdump 0.0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e6"
OPTIONS+=" -nu 0.00008"
#OPTIONS+=" -nu 0.000015625"
#OPTIONS+=" -nu 0.00015625"
OPTIONS+=" -tend 1 -tdump 0.01"
