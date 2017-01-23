#!/bin/bash

BASENAME=2Stefans
NNODE=32
NNODEX=16
NNODEY=2
#WCLOCK=48:00
#WSECS=172000
WCLOCK=12:00
WSECS=43000
FFACTORY=factory2Stefans

OPTIONS=
OPTIONS+=" -nActions 0"
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.15"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.00005"
OPTIONS+=" -tend 80 -tdump 0"
