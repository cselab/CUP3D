#!/bin/bash

#BASENAME=2Stefans
#NNODE=8
#NNODEX=8
#NNODEY=1
#WCLOCK=48:00
#WSECS=172000
#WCLOCK=12:00
#WSECS=43000
FFACTORY=factory2Stefans

OPTIONS=
OPTIONS+=" -nActions 0"
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 8"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.000008"
OPTIONS+=" -tend 0"
OPTIONS+=" -tdump 0.0"
OPTIONS+=" -fdump 0"
OPTIONS+=" -nsteps 101"
