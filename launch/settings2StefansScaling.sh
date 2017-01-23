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
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 16"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.15"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.00005"
OPTIONS+=" -tend 0.051"
OPTIONS+=" -tdump 0.04"
OPTIONS+=" -fdump 0"
OPTIONS+=" -nsteps 101"
