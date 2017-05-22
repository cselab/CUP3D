#!/bin/bash
BASENAME=Stefans_PID_highRes_$FACTORYNAME
NNODE=64
NNODEX=32
NNODEY=2
FFACTORY=factoryStefan_PID_$FACTORYNAME
echo $BASENAME $FFACTORY
OPTIONS=
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 16"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.000008"
OPTIONS+=" -tend 20 -tdump 0.05"
