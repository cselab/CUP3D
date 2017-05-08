#!/bin/bash
BASENAME=Stefans_PID_10_newestCorrectionMethod_amp100_ang10_10
NNODE=4
NNODEX=4
NNODEY=1
FFACTORY=factoryStefan_PID

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 8"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy 1"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.000008"
OPTIONS+=" -tend 20 -tdump 0.05"
