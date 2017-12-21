#!/bin/bash
BASENAME=Stefans_PID_lowRes_xP29_yP00
NNODE=32
NNODEX=32
NNODEY=1
FFACTORY=factoryStefan_PID_xP290_yP000
WCLOCK=06:00:00

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 8"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e6"
OPTIONS+=" -nu 0.000008"
OPTIONS+=" -tend 20 -tdump 0.05"
