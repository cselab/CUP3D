#!/bin/bash
BASENAME=Dcyl_forced_07
NNODE=4
NNODEX=4
NNODEY=1
WCLOCK=12:00
FFACTORY=factoryDcyl

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 8"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.1"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00001"
OPTIONS+=" -tend 20 -tdump 0.1"


