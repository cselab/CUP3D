#!/bin/bash
BASENAME=Stefan_Dcyl_07
NNODE=32
NNODEX=32
NNODEY=1
WCLOCK=12:00
FFACTORY=factoryDcyl_Fish

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 16"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.3"
OPTIONS+=" -uinfx 0.3"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00001"
OPTIONS+=" -tend 20 -tdump 0.5"


