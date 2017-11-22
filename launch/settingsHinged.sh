#!/bin/bash

BASENAME=hinged_Nu1by20
NNODE=64
NNODEX=64
NNODEY=1
WCLOCK=24:00:00
WSECS=43000
FFACTORY=factoryCarlingHinged

OPTIONS=
#OPTIONS+=" -bpdx 144 -bpdy 108 -bpdz 36"
#OPTIONS+=" -bpdx 160 -bpdy 64 -bpdz 40"
#OPTIONS+=" -bpdx 128 -bpdy 128 -bpdz 32"
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 32"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e6"
OPTIONS+=" -nu 0.000004"
OPTIONS+=" -tend 8 -tdump 0.05"
