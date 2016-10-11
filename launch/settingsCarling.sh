#!/bin/bash

BASENAME=FlowPastCarlingFishRe0550_CFL001_test2ndOrderDiv
NNODE=8
NNODEX=8
NNODEY=1
FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 1 -restart 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.01"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0001136363636"
OPTIONS+=" -tend 8"
