#!/bin/bash

BASENAME=DeadFishRe0550_CFL001_test2ndOrderDiv
NNODE=8
NNODEX=8
NNODEY=1
FFACTORY=factoryDeadFish
WCLOCK=12:00
# WSECS=43000

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 1 restart 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -uinfx 0.3"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.00008"
OPTIONS+=" -tend 10"
