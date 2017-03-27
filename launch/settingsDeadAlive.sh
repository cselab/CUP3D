#!/bin/bash

BASENAME=DeadAlive_00
NNODE=8
NNODEX=8
NNODEY=1
FFACTORY=factoryDeadAlive

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 8"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1 -lambda 1e5"
OPTIONS+=" -nu 0.000008 -length 0.2"
OPTIONS+=" -tend 20 -tdump 0.05"
