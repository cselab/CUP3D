#!/bin/bash

BASENAME=Naca_Vortex_09
NNODE=32
NNODEX=16
NNODEY=2
FFACTORY=factoryNaca
#FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -uinfx 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 8e-6"
OPTIONS+=" -tend 20 -tdump 0.05"
