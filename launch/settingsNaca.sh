#!/bin/bash

PARTITION=debug
WCLOCK=00:30:00
BASENAME=Naca_Forced_17
NNODE=4
NNODEX=4
NNODEY=1
FFACTORY=factoryNaca
#FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 8"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -3Ddump 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -uinfx 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 8e-6"
OPTIONS+=" -tend 20 -tdump 0.05"
