#!/bin/bash

#PARTITION=debug
#WCLOCK=00:30:00
BASENAME=Naca_Heaving_05
NNODE=8
NNODEX=8
NNODEY=1
FFACTORY=factoryNacaHeaving

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -3Ddump 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -uinfx 0.0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 8e-6"
OPTIONS+=" -tend 20 -tdump 0.05"
