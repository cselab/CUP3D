#!/bin/bash

BASENAME=2CarlingsRe5k_CFL01
NNODE=8
NNODEX=8
NNODEY=1
WCLOCK=12:00
WSECS=43000
FFACTORY=factory2Carlings

OPTIONS=
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 32"
OPTIONS+=" -2Ddump 1 -restart 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.15"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0000045"
OPTIONS+=" -tend 20"
