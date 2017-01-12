#!/bin/bash

BASENAME=CarlingFishRe4000_BurstCoast
NNODE=8
#WCLOCK=48:00
#WSECS=172500
WCLOCK=12:00
WSECS=43000
FFACTORY=factoryCarling
mkdir -p ${BASEPATH}${BASENAME}
cp burst_coast_carling_params.txt ${BASEPATH}${BASENAME}/burst_coast_carling_params.txt

OPTIONS=
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 64"
OPTIONS+=" -2Ddump 1 -restart 1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.000015625"
