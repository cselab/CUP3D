#!/bin/bash

BASENAME=CarlingFishRe400_BurstCoast_finned
NNODE=4
#WCLOCK=48:00
#WSECS=172500
WCLOCK=12:00
#WSECS=43000
FFACTORY=factoryCarlingBurstCoast
mkdir -p ${BASEPATH}${BASENAME}
cp burst_coast_carling_params_re400_eff.txt ${BASEPATH}${BASENAME}/burst_coast_carling_params.txt

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 1 -restart 0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.2"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00015625"
