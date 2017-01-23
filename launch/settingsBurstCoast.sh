#!/bin/bash

#BASENAME=CarlingFishRe4000_BurstCoast_efficiency_fixed
BASENAME=CarlingFishRe4000_BurstCoast_velocity_fixed

NNODE=16
NNODEX=8
NNODEY=2
#WCLOCK=48:00
#WSECS=172500
WCLOCK=12:00
#WSECS=43000
FFACTORY=factoryCarlingBurstCoast
mkdir -p ${BASEPATH}${BASENAME}
#cp burst_coast_carling_params_re4000_eff.txt ${BASEPATH}${BASENAME}/burst_coast_carling_params.txt
cp burst_coast_carling_params_re4000_vel.txt ${BASEPATH}${BASENAME}/burst_coast_carling_params.txt

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 0 -restart 0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -tend 20"
