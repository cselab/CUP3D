#!/bin/bash

BASENAME=CarlingFishRe4000_BurstCoast
NNODE=2
FFACTORY=factoryCarling
cp burst_coast_carling_params.txt ${FOLDER}/burst_coast_carling_params.txt

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 1"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.000015625"
