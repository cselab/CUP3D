#!/bin/bash

CFL=0.1
WCLOCK=00:30:00
WSECS=43200
BASENAME=plate_Re100
NNODE=2
FFACTORY=factoryPlate

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL ${CFL}"
OPTIONS+=" -uinfx 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0002"
OPTIONS+=" -tdump 0.01"
OPTIONS+=" -tend 0.1"
