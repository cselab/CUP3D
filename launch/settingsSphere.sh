#!/bin/bash

WCLOCK=12:00
WSECS=43200
BASENAME=FlowPastFixedCyl_FFTW_000
NNODE=4
FFACTORY=factoryFixedSphere

OPTIONS=
#OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 1 -3Ddump 0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -uinfx 0.2"
OPTIONS+=" -length 0.1"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00001"
OPTIONS+=" -tend 10"
