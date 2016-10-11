#!/bin/bash

CFL=0.5
WCLOCK=12:00
WSECS=43200
BASENAME=FlowPastFixedCylRe100_CFL${CFL}
NNODE=2
FFACTORY=factoryFixedSphere

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -2Ddump 1"
OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL ${CFL}"
OPTIONS+=" -uinfx 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0002"
OPTIONS+=" -tend 6"
