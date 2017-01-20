#!/bin/bash

CFL=0.1
WCLOCK=12:00
WSECS=43200
#BASENAME=FlowPastFixedCylRe050_CFL${CFL}_testScaling
#NNODE=1
FFACTORY=factoryFixedSphere

OPTIONS=
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 16"
OPTIONS+=" -2Ddump 1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL ${CFL}"
OPTIONS+=" -uinfx 0.1"
OPTIONS+=" -length 0.1"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0002"
OPTIONS+=" -tend 10"
OPTIONS+=" -tdump 1"
OPTIONS+=" -nsteps 1"
