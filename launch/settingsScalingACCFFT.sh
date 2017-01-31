#!/bin/bash

CFL=0.1
WCLOCK=12:00
WSECS=43200
#BASENAME=ScalingCyl_Nx64_Ny32_Nz8
#NNODE=1
FFACTORY=factoryFixedSphere

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 16"
OPTIONS+=" -2Ddump 0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL ${CFL}"
OPTIONS+=" -uinfx 0.1"
OPTIONS+=" -length 0.1"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.0002"
OPTIONS+=" -tend 10"
OPTIONS+=" -tdump 0"
OPTIONS+=" -nsteps 10001"
