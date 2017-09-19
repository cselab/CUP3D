#!/bin/bash

BASENAME=ForcedSphereStartup_15_04
NNODE=4
FFACTORY=factoryForcedSphereStartup

OPTIONS=
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 32"
#OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 16"
OPTIONS+=" -2Ddump 0"
OPTIONS+=" -tdump 0.0"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -uinfx 0.0"
OPTIONS+=" -length 0.1"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 5e-5"
OPTIONS+=" -nsteps 101"
OPTIONS+=" -tend 8"
