#!/bin/bash

BASENAME=FlowPastFixedCyl
FFACTORY=factoryFixedSphere

OPTIONS=
OPTIONS+=" -bpdx 16 -bpdy 8 -bpdz 8"
OPTIONS+=" -2Ddump 0"
OPTIONS+=" -tdump 10"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -nsteps 101"
OPTIONS+=" -uinfx 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0004"
OPTIONS+=" -tend 8"
