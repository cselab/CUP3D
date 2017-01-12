#!/bin/bash

BASENAME=FlowPastFixedCylRe050_scaling08node${_MY_CC_}
NNODE=8
FFACTORY=factoryFixedSphere

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
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
