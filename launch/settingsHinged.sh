#!/bin/bash

BASENAME=prod_1_effic_equal_try2
NNODE=4
NNODEX=2
NNODEY=2
WCLOCK=24:00:00
WSECS=43000
FFACTORY=factoryCarlingHinged

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 64 -bpdz 16"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e6"
OPTIONS+=" -nu 0.00008"
OPTIONS+=" -tend 2 -tdump 0.005"
