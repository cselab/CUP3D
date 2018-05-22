#!/bin/bash
BASENAME=FlowPastStefan
FFACTORY=factoryStefanSmall

OPTIONS=
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 4"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.5"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.0000625"
OPTIONS+=" -tend 20 -tdump 1"
