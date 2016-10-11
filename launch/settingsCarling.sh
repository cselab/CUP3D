#!/bin/bash

BASENAME=FlowPastCarlingFishRe0550_CFL001_test2ndOrderDiv
NNODE=8
NNODEX=8
NNODEY=1
WCLOCK=12:00
WSECS=43200
FFACTORY=factoryCarling

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
<<<<<<< HEAD
OPTIONS+=" -2Ddump 1 -restart 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.01"
=======
OPTIONS+=" -2Ddump 1"
OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
>>>>>>> 0a755e62f8df1b459aa9fc998e27259cdd50918d
OPTIONS+=" -length 0.25"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0001136363636"
OPTIONS+=" -tend 8"
