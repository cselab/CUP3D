#!/bin/bash

CFL=0.1
WCLOCK=00:30:00
WSECS=43200
BASENAME=plate_Re100
NNODE=2

# removed ypos and zpos. cose should put object in the middle of the domain
# removed uinf from options, moved to factory and forced/fixed obstacle
FACTORY='IF3D_PlateObstacle L=0.2 xpos=0.25 xvel=0.1 a=0.2 b=0.4 thickness=0.05 alpha=30 bForcedInSimFrame=1 bFixFrameOfRef=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 1"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL ${CFL}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e4"
OPTIONS+=" -nu 0.0002"
OPTIONS+=" -tdump 0.01"
OPTIONS+=" -tend 0.1"
