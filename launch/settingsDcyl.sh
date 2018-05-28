#!/bin/bash
BASENAME=Dcyl_forced_13
NNODE=4
NNODEX=4
NNODEY=1

# L is the diameter. 
# Obstacle is both fixed and forced. We impose an uniform velocity and the sim box follows it
FACTORY='IF3D_DCylinder L=0.1 xpos=0.2 xvel=0.1 bFixFrameOfRef=1 bForcedInSimFrame=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 8"
OPTIONS+=" -dump2D 1 -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.1"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00001"
OPTIONS+=" -tend 20 -tdump 0.1"


