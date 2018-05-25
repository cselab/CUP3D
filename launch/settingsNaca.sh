#!/bin/bash

BASENAME=Naca_Forced_highres_00
NNODE=32
NNODEX=32
NNODEY=1

# Naca in this case is forced. No heaving / pitching.
# airfoil moves with speed xvel and sim box moves with airfoil if bFixFrameOfRef_x
FACTORY='IF3D_NacaOperator L=0.2 thickness=0.15 xpos=0.25 xvel=0.2 bFixFrameOfRef_x=1 bForcedInSimFrame=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 16"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -3Ddump 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5 -DLM 10"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 0.000008"
OPTIONS+=" -tend 20 -tdump 0.05"
