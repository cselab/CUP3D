#!/bin/bash

BASENAME=Naca_Forced_18
NNODE=4
NNODEX=4
NNODEY=1

# Naca in this case is forced. No heaving / pitching.
# airfoil moves with speed xvel and sim box moves with airfoil if bFixFrameOfRef_x
FACTORY='IF3D_NacaOperator L=0.2 thickness=0.12 xpos=0.25 xvel=0.1 bFixFrameOfRef_x=1 bFixFrameOfRef=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 8"
OPTIONS+=" -2Ddump 1 -restart 0"
OPTIONS+=" -3Ddump 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 8e-6"
OPTIONS+=" -tend 20 -tdump 0.05"
