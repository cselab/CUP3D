#!/bin/bash

#PARTITION=debug
#WCLOCK=00:30:00
BASENAME=Naca_Heaving_highres_02
NNODE=32
NNODEX=32
NNODEY=1

# this is a heaving and pitching naca, free to move in the x direction from the fluid forces
# sim box moves with the x center of mass of the airfoil
# Fheave and Fpitch are frequencies of the two motions
# Apitch and Aheave (non dimensional) are the amplitudes of the two motions
# Ppitch is the phase, divided by 2*M_PI (therefore if Ppitch=0.5 they are antiphase)
# Mpitch is the mean pitch. if Mpitch=0.1 the pitching motion is sinusoidal around 0.1 radians
FACTORY='IF3D_NacaOperator L=0.2 thickness=0.12 xpos=0.25 xvel=0.15155 Fheave=0.16165 Aheave=0.75 Fpitch=0.16165 Apitch=0.5235987 Ppitch=0.208333 Mpitch=0 bForcedInSimFrame=1 bFixFrameOfRef_x=1 bFixFrameOfRef_y=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 128 -bpdy 64 -bpdz 16"
OPTIONS+=" -dump2D 1 -restart 0"
OPTIONS+=" -dump3D 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5 -DLM 10"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 0.00002755428299"
OPTIONS+=" -tend 8 -tdump 0.008"
