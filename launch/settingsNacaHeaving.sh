#!/bin/bash

#PARTITION=debug
#WCLOCK=00:30:00
BASENAME=Naca_Heaving_05
NNODE=8
NNODEX=8
NNODEY=1

# this is a heaving and pitching naca, free to move in the x direction from the fluid forces
# sim box moves with the x center of mass of the airfoil
# Fheave and Fpitch are frequencies of the two motions
# Apitch and Aheave (non dimensional) are the amplitudes of the two motions
# Ppitch is the phase, divided by 2*M_PI (therefore if Ppitch=0.5 they are antiphase)
# Mpitch is the mean pitch. if Mpitch=0.1 the pitching motion is sinusoidal around 0.1 radians
FACTORY='IF3D_NacaOperator L=0.2 thickness=0.12 xpos=0.25 Fheave=1 Aheave=0.05 Fpitch=0 Apitch=0 Ppitch=0 Mpitch=0 bForcedInSimFrame_y=1 bForcedInSimFrame_z=1 bFixFrameOfRef_x=1
'

OPTIONS=
OPTIONS+=$(printf ' -factory-content %q' "$FACTORY")
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -2Ddump 0 -restart 0"
OPTIONS+=" -3Ddump 1"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -uinfx 0.0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -length 0.2"
OPTIONS+=" -lambda 1e5"
#OPTIONS+=" -nu 0.000015625"
OPTIONS+=" -nu 8e-6"
OPTIONS+=" -tend 20 -tdump 0.05"
