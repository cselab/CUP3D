#!/bin/bash
NNODE=${NNODE:-64}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-16}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}
LEVELS=${LEVELS:-8}
CFL=${CFL:-0.7} # if 0, DT is used
PT=${PT:-1e-7}
PTR=${PTR:-1e-4}
NU=${NU:-0.00000625} #Re=10,000
FACTORY='ExternalObstacle L=0.5 xpos=7.5 xvel=-0.125 bForcedInSimFrame=1 bFixFrameOfRef=1 externalObstaclePath=/users/mchatzim/AMR_codes/3D_CUP/launch/externalObstacles/Robot.ply
'
OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0.1 -tend 1000.0 "
OPTIONS+=" -CFL ${CFL} -lambda ${LAMBDA} -nu ${NU}"
OPTIONS+=" -levelMax ${LEVELS} -Rtol 4.00 -Ctol 0.05 -levelStart 3"
OPTIONS+=" -extentx 8.0 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -dumpOmegaX 0 -dumpOmegaY 0 dumpOmegaZ 0"
OPTIONS+=" -dumpVelocityX 1 -dumpVelocityY 1 -dumpVelocityZ 1"
