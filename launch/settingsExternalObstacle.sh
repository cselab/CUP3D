#!/bin/bash
NNODE=${NNODE:-64}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-16}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}
LEVELS=${LEVELS:-8}
CFL=${CFL:-0.2} # if 0, DT is used
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}
NU=${NU:-0.00000625} #Re=10,000
BC=${BC:-dirichlet}
FACTORY='ExternalObstacle L=0.5 xpos=0.5 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1 externalObstaclePath=/users/mchatzim/AMR_codes/3D_CUP/launch/externalObstacles/Robot.ply
'
OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 1000.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL ${CFL} -lambda ${LAMBDA} -nu ${NU}"
OPTIONS+=" -ImplicitPenalization 1"
OPTIONS+=" -levelMax ${LEVELS} -Rtol 0.50 -Ctol 0.05 -levelStart 3"
OPTIONS+=" -extentx 8.0 "
OPTIONS+=" -TimeOrder 2 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -dumpOmegaX 0 -dumpOmegaY 0 dumpOmegaZ 0"
OPTIONS+=" -dumpVelocityX 1 -dumpVelocityY 1 -dumpVelocityZ 1"
