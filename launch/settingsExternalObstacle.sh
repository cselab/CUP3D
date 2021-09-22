#!/bin/bash
NNODE=${NNODE:-64}
DLM=${DLM:--1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-2} #16
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}
LEVELS=${LEVELS:-7}
CFL=${CFL:-0.2} # if 0, DT is used
DT=${DT:-1e-4}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}
# Re=100 <-> 0.00125
NU=${NU:-0.000625}
BC=${BC:-dirichlet}
FACTORY='ExternalObstacle L=0.5 xpos=0.5 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1 externalObstaclePath=/users/pweber/korali/examples/study.cases/CUP3D/_deps/CUP-3D/launch/externalObstacles/Robot.ply
'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0 -tend 0 -nsteps 1" # -fdump 1
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL ${CFL} -dt ${DT} -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -ImplicitPenalization 1"
OPTIONS+=" -levelMax ${LEVELS} -Rtol 4.00 -Ctol 1.00" #-levelStart 3 
OPTIONS+=" -extentx 1.0 " # -extentx 2.0
OPTIONS+=" -TimeOrder 2 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -dumpOmegaX 1 -dumpOmegaY 1 dumpOmegaZ 1"
OPTIONS+=" -dumpVelocityX 1 -dumpVelocityY 1 -dumpVelocityZ 1"
