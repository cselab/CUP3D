#!/bin/bash
NNODE=64
DLM=${DLM:-1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-32}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
LEVELS=${LEVELS:-5}
CFL=${CFL:-0.2} # if 0, DT is used
DT=${DT:-1e-4}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}
BC=${BC:-dirichlet}

# L=0.2 stefanfish Re=1'000 <-> NU=0.00004
NU=${NU:-0.00004}

FACTORY='StefanFish L=0.2 T=1 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 5"
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL ${CFL} -dt ${DT} -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -ImplicitPenalization 1"
OPTIONS+=" -levelMax ${LEVELS} -levelStart 1 -Rtol 4.00 -Ctol 1.00"
OPTIONS+=" -extentx 4.0 " # -extentx 2.0
OPTIONS+=" -TimeOrder 2 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -dumpOmegaX 1 -dumpOmegaY 1 dumpOmegaZ 1"
OPTIONS+=" -dumpVelocityX 1 -dumpVelocityY 1 -dumpVelocityZ 1"