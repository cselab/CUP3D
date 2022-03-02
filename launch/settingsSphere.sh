#!/bin/bash
NNODE=${NNODE:-64}
DLM=${DLM:-1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-32} #16
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
LEVELS=${LEVELS:-7}
CFL=${CFL:-0.1} # if 0, DT is used
DT=${DT:-0}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}

# Re=300 <-> NU=0.00005208333333; Re=500 <-> NU=0.00003125; Re=1000 <-> NU=0.000015625; Re=10'000 <-> NU=0.0000015625; Re=420'000 <-> NU=0.000000037202381; Re=1'140'000 <-> NU=0.00000001370614
NU=${NU:-0.000015625}

FACTORY='Sphere L=0.125 xpos=0.6 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1 bBreakSymmetry=1'

OPTIONS=
OPTIONS+=" -tdump 0 -tend 10" # -fdump 1
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL ${CFL} -dt ${DT} -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -levelMax ${LEVELS} -levelStart 1 -Rtol 4.00 -Ctol 1.00"
OPTIONS+=" -extentx 4.0 " # -extentx 2.0
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -dumpOmegaX 1 -dumpOmegaY 1 dumpOmegaZ 1"
OPTIONS+=" -dumpVelocityX 1 -dumpVelocityY 1 -dumpVelocityZ 1"
