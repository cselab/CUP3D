#!/bin/bash
NNODE=${NNODE:-64}
DLM=${DLM:-0}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-4}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}
LEVELS=${LEVELS:-8}
CFL=${CFL:-0.7}
PT=${PT:-1e-6}
PTR=${PTR:-1e-4}

# L=0.2 stefanfish Re=1'000 <-> NU=0.00004
NU=${NU:-0.00001} #Re=4000

# plane 1 with 4 fish
FACTORY+="StefanFish xvel=0.2  L=0.2 T=1 xpos=0.90 ypos=1.00 bForcedInSimFrame=1 bFixFrameOfRef=1  heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish xvel=0.2  L=0.2 T=1 xpos=1.20 ypos=0.90 bForcedInSimFrame=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish xvel=0.2  L=0.2 T=1 xpos=1.20 ypos=1.10 bForcedInSimFrame=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish xvel=0.2  L=0.2 T=1 xpos=1.50 ypos=1.00 bForcedInSimFrame=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0.1 -tend 50"
OPTIONS+=" -CFL ${CFL} -dt ${DT} -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -levelMax ${LEVELS} -levelStart 4 -Rtol 4.00 -Ctol 1.00"
OPTIONS+=" -extentx 4.0 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
