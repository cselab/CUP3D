#!/bin/bash
NNODE=${NNODE:-64}
BPDX=${BPDX:-16}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
LEVELS=${LEVELS:-5}
CFL=${CFL:-0.2} # if 0, DT is used
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}

# L=0.2 stefanfish Re=1'000 <-> NU=0.00004
NU=${NU:-0.00004}

FACTORY='StefanFish L=0.2 T=1 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan'

OPTIONS=
OPTIONS+=" -tdump 0.1 -tend 5"
OPTIONS+=" -CFL ${CFL} -nu ${NU}"
OPTIONS+=" -levelMax ${LEVELS} -levelStart 1 -Rtol 4.00 -Ctol 1.00"
OPTIONS+=" -extentx 2.0 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
