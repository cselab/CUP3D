#!/bin/bash
NNODE=32
DLM=${DLM:-1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-4}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}
LEVELS=${LEVELS:-6}
CFL=${CFL:-0.4} # if 0, DT is used
DT=${DT:-1e-4}
PT=${PT:-1e-5}
PTR=${PTR:-1e-2}
BC=${BC:-dirichlet}

# L=0.2 stefanfish Re=1'000 <-> NU=0.00004
NU=${NU:-0.00004}

# bottom fish
FACTORY="StefanFish L=0.2 T=1 xpos=1.20 ypos=1.00 zpos=0.80 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan
"

# plane 1 with 4 fish
FACTORY+="StefanFish L=0.2 T=1 xpos=0.90 ypos=1.00 zpos=0.90 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=0.90 zpos=0.90 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=1.10 zpos=0.90 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.50 ypos=1.00 zpos=0.90 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"

# plane 2 with 9 fish
FACTORY+="StefanFish L=0.2 T=1 xpos=0.60 ypos=1.00 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=0.90 ypos=0.90 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=0.90 ypos=1.10 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=0.80 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=1.00 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=1.20 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.50 ypos=0.90 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.50 ypos=1.10 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.80 ypos=1.00 zpos=1.00 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"

# plane 3 with 4 fish
FACTORY+="StefanFish L=0.2 T=1 xpos=0.90 ypos=1.00 zpos=1.10 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=0.90 zpos=1.10 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=1.10 zpos=1.10 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1 xpos=1.50 ypos=1.00 zpos=1.10 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan
"

# top fish
FACTORY+="StefanFish L=0.2 T=1 xpos=1.20 ypos=1.00 zpos=1.20 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1 heightProfile=danio widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 5"
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL ${CFL} -dt ${DT} -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -ImplicitPenalization 1"
OPTIONS+=" -levelMax ${LEVELS} -levelStart 4 -Rtol 4.00 -Ctol 1.00"
OPTIONS+=" -extentx 4.0 " # -extentx 2.0
OPTIONS+=" -TimeOrder 2 "
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -dumpOmegaX 1 -dumpOmegaY 1 dumpOmegaZ 1"
OPTIONS+=" -dumpVelocityX 1 -dumpVelocityY 1 -dumpVelocityZ 1"
