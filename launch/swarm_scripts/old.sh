#!/bin/bash
NNODE=32

BPDX=${BPDX:-8}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}

NU=${NU:-0.00004}
BC=${BC:-freespace}

L=0.2

FACTORY=
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=0.5 zpos=0.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=0.5 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=0.75 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=0.75 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=0.5 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=0.5 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=0.75 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=0.75 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=0.5 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=0.75 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.0 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.0 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.25 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.25 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.0 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.0 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.25 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.25 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.0 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.25 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.5 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.5 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.5 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5 ypos=1.5 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.25 ypos=1.5 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=0.5 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.5 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.75 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=0.75 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=0.5 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.5 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.75 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=0.75 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=0.5 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.75 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.0 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.0 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.25 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.25 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.0 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.0 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.25 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.25 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.0 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.25 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.5 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.5 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.5 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.5 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.75 ypos=1.5 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=0.5 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=0.75 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=0.5 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=0.75 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=0.5 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.0 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.25 zpos=0.75 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.0 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.25 zpos=1.25 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.0 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.5 zpos=0.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.5 zpos=1.0 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.25 ypos=1.5 zpos=1.5 bFixToPlanar=1 heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 2.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 100.0 "
OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.3 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 5 -levelStart 3 -Rtol 0.1 -Ctol 0.01"
OPTIONS+=" -Advection3rdOrder=true"
