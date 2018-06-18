#!/bin/bash
BASENAME=carling_Re750
NNODEX=16
NNODEY=1
NNODE=$(($NNODEX * $NNODEY))

FACTORY='IF3D_CarlingFish L=0.2 T=1.0 xpos=0.35 ypos=0.125 zpos=0.125 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=danio
'

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 16 -bpdz 16"
OPTIONS+=" -restart 0"
OPTIONS+=" -nprocsx ${NNODEX}"
OPTIONS+=" -nprocsy ${NNODEY}"
OPTIONS+=" -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -length 0.2"
OPTIONS+=" -nu 0.00005333333333"
OPTIONS+=" -tend 15 -tdump 0.05"
