#!/bin/bash
NNODE=64

BPDX=${BPDX:-8}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
NU=${NU:-0.0002} #Re=U*L/nu = 1*0.1/nu = 500

# This follows the paper "Iterative Brinkman penalization for simulation of impulsively started flow past a sphere and a circular disc" (2017)
# L is the diameter.
# We set the thickness equal to L/16 (halflength = thickness/2)
FACTORY='Cylinder L=0.1 xpos=1.0 ypos=1.0 zpos=1.0 xvel=0.0 yvel=0.0 zvel=1.0 bFixFrameOfRef=1 bForcedInSimFrame=1 halflength=0.003125
'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -extentx 2.0"
OPTIONS+=" -tdump 0.1 -fdump 1 -tend 0.01 "
OPTIONS+=" -CFL 0.1 -nu ${NU}"
OPTIONS+=" -levelMax 7 -levelStart 4 -Rtol 0.4 -Ctol 0.1"
OPTIONS+=" -lambda 1e12"
OPTIONS+=" -dumpP 1"
OPTIONS+=" -dumpChi 1"
OPTIONS+=" -dumpOmega 1"
OPTIONS+=" -dumpOmegaX 1"
OPTIONS+=" -dumpOmegaY 1"
OPTIONS+=" -dumpOmegaZ 1"
OPTIONS+=" -dumpVelocity 1"
OPTIONS+=" -dumpVelocityX 1"
OPTIONS+=" -dumpVelocityY 1"
OPTIONS+=" -dumpVelocityZ 1"
