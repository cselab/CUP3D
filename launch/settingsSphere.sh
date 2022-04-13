#!/bin/bash
NNODE=${NNODE:-32}
LEVELS=${LEVELS:-7}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}
NU=${NU:-0.00005208333333} 
# Re=300   <-> NU=0.00005208333333
# Re=500   <-> NU=0.00003125
# Re=1000  <-> NU=0.000015625
# Re=10000 <-> NU=0.0000015625

FACTORY='Sphere L=0.125 xpos=0.6 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1 bBreakSymmetry=1'

OPTIONS=
OPTIONS+=" -nu ${NU} -levelMax ${LEVELS} -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -tdump 0 -tend 10"
OPTIONS+=" -bpdx 32 -bpdy 8 -bpdz 8 -extentx 4.0"
OPTIONS+=" -CFL 0.4 -lambda 1e8 -levelStart 1 -Rtol 5.00 -Ctol 0.10"
OPTIONS+=" -dumpOmegaX 0 -dumpOmegaY 0 dumpOmegaZ 0"
OPTIONS+=" -dumpVelocityX 0 -dumpVelocityY 0 -dumpVelocityZ 0"
