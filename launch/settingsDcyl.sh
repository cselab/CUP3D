#!/bin/bash
NNODE=32

BPDX=${BPDX:-8}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
NU=${NU:-0.0002} #Re=U*L/nu = 1*0.1/nu = 500

# L is the diameter.
# We set the thickness equal to L/16 (halflength = thickness/2)
FACTORY='Cylinder L=0.1 xpos=2.0 ypos=2.0 zpos=2.0 xvel=0.0 yvel=0.0 zvel=1.0 bFixFrameOfRef=1 bForcedInSimFrame=1 halflength=0.003125
'

OPTIONS=
OPTIONS+=" -extentx 4.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 20.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.3 -use-dlm -10 -nu ${NU}"
OPTIONS+=" -levelMax 6 -levelStart 3 -Rtol 0.5 -Ctol 0.05"
OPTIONS+=" -implicitPenalization 1"
OPTIONS+=" -TimeOrder 2"
