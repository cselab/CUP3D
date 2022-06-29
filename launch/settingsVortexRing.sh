#!/bin/bash
NNODE=${NNODE:-32}
PT=${PT:-1e-8}
PTR=${PTR:-1e-4}

PSOLVER="cuda_iterative"
COND="vorticity"
BC="periodic"
FACTORY= #no solid bodies

OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -bpdx 2 -bpdy 2 -bpdz 16"
OPTIONS+=" -initCond ${COND}"
OPTIONS+=" -tdump 0. -fdump 50 -tend 0.0 -nsteps 10000 -extent 16.0 -bMeanConstraint 0 "
OPTIONS+=" -CFL 0.4 -nu ${NU} -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -levelMax 3 -levelStart 2 -Rtol 1000.00 -Ctol 100.0"
#OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
#OPTIONS+=" -dumpOmegaX 1 "
#OPTIONS+=" -dumpOmegaY 1 "
#OPTIONS+=" -dumpOmegaZ 1 "
