#!/bin/bash
NNODE=${NNODE:-4}
PSOLVER="iterative" #CPU Poisson solver
#PSOLVER="cuda_iterative" #GPU Poisson solver

#Simulation of a single fish
#
# CarlingFish : prescribed midline motion
# StefanFish  : prescribed midline curvature, midline motion computed from Frenet equations
#
# Re = 1000 = L^2/T/NU => NU=0.000040, for L=0.2 and T=1.0
# Re = 5000 = L^2/T/NU => NU=0.000008, for L=0.2 and T=1.0
#
# Other options:
# CorrectPosition  = {0,1} : control yaw angle and x,y-position (available for StefanFish only)
# CorrectPositionZ = {0,1} : control pitch angle and z-position (available for StefanFish only)
# CorrectRoll      = {0,1} : control roll angle by imposing angular velocity (available for StefanFish only)
#
# heightProfile = baseline,danio,stefan
# widthProfile = baseline,danio,stefan
#

FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.5 bFixFrameOfRef=1 CorrectPosition=0 CorrectPositionZ=0 CorrectRoll=0 heightProfile=baseline widthProfile=baseline
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -extent 2.0 -bpdx 4 -bpdy 2 -bpdz 2"
OPTIONS+=" -tdump 0.1 -tend 100.0"
OPTIONS+=" -CFL 0.3 -nu 0.000008 -lambda 1e10 "
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -levelMax 7 -levelStart 3 -levelMaxVorticity 6 -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -restart 0 -checkpointsteps 1000 "

