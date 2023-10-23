#!/bin/bash
NNODE=64 #number of nodes to use for a simulation

#----------------------------------
#Settings for cylinder
#----------------------------------
XPOS=${XPOS:-1.0}      #cylinder center (x-coordinate)
XVEL=${XVEL:-0.2}      #cylinder velocity (x-component)
DIAMETER=${DIAMETER:-0.1}  #cylinder diameter
HALFLENGTH=${HALFLENGTH:-1.25}  #cylinder halflength

#----------------------------------
#Settings for pressure equation 
#----------------------------------
#PSOLVER="iterative"       #CPU solver
PSOLVER="cuda_iterative"  #GPU solver
PT=${PT:-1e-8}            #absolute error tolerance
PTR=${PTR:-0}             #relative error tolerance

#----------------------------------
#Settings for simulation domain
#----------------------------------
EXTENT=${EXTENT:-4}    #length of largest side
BPDX=${BPDX:-4}        #number of blocks in x-side, at refinement level = 0
BPDY=${BPDY:-2}        #number of blocks in y-side, at refinement level = 0
BPDZ=${BPDZ:-2}        #number of blocks in z-side, at refinement level = 0
#
# Coarsest possible mesh (at refinement level = 0) is a
# (BPDX * BS) x (BPDY * BS) x (BPDZ * BS) grid, where BS is the number 
# of grid points per block (equal to 8, by default).

#--------------------------------------
#Settings for Adaptive Mesh Refinement
#--------------------------------------
RTOL=${RTOL:-10.}              #grid is refined when curl(u) > Rtol (u: fluid velocity)
CTOL=${CTOL:-1.}               #grid is compressed when curl(u) < Ctol (u: fluid velocity)
LEVELS=${LEVELS:-8}            #maximum number of refinement levels allowed
LEVELSSTART=${LEVELSSTART:-4}  #at t=0 the grid is uniform and at this refinement level. Must be strictly less than LEVELS.

#--------------------------------------
#Other settings
#--------------------------------------
NU=${NU:-0.000010} #fluid viscosity
#The Reynolds number is defined as Re = XVEL * DIAMETER / NU and can be controlled by
#modifying the values of NU. Here are some examples:
# Re=1000  <-> NU=0.000040
# Re=2000  <-> NU=0.000020
# Re=4000  <-> NU=0.000010

#--------------------------------------
#Timestep and file saving
#--------------------------------------
TDUMP=${TDUMP:-0.1}   #Save files for t = i*TDUMP, i=0,1,...
TEND=${TEND:-100.}    #Perform simulation until t=TEND
CFL=${CFL:-0.50}      #Courant number: controls timestep size (should not exceed 1.0).
VERBOSE=${VERBOSE:-1} #Set to 1 for more verbose screen output.

OPTIONS="-bpdx $BPDX -bpdy $BPDY -bpdz $BPDZ -levelMax $LEVELS -levelStart $LEVELSSTART -Rtol $RTOL -Ctol $CTOL -extentx $EXTENT -CFL $CFL -tdump $TDUMP -nu $NU -tend $TEND -verbose $VERBOSE -poissonTol $PT -poissonTolRel $PTR -poissonSolver $PSOLVER"
OPTIONS+=" -bMeanConstraint 0 -BC_z periodic " #hardcoded periodic BCs in spanwise direction
FACTORY="CylinderNozzle xpos=$XPOS bFixFrameOfRef=1 bForcedInSimFrame=1 xvel=$XVEL L=$DIAMETER halflength=$HALFLENGTH ";
