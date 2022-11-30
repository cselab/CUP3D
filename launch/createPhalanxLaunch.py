import numpy as np
import argparse
import math

def launch_script(dx,dy,levels,fish_per_level,name):

	f = open(name + ".sh", "w")
	f.write(
'#!/bin/bash\n\
NNODE=16\n\
PSOLVER=\"iterative\"\n\
FACTORY=\n\
FACTORY+=\"StefanFish L=0.2 T=1.0 xpos=0.5 ypos=2.0 zpos=2.0 CorrectPosition=true CorrectPositionZ=true CorrectRoll=true heightProfile=baseline widthProfile=baseline bFixFrameOfRef=1\n\"\n')

	FishOfThisLevel = 1
	xstart = 0.5
	ystart = 2.0
	zstart = 2.0
	for level in range(1,levels+1):
		for fish in range(1,fish_per_level):
			sign = 1
			if fish % 2 == 0:
				sign = -1
			x = xstart 
			#if (fish-1) % 4 <= 1:
			#	x = xstart + 0.05
			y = ystart + sign * int((fish+1)/2)*dy
			z = zstart
			f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectPositionZ=true CorrectRoll=true  heightProfile=baseline widthProfile=baseline bFixFrameOfRef=0\n\"\n'.format(x,y,z))
		xstart += dx

	# WRITE SOLVER SETTINGS
	f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx 4.0"\n\
OPTIONS+=" -checkpointsteps 10000000 "\n\
OPTIONS+=" -bpdx 6 -bpdy 6 -bpdz 6"\n\
OPTIONS+=" -tdump 0.25 -tend 100.0"\n\
OPTIONS+=" -CFL 0.4 -nu 0.000008 -lambda 1e10"\n\
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4"\n\
OPTIONS+=" -levelMax 8 -levelStart 4 -levelMaxVorticity 7 -Rtol 1.0 -Ctol 0.1"\n\
OPTIONS+=" -poissonSolver ${PSOLVER}"')


dx=1.0
dy=0.05
levels=1
launch_script(dx,dy,levels,1,"phalanx1")
launch_script(dx,dy,levels,3,"phalanx3")
launch_script(dx,dy,levels,5,"phalanx5")
launch_script(dx,dy,levels,7,"phalanx7")
