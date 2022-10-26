import numpy as np
import argparse
import math

def launch_script(dx,dy,dz,name,levels):

	f = open(name + ".sh", "w")
	f.write(
'#!/bin/bash\n\
NNODE=16\n\
PSOLVER=\"iterative\"\n\
FACTORY=\n\
FACTORY+=\"StefanFish L=0.2 T=1.0 xpos=0.5 ypos=2.0 zpos=2.0 CorrectPosition=true CorrectPositionZ=true CorrectRoll=true heightProfile=baseline widthProfile=baseline bFixFrameOfRef=1\n\"\n')

	FishOfThisLevel = 1
	xstart = 0.5
	zstart = 2.0
	kount = -1
	for level in range(1,levels):
	  FishOfThisLevel += 1
	  xstart += dx
	  ystart = 2.0 - level*dy
	  for fish in range(FishOfThisLevel):
	    x = xstart
	    y = ystart + fish * 2* dy
	    z = zstart
	    f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectPositionZ=true CorrectRoll=true heightProfile=baseline widthProfile=baseline bFixFrameOfRef=0\n\"\n'.format(x,y,z))
 
	# WRITE SOLVER SETTINGS
	f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx 4.0"\n\
OPTIONS+=" -bpdx 4 -bpdy 4 -bpdz 4"\n\
OPTIONS+=" -tdump 0.2 -tend 50.0"\n\
OPTIONS+=" -CFL 0.5 -nu 0.000008 -lambda 1e12"\n\
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4"\n\
OPTIONS+=" -levelMax 9 -levelStart 4 -levelMaxVorticity 8 -Rtol 5.0 -Ctol 0.2"\n\
OPTIONS+=" -poissonSolver ${PSOLVER}"')

'''
parser = argparse.ArgumentParser()
parser.add_argument('--levels', help='number of levels'      , default=2,     type=int  )
parser.add_argument('--dx'    , help='x spacing between fish', required=True, type=float)
parser.add_argument('--dy'    , help='y spacing between fish', required=True, type=float)
parser.add_argument('--dz'    , help='z spacing between fish', default=0.   , type=float)
parser.add_argument('--name'  , help='name of settings file' , required=True, type=str  )
args = vars(parser.parse_args())
dx   = args['dx']
dy   = args['dy']
dz   = args['dz']
name = args['name']
levels = args['levels']
launch_script(dx,dy,dz,name,levels)
'''

'''
DeltaX = np.linspace(0.30,0.50,8)
DeltaY = np.linspace(0.05,0.20,8)
dz     = 0.0
levels = 4
kount = 1
dx_array = []
dy_array = []
dx_array.append(0.3285714285714285)
dy_array.append(0.05)
dx_array.append(0.3857142857142857)
dy_array.append(0.05)
dx_array.append(0)
dy_array.append(0)
for dx in DeltaX:
	for dy in DeltaY:
		name = 'test'+str(kount)
		dx_array[2] = dx
		dy_array[2] = dy
		print(kount,dx,dy)
		launch_script(dx_array,dy_array,dz,name,levels)
		kount += 1
'''
dx=0.45
dy=0.075
dz=0.0
launch_script(dx,dy,dz,"levels2",2)
launch_script(dx,dy,dz,"levels3",3)
launch_script(dx,dy,dz,"levels4",4)
launch_script(dx,dy,dz,"levels5",5)
