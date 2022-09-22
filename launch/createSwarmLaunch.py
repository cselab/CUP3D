import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('--levels', help='number of levels'      , required=True, type=int  )
parser.add_argument('--dx'    , help='x spacing between fish', required=True, type=float)
parser.add_argument('--dy'    , help='y spacing between fish', required=True, type=float)
parser.add_argument('--dz'    , help='z spacing between fish', required=True, type=float)
parser.add_argument('--name'  , help='name of settings file' , required=True, type=str  )
args = vars(parser.parse_args())
dx   = args['dx']
dy   = args['dy']
dz   = args['dz']
name = args['name']
levels = args['levels']

f = open(name + ".sh", "w")

f.write(
'#!/bin/bash\n\
NNODE=128\n\
PSOLVER=\"cuda_iterative\"\n\
FACTORY=\n\
FACTORY+=\"StefanFish L=0.2 T=1.0 xpos=0.5 ypos=2.0 zpos=2.0 CorrectPosition=true CorrectPositionZ=true CorrectRoll=true heightProfile=baseline widthProfile=baseline bFixFrameOfRef=1\n\"\n')

FishOfThisLevel = 1

xstart = 0.5
zstart = 2.0

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
OPTIONS+=" -extentx 8.0"\n\
OPTIONS+=" -bpdx 8 -bpdy 4 -bpdz 4"\n\
OPTIONS+=" -tdump 0.1 -tend 100.0"\n\
OPTIONS+=" -CFL 0.8 -nu 0.000008"\n\
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4"\n\
OPTIONS+=" -levelMax 9 -levelStart 4 -Rtol 5.0 -Ctol 0.1"\n\
OPTIONS+=" -poissonSolver ${PSOLVER}"')
