import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('--numFish', help='number of fish on the edges of the cube', required=True, type=int)
parser.add_argument('--delta', help='x/y spacing between the fish', required=True, type=float)
args = vars(parser.parse_args())
numFish = args['numFish']
delta = args['delta']

f = open("settingsSchool.sh", "w")

f.write(
'#!/bin/bash\n\
NNODE=256\n\
\n\
PSOLVER=\"cuda_iterative\"\n\
FACTORY=\n'
)

x0 = 1.0
y0 = 2.0
z0 = 2.0

a = 0
First = True
for i in np.arange(numFish):
  for j in np.arange(numFish):
    for k in np.arange(numFish):
      xOuter = x0 + i*delta
      yOuter = y0 + j*delta
      zOuter = z0 + k*delta
      xInner = x0 + delta/2 + i*delta
      yInner = y0 + delta/2 + j*delta
      zInner = z0 + delta/2 + k*delta

      # ADD FISH ON EDGE
      if First:
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan bFixFrameOfRef=1\n\"\n'.format(xOuter,yOuter,zOuter))
        First = False
        a += 1
      else:   
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan\n\"\n'.format(xOuter,yOuter,zOuter))
        a += 1

      if ( i != numFish-1 ) and ( k!= numFish-1 ):# ADD FISH ON X-Z FACE
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan\n\"\n'.format(xInner,yOuter,zInner))
        a += 1
      if ( i != numFish-1 ) and ( j!= numFish-1 ):# ADD FISH ON X-Y FACE
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan\n\"\n'.format(xInner,yInner,zOuter))
        a += 1
      if ( j != numFish-1 ) and ( k!= numFish-1 ):# ADD FISH ON Y-Z FACE
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan\n\"\n'.format(xOuter,yInner,zInner))
        a += 1
 
# WRITE SOLVER SETTINGS
f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx 8.0"\n\
OPTIONS+=" -bpdx 8 -bpdy 4 -bpdz 4"\n\
OPTIONS+=" -tdump 0.5 -tend 100.0 "\n\
OPTIONS+=" -CFL 0.7 -nu 0.00001 "\n\
OPTIONS+=" -poissonTol 1e-5 -poissonTolRel 1e-4 "\n\
OPTIONS+=" -levelMax 7 -levelStart 4 -Rtol 5.0 -Ctol 0.1"\n\
OPTIONS+=" -poissonSolver ${PSOLVER}"')

print (a)
