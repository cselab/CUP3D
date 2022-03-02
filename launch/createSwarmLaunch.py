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
"#!/bin/bash\n\
NNODE=256\n\
BPDX=${BPDX:-8}\n\
BPDY=${BPDY:-4}\n\
BPDZ=${BPDZ:-4}\n\
NU=${NU:-0.000005}\n\
BC=${BC:-freespace}\n\
\n\
FACTORY=\n"
)

x0 = 1.0
y0 = 2.0
z0 = 2.0

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
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} bForcedInSimFrame=1 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=stefan\n\"\n'.format(xOuter,yOuter,zOuter))
        First = False
      else:   
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} bForcedInSimFrame=1 bFixToPlanar=1 heightProfile=danio widthProfile=stefan\n\"\n'.format(xOuter,yOuter,zOuter))

      if ( i != numFish-1 ) and ( k!= numFish-1 ):# ADD FISH ON X-Z FACE
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} bForcedInSimFrame=1 bFixToPlanar=1 heightProfile=danio widthProfile=stefan\n\"\n'.format(xInner,yOuter,zInner))
      if ( i != numFish-1 ) and ( j!= numFish-1 ):# ADD FISH ON X-Y FACE
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} bForcedInSimFrame=1 bFixToPlanar=1 heightProfile=danio widthProfile=stefan\n\"\n'.format(xInner,yInner,zOuter))
      if ( j != numFish-1 ) and ( k!= numFish-1 ):# ADD FISH ON Y-Z FACE
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} bForcedInSimFrame=1 bFixToPlanar=1 heightProfile=danio widthProfile=stefan\n\"\n'.format(xOuter,yInner,zInner))
 
# WRITE SOLVER SETTINGS
f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx 8.0"\n\
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"\n\
OPTIONS+=" -tdump 0.025 -tend 100.0 "\n\
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"\n\
OPTIONS+=" -CFL 0.6 -nu ${NU}"\n\
OPTIONS+=" -poissonTol 5e-7 -poissonTolRel 1e-4"\n\
OPTIONS+=" -levelMax 8 -levelStart 5 -Rtol 1.0 -Ctol 0.01"\n\
')
