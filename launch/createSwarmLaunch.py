import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('--numFish', help='number of fish on the edges of the cube', required=True, type=int)
parser.add_argument('--x0', help='first position in x', required=True, type=float)
parser.add_argument('--y0', help='first position in y', required=True, type=float)
parser.add_argument('--z0', help='first position in z', required=True, type=float)
parser.add_argument('--delta', help='x/y spacing between the fish', required=True, type=float)
args = vars(parser.parse_args())
numFish = args['numFish']
x0 = args['x0']
y0 = args['y0']
z0 = args['z0']
delta = args['delta']

f = open("settingsCarlingSwarm.sh", "w")
## WRITE GRID CONFIGURATION
## FOR AMR CODE
# f.write(
# "#!/bin/bash\n\
# NNODE=8\n\
# \n\
# BPDX=${BPDX:-64}\n\
# BPDY=${BPDY:-64}\n\
# BPDZ=${BPDZ:-64}\n\
# \n\
# NU=${NU:-0.00004}\n\
# BC=${BC:-freespace}\n\
# \n\
# L=0.2\n\
# \n\
# FACTORY=\n"
# )

## FOR UNIFORM CODE
f.write(
"#!/bin/bash\n\
NNODEX=${NNODEX:-2}\n\
NNODEY=${NNODEY:-2}\n\
NNODEZ=${NNODEZ:-2}\n\
NNODE=8\n\
\n\
BPDX=${BPDX:-64}\n\
BPDY=${BPDY:-64}\n\
BPDZ=${BPDZ:-64}\n\
\n\
NU=${NU:-0.00001}\n\
BC=${BC:-freespace}\n\
\n\
L=0.2\n\
\n\
FACTORY=\n"
)

counter = 0
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
      if counter == 0:
        f.write('FACTORY+="CarlingFish L=0.1 T=1.0 xpos={} ypos={} zpos={} bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan\n\
"\n'.format(xOuter,yOuter,zOuter))
      else:   
        f.write('FACTORY+="CarlingFish L=0.1 T=1.0 xpos={} ypos={} zpos={} bFixToPlanar=1 heightProfile=stefan widthProfile=stefan\n\
"\n'.format(xOuter,yOuter,zOuter))
      counter+=1
      # ADD FISH ON X-Z FACE
      if ( i != numFish-1 ) and ( k!= numFish-1 ):
        f.write('FACTORY+="CarlingFish L=0.1 T=1.0 xpos={} ypos={} zpos={} bFixToPlanar=1 heightProfile=stefan widthProfile=stefan\n\
"\n'.format(xInner,yOuter,zInner))
        counter+=1
      # ADD FISH ON X-Y FACE
      if ( i != numFish-1 ) and ( j!= numFish-1 ):
        f.write('FACTORY+="CarlingFish L=0.1 T=1.0 xpos={} ypos={} zpos={} bFixToPlanar=1 heightProfile=stefan widthProfile=stefan\n\
"\n'.format(xInner,yInner,zOuter))
        counter+=1
      # ADD FISH ON Y-Z FACE
      if ( j != numFish-1 ) and ( k!= numFish-1 ):
        f.write('FACTORY+="CarlingFish L=0.1 T=1.0 xpos={} ypos={} zpos={} bFixToPlanar=1 heightProfile=stefan widthProfile=stefan\n\
"\n'.format(xOuter,yInner,zInner))
        counter+=1
 
# WRITE SOLVER SETTINGS
f.write('\nOPTIONS=\n\
OPTIONS+="-extentx 4.0"\n\
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"\n\
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 50.0 "\n\
OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "\n\
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"\n\
OPTIONS+=" -CFL 0.3 -use-dlm 10 -nu ${NU}"\n\
OPTIONS+=" -levelMax 7 -levelStart 5 -Rtol 1.0 -Ctol 0.1"\n\
OPTIONS+=" -Advection3rdOrder=true"')

# PRIN THE COUNT OF FISH PLACED
print("wrote settingsCarlingSwarm.sh for {} fish".format(counter))
