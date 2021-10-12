import numpy as np
import argparse
import math


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--nodes', required=True, type=int)
  args = vars(parser.parse_args())
  nodes = args['nodes']

  f = open("settingsWeakScaling-nodes" + str(nodes) + ".sh", "w")
  f.write(\
"#!/bin/bash\n\
NNODE="+str(nodes*nodes*nodes)+"\n\
BPDX=${BPDX:-"+str(4*nodes)+"}\n\
BPDY=${BPDY:-"+str(4*nodes)+"}\n\
BPDZ=${BPDZ:-"+str(4*nodes)+"}\n\
NU=${NU:-0.00001}\n\
BC=${BC:-freespace}\n\
\n\
\n\
FACTORY=\n")
  for k in range(nodes):
    for j in range(nodes):
      for i in range(nodes):
        f.write('FACTORY+="StefanFish L=0.2 T=1.0 xpos={} ypos={} zpos={} bForcedInSimFrame=1 bFixFrameOfRef=1 Profile=stefan widthProfile=stefan\n\"\n'.format(0.5+i,0.5+j,0.5+k))

  # WRITE SOLVER SETTINGS
  f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx '+str(nodes)+'"\n\
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"\n\
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0 -tend 1.0 "\n\
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"\n\
OPTIONS+=" -CFL 0.7 -use-dlm -1 -nu ${NU}"\n\
OPTIONS+=" -levelMax 5 -levelStart 3 -Rtol 1.0 -Ctol 0.01"\n\
OPTIONS+=" -TimeOrder 2"\n')
