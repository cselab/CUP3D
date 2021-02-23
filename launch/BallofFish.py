import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky # computes upper triangle by default, matches paper
import argparse
import math

def sample(S, z_hat, m_FA):
    '''
    Samples points uniformly in ellispoid in n-dimensions.
    x^2/a^2+y^2/b^2+z^2/c^2 = 1 
    S.shape = (n,n)
    S = diag(a^2,b^2,c^2)
    m_Fa = number of points
    '''
    Gamma_Threshold = 1.0
    nz = S.shape[0]
    z_hat = z_hat.reshape(nz,1)
    X_Cnz = np.random.normal(size=(nz, m_FA))
    rss_array = np.sqrt(np.sum(np.square(X_Cnz),axis=0))
    kron_prod = np.kron( np.ones((nz,1)), rss_array)
    X_Cnz = X_Cnz / kron_prod       # Points uniformly distributed on hypersphere surface
    R = np.ones((nz,1))*( np.power( np.random.rand(1,m_FA), (1./nz)))
    unif_sph=R*X_Cnz;               # m_FA points within the hypersphere
    T = np.asmatrix(cholesky(S))    # Cholesky factorization of S => S=T’T
    unif_ell = T.H*unif_sph ; # Hypersphere to hyperellipsoid mapping
    # Translation and scaling about the center
    z_fa=(unif_ell * np.sqrt(Gamma_Threshold)+(z_hat * np.ones((1,m_FA))))
    return np.array(z_fa)

def FishSamples(a,b,c,fish,L):
  S = np.eye(3)
  S[0][0] = a**2
  S[1][1] = b**2
  S[2][2] = c**2
  z_hat = np.zeros(3)
  xyz = sample(S,z_hat,12345*fish)
  xvalid=[]
  yvalid=[]
  zvalid=[]
  xL = 1.10
  yL = 0.25
  zL = 0.25
  for i in range(xyz.shape[1]):
    xtest = xyz[0,i]
    ytest = xyz[1,i]
    ztest = xyz[2,i]
    valid = True
    for j in range(len(xvalid)):
       r = np.sqrt( ((xtest-xvalid[j])/xL)**2 + ((ytest-yvalid[j])/yL)**2 + ((ztest-zvalid[j])/zL)**2 )
       if r < 1.0*L:
          valid=False
          break
    if valid == True:
       xvalid.append(xtest)
       yvalid.append(ytest)
       zvalid.append(ztest)
       print("Valid = ", len(xvalid))
    if (len(xvalid)==fish):
       break

  return xvalid,yvalid,zvalid


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--fish', required=True, type=int)
  args = vars(parser.parse_args())
  fish = args['fish']

  L = 0.2
  x,y,z = FishSamples(2.0,1.0,1.0,fish,0.2)
  x = 3.0 + np.asarray(x) + 0.5*L
  y = 2.0 + np.asarray(y)
  z = 2.0 + np.asarray(z)

  f = open("settingsEllipsoidSwarm.sh", "w")
  f.write(\
"#!/bin/bash\n\
NNODE=32\n\
BPDX=${BPDX:-8}\n\
BPDY=${BPDY:-4}\n\
BPDZ=${BPDZ:-4}\n\
NU=${NU:-0.00001}\n\
BC=${BC:-freespace}\n\
\n\
\n\
FACTORY=\n")
  for j in range(fish):
    if j==0:
      f.write('FACTORY+="StefanFish L='+str(L)+' T=1.0 xpos={} ypos={} zpos={} bCorrectPosition=true heightProfile=danio widthProfile=danio bFixFrameOfRef=1 \n\"\n'.format(x[j],y[j],z[j]))
    else:
      f.write('FACTORY+="StefanFish L='+str(L)+' T=1.0 xpos={} ypos={} zpos={} bCorrectPosition=true heightProfile=danio widthProfile=danio\n\"\n'.format(x[j],y[j],z[j]))

  # WRITE SOLVER SETTINGS
  f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx 8.0"\n\
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"\n\
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 100.0 "\n\
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"\n\
OPTIONS+=" -CFL 0.5 -use-dlm -1 -nu ${NU}"\n\
OPTIONS+=" -levelMax 7 -levelStart 4 -Rtol 0.1 -Ctol 0.01"\n\
OPTIONS+=" -TimeOrder 2"\n')
