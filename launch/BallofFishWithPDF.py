import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
import argparse
import math
import pandas as pd
from scipy import interpolate 
import random

def sample(S, samples):
    '''
    Samples points uniformly in ellispoid in n-dimensions.
    x^2/a^2+y^2/b^2+z^2/c^2 = 1 
    S.shape = (n,n)
    S = diag(a^2,b^2,c^2)
    samples = number of points
    '''
    Gamma_Threshold = 1.0
    nz = S.shape[0]
    z_hat = np.zeros((nz,1))#z_hat.reshape(nz,1)
    X_Cnz = np.random.normal(size=(nz, samples))
    rss_array = np.sqrt(np.sum(np.square(X_Cnz),axis=0))
    kron_prod = np.kron( np.ones((nz,1)), rss_array)
    X_Cnz = X_Cnz / kron_prod       # Points uniformly distributed on hypersphere surface
    R = np.ones((nz,1))*( np.power( np.random.rand(1,samples), (1./nz)))
    unif_sph=R*X_Cnz;               # "samples" points within the hypersphere
    T = np.asmatrix(cholesky(S))    # Cholesky factorization of S => S=T’T
    unif_ell = T.H*unif_sph ; # Hypersphere to hyperellipsoid mapping
    # Translation and scaling about the center
    z_fa=(unif_ell * np.sqrt(Gamma_Threshold)+(z_hat * np.ones((1,samples))))
    return np.array(z_fa)


def FishSamples(a,b,c,L):
  S = np.eye(3)
  S[0][0] = a**2
  S[1][1] = b**2
  S[2][2] = c**2
  xyz = sample(S,100000)

  #read histogram of nearest neighbor distance
  distribution = pd.read_csv("FishPDF.csv",delimiter=',')
  fishCM = 30.0 #centimetres
  dist_x = (1.0/fishCM)*np.asarray(distribution["distance"])
  dist_y =              np.asarray(distribution["value"]   )
  dist_y = 1.0/(np.trapz(dist_y,dist_x)) *  dist_y #normalize to create PDF
  maxf = np.max(dist_y)

  dist_x_new = np.zeros(len(dist_x)+1)
  dist_x_new[1:] = dist_x
  dist_y_new = np.zeros(len(dist_y)+1)
  dist_y_new[1:] = dist_y
  dist_x_new = np.append(dist_x_new,100.0)  
  dist_y_new = np.append(dist_y_new,  0.0)  

  f = interpolate.interp1d(dist_x_new, dist_y_new, kind='linear')
  V = 4.*np.pi/3.*a*b*c

  xvalid=[]
  yvalid=[]
  zvalid=[]
  dnn=[]
  for i in range(xyz.shape[1]):
    xtest = xyz[0,i]
    ytest = xyz[1,i]
    ztest = xyz[2,i]
    valid = True

    for j in range(len(xvalid)):
        r = (xtest-xvalid[j])**2 +(ytest-yvalid[j])**2+(ztest-zvalid[j])**2
        r = r / L / L
        if r < dnn[j]*dnn[j]:
           valid = False
           break

    if valid == True:
       xvalid.append(xtest)
       yvalid.append(ytest)
       zvalid.append(ztest)
       for ijk in range(100000):
          uu  = random.uniform(0, 1.0)
          rnn = random.uniform(0, 3.0)
          if uu < f(rnn)/maxf:
             dnn.append(rnn)
             break
       NumberOfFish = len(xvalid)
       print(i,"density=",NumberOfFish*L*L*L/V,"fish=",NumberOfFish)

  fish = len(xvalid)
  d = np.zeros(fish)
  for i in range(fish):
      r = (xvalid[i]-xvalid)**2 + (yvalid[i]-yvalid)**2 + (zvalid[i]-zvalid)**2 
      r[i] = 1e6
      d[i] = np.min(r)**0.5/L

  fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
  axs.plot(dist_x,dist_y)
  axs.hist(d, bins=20, density=True)
  plt.show()


  return xvalid,yvalid,zvalid


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--a'   , required=True, type=float)
  parser.add_argument('--b'   , required=True, type=float)
  parser.add_argument('--c'   , required=True, type=float)
  parser.add_argument('--name', required=True, type=str  )
  parser.add_argument('--leader', required=True, type=int)
  args = vars(parser.parse_args())
  a = args['a']
  b = args['b']
  c = args['c']
  name = args['name']
  leader = args['leader']

  random.seed(0)
  np.random.seed(0)

  L = 0.2
  x,y,z = FishSamples(L*a,L*b,L*c,L)
  fish = len(x)

  x = 3.0 + np.asarray(x)
  y = 2.0 + np.asarray(y)
  z = 2.0 + np.asarray(z)
  idx = np.argmin(x)

  f = open(name+".sh", "w")
  f.write(\
"#!/bin/bash\n\
NNODE=16\n\
PSOLVER=\"iterative\"\n\
\n\
\n\
FACTORY=\n")
  for j in range(fish):
      leader1 = (j == idx) 
      if leader != 1:
         leader1 = True
      f.write('FACTORY+="StefanFish L='+str(L)+' T=1.0 xpos={} ypos={} zpos={} CorrectPosition=true CorrectPositionZ=true CorrectRoll=true bFixFrameOfRef={} \n\"\n'.format(x[j],y[j],z[j],leader1))

  # WRITE SOLVER SETTINGS
  f.write('\nOPTIONS=\n\
OPTIONS+=" -extentx 8.0"\n\
OPTIONS+=" -bpdx 8 -bpdy 4 -bpdz 4"\n\
OPTIONS+=" -tdump 0.1 -tend 100.0 "\n\
OPTIONS+=" -CFL 0.3 -nu 0.000008"\n\
OPTIONS+=" -levelMax 8 -levelMaxVorticity 7 -levelStart 3 -Rtol 1.0 -Ctol 0.1"\n\
OPTIONS+=" -bMeanConstraint 2 "\n\
OPTIONS+=" -poissonSolver ${PSOLVER}"\n')
