import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
import seaborn as sns
import math
import scipy
sns.set_theme()
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']

def plotDrag( root, runname, speed, radius, i, label , name):
  data = np.loadtxt(root+runname+"/forceValues_surface_0.dat", skiprows=1)
  t = data[:,1]*speed/(radius)
  cd =  2*data[:,3] /(radius*radius*np.pi*speed*speed)
  plt.plot(t, cd, color=colors[i], label=name+" (S)")
  data = np.loadtxt(root+runname+"/forceValues_0.dat", skiprows=1)
  t = data[:,1]*speed/(radius)
  cd =  -2*data[:,3] /(radius*radius*np.pi*speed*speed)
  #plt.plot(t, cd, color=colors[i],linestyle="--", label=name+" (N)")

def dragWang( Re, t ):
  Re = int(Re)/2 # Wang uses radius-based Re
  return 12*np.sqrt(1/(np.pi*t*Re))+12/Re

if __name__ == '__main__':
  rootSCRATCH = "/scratch/snx3000/mchatzim/CubismUP3D/"
  '''
  cases= ["300","500","1000"]
  runname = ["Ploumhans-Sphere/re300","Ploumhans-Sphere/re500","Ploumhans-Sphere/re1000_fine"]
  name  = ["Re=300","Re=500","Re=1000"]
  speed = 0.125
  radius = 0.0625
  '''
  cases= ["500","1000"]
  runname = ["Re500_2","Re1000_3"]
  name  = ["Re=500","Re=1000"]
  speed = 1.
  radius = 0.5
  time = np.linspace( 1e-10, 0.1, 1000001 )

  for i in range( len(runname) ):
    plt.plot(time, dragWang( cases[i], time ), linestyle="--", color=colors[i], label="Wang (1969)")
    plotDrag( rootSCRATCH, runname[i], speed, radius, i, cases[i], name[i])
  plt.yscale("log")
  plt.xscale("log")
  plt.ylim(1e-6,1e+7)
  #plt.ylim(1e-6,2.0)
  plt.xlim(1e-5,0.1)
  plt.xlabel("Time")
  plt.ylabel("Drag Coefficient")
  plt.legend(facecolor="white", edgecolor="white", ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.3))
  plt.grid(b=True, which='major', color="white", linestyle='-')
  plt.tight_layout()
  plt.show()
