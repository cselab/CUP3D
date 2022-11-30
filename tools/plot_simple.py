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
  data = np.loadtxt(root+runname+"/forceValues_0.dat", skiprows=1)
  t = data[:,1]*speed/(radius)
  cd =  2*data[:,2] /(radius*radius*np.pi*speed*speed)
  plt.plot(t, cd, color='green',linestyle="--", label=name)

  data = np.loadtxt(root+runname+"/forceValues_penalization_0.dat", skiprows=1)
  t = data[:,1]*speed/(radius)
  cd =  -2*data[:,3] /(radius*radius*np.pi*speed*speed)
  plt.plot(t, cd, color='red',linestyle="--", label=name+" (pen)")

def dragWang( Re, t ):
  Re = int(Re)/2 # Wang uses radius-based Re
  return 12*np.sqrt(1/(np.pi*t*Re))+12/Re

if __name__ == '__main__':
  rootSCRATCH = "/scratch/project_465000158/chatzima/CubismUP3D/"
  cases= ["500"]
  runname = ["sphere500"]
  name  = ["Re=500"]
  speed = 0.125
  radius = 0.0625
  time = np.linspace( 1e-10, 1.0, 1001 )
  for i in range( len(runname) ):
    plt.plot(time, dragWang( cases[i], time ), linestyle="--", color=colors[i], label="Wang (1969)")
    plotDrag( rootSCRATCH, runname[i], speed, radius, i, cases[i], name[i])
  #plt.yscale("log")
  #plt.xscale("log")
  #plt.ylim(1e-6,1e+7)
  plt.ylim(1e-6,1e1)
  plt.xlim(1e-5,1.0)
  plt.xlabel("Time")
  plt.ylabel("Drag Coefficient")
  plt.legend(facecolor="white", edgecolor="white", ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.3))
  plt.grid(b=True, which='major', color="white", linestyle='-')
  plt.tight_layout()
  plt.savefig("s.png")
  plt.show()
