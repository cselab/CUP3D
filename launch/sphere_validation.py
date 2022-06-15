import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
import seaborn as sns
import math
import scipy
sns.set_theme()

speed = 0.125
radius = 0.0625
Re = 500

def plotDrag(root):
  data = np.loadtxt(root+"/forceValues_0.dat", skiprows=1)
  t = data[:,1]*speed/(radius)
  cd =  2*data[:,2] /(radius*radius*np.pi*speed*speed)
  #cd =  2*data[:,8] /(radius*radius*np.pi*speed*speed)
  #cd =  2*data[:,11] /(radius*radius*np.pi*speed*speed)
  plt.plot(t, cd)

  #data = np.loadtxt(root+"/forceValues_surface_0.dat", skiprows=1)
  #t = data[:,1]*speed/(radius)
  #cd =  2*data[:,3] /(radius*radius*np.pi*speed*speed)
  #cd =  2*data[:,9] /(radius*radius*np.pi*speed*speed)
  #cd =  2*data[:,12] /(radius*radius*np.pi*speed*speed)
  plt.plot(t, cd)

def dragWang( Re, t ):
  Re = int(Re)/2 # Wang uses radius-based Re
  return (12*np.sqrt(1/(np.pi*t*Re))+12/Re)*(3.0/3.0)

if __name__ == '__main__':
  root = "/project/s929/mchatzim/SphereValidation/Re1000"
  root = "/scratch/snx3000/mchatzim/CubismUP3D/sphere/fine"
  time = np.linspace(1e-10, 1.0, 1001)
  plt.plot(time, dragWang( Re, time ), label="Wang (1969)")
  plotDrag(root)
  plt.ylim(0,3.0)
  plt.xlim(0,1.0)
  plt.xlabel("Time")
  plt.ylabel("Drag Coefficient")
  plt.legend()
  plt.tight_layout()
  plt.show()
