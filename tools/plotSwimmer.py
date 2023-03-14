import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
import seaborn as sns
import math
import scipy
sns.set_theme()
sns.set_style("whitegrid")

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']

def plotSpeed( root, runname ):
  data = np.loadtxt(root+runname+"/velocity_0.dat", skiprows=1)
  t = data[:,1]
  vel =  -data[:,9]
  plt.plot(t, vel, label=runname)

if __name__ == '__main__':
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  levels= np.arange(4,10)
  runnames = ["larvalFish-L9-LV{}".format(l) for l in levels]
  for runname in runnames:
    plotSpeed( rootSCRATCH, runname)

  # plt.ylim(1e-6,1e1)
  # plt.xlim(1e-5,1.0)
  plt.xlabel("Time")
  plt.ylabel("Velocity")
  plt.legend(facecolor="white", edgecolor="white", ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.3))
  plt.grid(b=True, which='major', color="white", linestyle='-')
  plt.tight_layout()
  plt.savefig("s.png")
  plt.show()
