import numpy as np
import matplotlib.pyplot as plt


def plotDrag( root, runname, speed, radius ):
  data = np.loadtxt(root+runname+"/forceValues_0.dat", skiprows=1)
  time = data[:,1]
  drag = -2*data[:,3]/(radius*radius*np.pi*speed*speed)
  plt.plot(time, drag, label=runname+", $C_D(t=100)=${}".format(drag[-1]))
  plt.ylim([0,2.5])
  # plt.xlim([0,10])
  plt.xlabel("time $t$")
  plt.ylabel("Drag Coefficient $C_D=2F/\pi r^2u_\infty^2$")
  plt.grid()
  # plt.title("$C_D(t=10)=${}".format(drag[-1]))

if __name__ == '__main__':
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  rootPROJECT = "/project/s929/pweber/"
  #runname = ["SphereRe300", "SphereRe300_levels5","SphereRe300_levels6"]
  runname = ["SphereRe500", "SphereRe500_levels5","SphereRe500_levels6"]
  #runname = ["SphereRe1000_levels4", "SphereRe1000_levels5","SphereRe1000_levels6"]
  speed = 0.125
  radius = 0.0625
  for i in range( len(runname) ):
    plotDrag( rootSCRATCH, runname[i], speed, radius )
  plt.legend()
  plt.show()
