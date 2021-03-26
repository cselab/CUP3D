import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys

colors = ['C0', 'C1', 'C2', 'C3']

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], min(1, amount * c[1]), c[2])

def plotDrag( root, runname, speed, radius, i ):
  data = np.loadtxt(root+runname+"/forceValues_surface_0.dat", skiprows=1)
  time = data[:,1]
  presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
  viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
  totDrag =  2*data[:,3] /(radius*radius*np.pi*speed*speed)

  # plt.plot(time, presDrag, color=lighten_color(colors[i],1.2) , label=runname+", $C_p$")
  # plt.plot(time, viscDrag, color=lighten_color(colors[i],1.4) , label=runname+", $C_v$")
  plt.plot(time, totDrag,  color=lighten_color(colors[i],1) ,   label=runname+", $C_D(t=100)=${}".format(totDrag[-1]))
  plt.ylim([0,1])
  # plt.xlim([0,1])
  plt.xlabel("time $t$")
  plt.ylabel("Drag Coefficients $C=2F/\pi r^2u_\infty^2$")
  plt.grid()
  # plt.title("$C_D(t=10)=${}".format(drag[-1]))

if __name__ == '__main__':
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  rootPROJECT = "/project/s929/pweber/"
  runname = ["SphereRe300", "SphereRe300_levels5","SphereRe300_levels6"]
  #runname = ["SphereRe500", "SphereRe500_levels5","SphereRe500_levels6"]
  #runname = ["SphereRe1000_levels4", "SphereRe1000_levels5","SphereRe1000_levels6"]
  speed = 0.125
  radius = 0.0625
  for i in range( len(runname) ):
    plotDrag( rootSCRATCH, runname[i], speed, radius, i )
  plt.legend()
  plt.show()
