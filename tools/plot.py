import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
import seaborn as sns
sns.set_theme()

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
  # data = np.loadtxt(root+runname+"/forceValues_surface_0.dat", skiprows=1)
  data = np.loadtxt(root+runname+"/forceValues_0.dat", skiprows=1)
  time = data[:,1]
  presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
  viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
  totDrag =  -2*data[:,3] /(radius*radius*np.pi*speed*speed)

  plt.plot(time, presDrag, color=lighten_color(colors[i],1.2))# , label=runname+", $C_p(t={})=${}".format(time[-1],presDrag[-1]))
  plt.plot(time, viscDrag, color=lighten_color(colors[i],1.4))# , label=runname+", $C_v(t={})=${}".format(time[-1],viscDrag[-1]))
  plt.plot(time, totDrag,  color=lighten_color(colors[i],1) ,   label=runname+", $C_D(t={:.1f})=${:.3f}".format(time[-1],totDrag[-1]))
  plt.ylim([0,1])
  # plt.xlim([0,1])
  plt.xlabel("time $t$")
  plt.ylabel("Drag Coefficients $C=2F/\pi r^2u_\infty^2$")
  plt.grid()
  # plt.title("$C_D(t=10)=${}".format(drag[-1]))

if __name__ == '__main__':
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/old/"
  rootPROJECT = "/project/s929/pweber/"
  case="Re1000"
  # Experimental Data from Roos et al. (1971)
  targetValues = { "Re300": 0.629, "Re500": 0.547, "Re1000": 0.472}
  runname = ["Sphere"+case+"_levels04", "Sphere"+case+"_levels05","Sphere"+case+"_levels06"] #, "Sphere"+case+"_levels07"]
  speed = 0.125
  radius = 0.0625
  for i in range( len(runname) ):
    plotDrag( rootSCRATCH, runname[i], speed, radius, i )
  plt.hlines(y=targetValues[case], xmin=0, xmax=100, color="black", linestyles="dashed", label="$C_D=${:.3f} [Roos et al. (1971)]".format(targetValues[case]), zorder=10)
  plt.legend()
  plt.grid()
  plt.show()
