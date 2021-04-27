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

def plotDrag( root, runname, speed, radius, i, label ):
  

  # data = np.loadtxt(root+runname+"/forceValues_surface_0.dat", skiprows=1)
  data = np.loadtxt(root+runname+"/forceValues_0.dat", skiprows=1)
  time = data[:,1]*speed/(radius)
  # presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
  # viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
  totDrag =  -2*data[:,3] /(radius*radius*np.pi*speed*speed)

  # plt.plot(time, presDrag, color=lighten_color(colors[i],1.2))# , label=runname+", $C_p(t={})=${}".format(time[-1],presDrag[-1]))
  # plt.plot(time, viscDrag, color=lighten_color(colors[i],1.4))# , label=runname+", $C_v(t={})=${}".format(time[-1],viscDrag[-1]))
  plt.plot(time, totDrag, color=lighten_color(colors[i],1), label="present")

  # plt.ylim([0,1.5])
  plt.ylim([0,4])
  # plt.xlim([0,40])
  plt.xlim([0,0.5])
  plt.xlabel("Time $T=tu/r$")
  plt.ylabel("Drag Coefficient $C_D=2|F_x|/\pi r^2u^2$")
  plt.grid()
  # plt.title("$C_D(t=10)=${}".format(drag[-1]))

if __name__ == '__main__':
  # "Numerical Solutions for Time-Dependent Flow Past an Impulsively Started Sphere", S. C. R. Dennis and J. D. A. Walker (1972)
  # def dragDennisWalker( Re, t ):
  #   k = 2*np.sqrt(2*t/Re)
  #   fricDrag = 16/(3*Re)*G1(0,t)
  #   presDrag = -8/(3*Re)*(G1(0,t)+dG1dx(0,t))

  # "The impulsive starting of a sphere", C.-Y. Wang
  def dragWang( Re, t ):
    Re = int(Re)/2 # Wang uses radius-based Re
    return 12*np.sqrt(1/(np.pi*t*Re))+12/Re

  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  rootPROJECT = "/project/s929/pweber/CUP3D/"
  cases= [ "300", "500", "1000" ] #"300", 
  # case = "1000"
  # for case in cases:
  levels = "4"
  poissonTol = "6"
  # Experimental Data from Roos et al. (1971)
  targetValues = { "300": 0.629, "500": 0.547, "1000": 0.472}
  # runname = ["sphere"+case+"_levels5_poissonTol"+poissonTol, "sphere"+case+"_levels6_poissonTol"+poissonTol,"sphere"+case+"_levels7_poissonTol"+poissonTol] #, "Sphere"+case+"_levels07"]
  # runname = ["sphere"+case+"_levels"+levels+"_poissonTol5","sphere"+case+"_levels"+levels+"_poissonTol6", "sphere"+case+"_levels"+levels+"_poissonTol7"]
  runname = ["sphereRe"+case+"_poissonTol"+poissonTol for case in cases]
  # runname = ["sphereRe"+case+"_poissonTol5", "sphereRe"+case+"_poissonTol6"]
  speed = 0.125
  radius = 0.0625
  for i in range( len(runname) ):
    time = np.linspace( 0, 0.5, 1001 )
    plt.plot(time, dragWang( cases[i], time ), linestyle="--", color=lighten_color(colors[i],0.75), label="Wang (1969), Re={}".format(cases[i]))
    # time = np.linspace( 0, 200, 1001 )
    # Cd   = np.full(time.shape, targetValues[cases[i]])
    # plt.plot(time, Cd, color=lighten_color(colors[i],0.75), linestyle="--", label="Roos et al. (1971), Re={}".format(cases[i]))
    plotDrag( rootSCRATCH, runname[i], speed, radius, i, cases[i] )
  
  plt.legend(facecolor="white", edgecolor="white", ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.3))
    # plt.title("Re={}".format(cases[i]))
  plt.grid()
  plt.tight_layout()
  plt.show()
