import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
import seaborn as sns
sns.set_theme()

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

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
  data = np.loadtxt(root+runname+"/forceValues_surface_0.dat", skiprows=1)
  time = data[:,1]*speed/(2*radius)
  # presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
  # viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
  totDrag =  2*data[:,3] /(radius*radius*np.pi*speed*speed)

  # plt.plot(time, presDrag, color=lighten_color(colors[i],1.2))# , label=runname+", $C_p(t={})=${}".format(time[-1],presDrag[-1]))
  # plt.plot(time, viscDrag, color=lighten_color(colors[i],1.4))# , label=runname+", $C_v(t={})=${}".format(time[-1],viscDrag[-1]))
  plt.plot(time, totDrag, color=lighten_color(colors[i],1), label="present ({} levels)".format(i+4))

  # plt.ylim([0,1.5])
  plt.ylim([0,4])
  # plt.xlim([0,65])
  plt.xlim([0,0.5])
  plt.xlabel("Time $T=tu/D$")
  plt.ylabel("Drag Coefficient $C_D=2|F_x|/\pi r^2u^2$")
  plt.grid()
  # plt.title("$C_D(t=10)=${}".format(drag[-1]))

def plotValidation():
  # "The impulsive starting of a sphere", C.-Y. Wang
  def dragWang( Re, t ):
    t = t*2 #Wang uses t*(u/r), not t*(u/D)
    Re = int(Re)/2 # Wang uses radius-based Re
    return 12*np.sqrt(1/(np.pi*t*Re))+12/Re

  # Experimental Data from Roos et al. (1971)
  targetValues = { "300": 0.629, "500": 0.547, "1000": 0.472}

  rootSCRATCH = "/scratch/snx3000/mchatzim/CubismUP3D/"
  rootPROJECT = "/project/s929/pweber/CUP3D/unmollified-chi/"
  rootVALID = "/project/s929/pweber/sphereValidationData/"
  cases= [ "300", "500", "1000" ]
  levels = "4"
  poissonTol = "6"

  # runnames = [ "sphereRe{}_levels{}_poissonTol{}".format(case, levels, poissonTol) for case in cases]
  # runnames = [ "sphereRe{}_levels{}".format(case, levels, poissonTol) for case in cases]
  runnames = [ "sphereRe{}_poissonTol{}".format(case, poissonTol) for case in cases]

  for i in range( len(cases) ):
    # Plot analytical solution at initial times
    # time = np.linspace( 0, 0.5, 1001 )
    # plt.plot(time, dragWang( cases[i], time ), linestyle="--", color=lighten_color(colors[i],0.75), label="Wang (1969)")

    # Plot experimental value
    time = np.linspace( 0, 100, 1001 )
    Cd   = np.full(time.shape, targetValues[cases[i]])
    plt.plot(time, Cd, color=lighten_color(colors[i],0.5), linestyle="--", label="Ross (1971)")

    # Plot Ploumhans simulation data
    ploumhansData = np.loadtxt(rootVALID+"Re"+cases[i]+".txt", delimiter=",")
    plt.plot(ploumhansData[:,0], ploumhansData[:,1], marker="2", linestyle='None', color=lighten_color(colors[i],0.5), label="Ploumhans, Winckelmans, Salmon, Leonard, and Warren (2002)")

    if cases[i] == "300":
      # Plot Johnson simulation data
      johnsonData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Johnson.txt", delimiter=",")
      plt.plot(johnsonData[8:,0], johnsonData[8:,1], marker="+", linestyle='None', color=lighten_color(colors[i],0.5), label="Johnson and Patel (1998)")
      # Plot simulation results
      speed = 0.125
      radius = 0.0625
      runname = "Ploumhans-Sphere/re300"
      plotDrag( rootSCRATCH, runname, speed, radius, i, cases[i] )
      plt.xlim([0,30])
      # plotDrag( rootPROJECT, runnames[i], speed, radius, i, cases[i] )
    elif cases[i] == "500":
      speed = 1
      radius = 0.5
      runname = "Re"+cases[i]+"_2"
      plotDrag( rootSCRATCH, runname, speed, radius, i, cases[i] )
      plt.xlim([0,50])
      # plotDrag( rootPROJECT, runnames[i], speed, radius, i, cases[i] )

    elif cases[i] == "1000":
      # Plot Mimeau simulation data
      mimeauData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Mimeau.txt", delimiter=",")
      plt.plot(mimeauData[:,0], mimeauData[:,1], marker="+", linestyle='None', color=lighten_color(colors[i],0.5), label="Mimeau, Cottet, and Mortazavi (2016)")
      # Plot Spietz simulation data
      mimeauData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Spietz.txt", delimiter=",")
      plt.plot(mimeauData[:,0], mimeauData[:,1], marker=".", linestyle='None', color=lighten_color(colors[i],0.5), label="Spietz, Hejlesen, and Walther (2017)")
      # Plot simulation results
      speed = 1
      radius = 0.5
      runname = "Re"+cases[i]+"_6"
      plotDrag( rootSCRATCH, runname, speed, radius, i, cases[i] )
      plt.xlim([0,20])
      # plotDrag( rootPROJECT, runnames[i], speed, radius, i, cases[i] )

  
    # plt.legend(facecolor="white", edgecolor="white", ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.3))
    plt.legend()
    plt.title("Re={}".format(cases[i]))
    plt.grid()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
  # "The impulsive starting of a sphere", C.-Y. Wang
  def dragWang( Re, t ):
    t = t*2 #Wang uses t*(u/r), not t*(u/D)
    Re = int(Re)/2 # Wang uses radius-based Re
    return 12*np.sqrt(1/(np.pi*t*Re))+12/Re

  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  cases= [ "1000" ]
  levels = np.arange(4,8)
  timesteps = [ "1e-4", "2e-4", "5e-4", "1e-3", "2e-3" ]

  # runnames = [ "sphereRe{}_levels{}".format(cases[0], level) for level in levels]
  runnames = [ "sphereRe{}_levels5_dt{}".format(cases[0], timestep) for timestep in timesteps]

  # Plot analytical solution at initial times
  time = np.linspace( 0, 0.5, 5001 )
  plt.plot(time, dragWang( cases[0], time ), linestyle="--", color="black", label="Wang (1969)")
  
  # Get target value for refinement study
  lowerBound = 0.0001
  upperBound = 0.05
  indices = (time > lowerBound) & (time < upperBound)
  target = np.mean( dragWang( cases[0], time )[indices] )

  # container for QoI
  meanPressureDrag = []
  meanViscousDrag = []
  meanDrag = []
  h = []
  for i in range( len(runnames) ):
    # Plot simulation results
    speed = 0.125
    radius = 0.0625
    plotDrag( rootSCRATCH, runnames[i], speed, radius, i, cases[0] )
    plt.xlim([0,0.5])
    # Compute QoI
    data = np.loadtxt(rootSCRATCH+runnames[i]+"/forceValues_surface_0.dat", skiprows=1)
    time = data[:,1]*speed/(2*radius)
    # presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
    # viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
    totDrag =  2*data[:,3] / (radius*radius*np.pi*speed*speed)

    indices = (time > lowerBound) & (time < upperBound)
    # meanPressureDrag.append( np.mean(presDrag[indices]) )
    # meanViscousDrag.append( np.mean(viscDrag[indices]) )
    print(len(totDrag[indices]))
    meanDrag.append( np.mean(totDrag[indices]) )

    #compute minimal gridspacing at that level
    # h.append( 4.0/(32*8*2**(int(levels[i])-1)) )

  plt.legend()
  plt.title("Re={}".format(cases[0]))
  plt.grid(b=True, which='minor', color="white", linestyle='-')
  plt.tight_layout()
  plt.show()


  h = np.array( timesteps ).astype(np.float)
  # plt.plot( h, np.abs( meanPressureDrag), "ko" )
  # plt.plot( h, np.abs( meanViscousDrag ), "go" )
  plt.plot( h, np.abs( meanDrag - target), "o" )
  # plt.plot( h, 10**1*h**(1), label="1st order", linestyle="--" )
  plt.plot( h, 7*10**2*h, label="1st order", linewidth=1, linestyle="--" )
  plt.plot( h, 7*10**4*h**(2), label="2nd order", linewidth=1, linestyle="--" )
  # plt.plot( h, 10**7*h**(3),   label="3rd order", linewidth=1, linestyle="--" )
  # plt.xlim( [5e4,1e6])
  # plt.xlabel("Gridspacing")
  plt.xlabel("Timestep")
  plt.ylabel("Error")
  plt.xscale("log")
  plt.yscale("log")
  plt.legend()
  plt.show()
