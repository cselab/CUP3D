import numpy as np
import scipy as sp
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
  data = np.loadtxt(root+runname+"/forceValues_surface_0.dat", skiprows=1)
  time = data[:,1]*speed/(2*radius)
  # time = data[:,1]*speed/(radius)
  presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
  viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
  totDrag  = 2*data[:,3] /(radius*radius*np.pi*speed*speed)

  # plt.plot(time, presDrag-1/3*totDrag, color=lighten_color(colors[i],1.2), label="$C_p-1/3C_D$" )# , label=runname+", $C_p(t={})=${}".format(time[-1],presDrag[-1]))
  # plt.plot(time, viscDrag-2/3*totDrag, color=lighten_color(colors[i],1.4), label="$C_v-2/3C_D$")# , label=runname+", $C_v(t={})=${}".format(time[-1],viscDrag[-1]))
  # plt.plot(time, totDrag, color=lighten_color(colors[i],1), label="present ({} levels)".format(7-i))
  plt.plot(time, totDrag, color=lighten_color(colors[i],1), label="present ({} levels)".format(i+3))#, label="present (dt = {})".format(i+1))

  plt.xlabel("Time $T=tu/D$")
  # plt.xlabel("Time $T=tu/r$")
  plt.ylabel("Drag Coefficient $C_D=2|F_x|/\pi r^2u^2$")


def plotValidation():
  # "The impulsive starting of a sphere", C.-Y. Wang
  def dragWang( Re, t ):
    # t = t*2 #Wang uses t*(u/r), not t*(u/D)
    Re = int(Re)/2 # Wang uses radius-based Re
    return 12*np.sqrt(1/(np.pi*t*Re))+12/Re

  # Experimental Data from Roos et al. (1971)
  targetValues = { "300": 0.629, "500": 0.547, "1000": 0.472}

  # rootSCRATCH = "/scratch/snx3000/mchatzim/CubismUP3D/"
  # rootSCRATCH = "/project/s929/mchatzim/SphereValidation/"
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  rootPROJECT = "/project/s929/pweber/CUP3D/unmollified-chi/"
  rootVALID = "/project/s929/pweber/sphereValidationData/"
  cases= [ "1000" ]
  levels = "4"
  poissonTol = "6"

  # runnames = [ "sphereRe{}_levels{}_poissonTol{}".format(case, levels, poissonTol) for case in cases]
  # runnames = [ "sphereRe{}_levels{}".format(case, levels, poissonTol) for case in cases]
  # runnames = [ "sphereRe{}_poissonTol{}".format(case, poissonTol) for case in cases]

  for i in range( len(cases) ):
    # Plot analytical solution at initial times
    # time = np.linspace( 0, 1, 1001 )
    # plt.plot(time, dragWang( cases[i], time ), linestyle="--", color=lighten_color(colors[i],0.75), label="Wang (1969)")

    # Plot experimental value
    # time = np.linspace( 0, 200, 1001 )
    # Cd   = np.full(time.shape, targetValues[cases[i]])
    # plt.plot(time, Cd, color=lighten_color(colors[i],0.5), linestyle="--", label="Ross (1971)")

    # Plot Ploumhans simulation data
    # ploumhansData = np.loadtxt(rootVALID+"Re"+cases[i]+".txt", delimiter=",")
    # plt.plot(ploumhansData[:,0], ploumhansData[:,1], marker="2", linestyle='None', color=lighten_color(colors[i],0.5), label="Ploumhans et al. (2002)") # , Winckelmans, Salmon, Leonard, and Warren

    if cases[i] == "300":
      # Plot Johnson simulation data
      # johnsonData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Johnson.txt", delimiter=",")
      # plt.plot(johnsonData[8:,0], johnsonData[8:,1], marker="+", linestyle='None', color=lighten_color(colors[i],0.5), label="Johnson and Patel (1998)")
      # Plot Spietz simulation data
      # spietzData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Spiez.txt", delimiter=",")
      # plt.plot(spietzData[8:,0], spietzData[8:,1], marker="--", linestyle='None', color=lighten_color(colors[i],0.5), label="Spietz et al. (2017)") # , Hejlesen, and Walther
      # Plot simulation results
      speed = 0.125
      radius = 0.0625
      runnames = [  "/Re300_long_levels3_big_domain"] #"diskRe300_levels4_breakSymmetry",
      for runname in runnames:
        plotDrag( rootSCRATCH, runname, speed, radius, i, cases[i] )
      plt.xlim([0,200])
      plt.ylim([0.65,0.75])
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
      # mimeauData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Mimeau.txt", delimiter=",")
      # plt.plot(mimeauData[:,0], mimeauData[:,1], marker="+", linestyle='None', color=lighten_color(colors[i],0.5), label="Mimeau et al. (2016)") #, Cottet, and Mortazavi
      # Plot Spietz simulation data
      # mimeauData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Spietz.txt", delimiter=",")
      # plt.plot(mimeauData[:,0], mimeauData[:,1], marker=".", linestyle='None', color=lighten_color(colors[i],0.5), label="Spietz et al. (2017)") # , Hejlesen, and Walther
      # Plot simulation results
      speed = 0.125
      radius = 0.0625
      # runname = "re"+cases[i]
      runnames = [ "sphereRe1000_levels{}".format(level) for level in np.arange(3,8) ]
      for j, runname in enumerate(runnames):
        plotDrag( rootSCRATCH, runname, speed, radius, j, cases[i] )
        plt.xlim([0,2])
        plt.ylim([0,3])
      # plotDrag( rootPROJECT, runnames[i], speed, radius, i, cases[i] )

  
    # plt.legend(facecolor="white", edgecolor="white", ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.3))
    plt.legend()
    # plt.title("Re={}".format(cases[i]))
    # plt.yscale("log")
    plt.grid(b=True, which='minor', color="white", linestyle='-')
    plt.tight_layout()
    plt.show()

import scipy.signal as signal

# def computeStrouhal():
#   speed = 0.125
#   radius = 0.0625

#   data = np.loadtxt("/scratch/snx3000/mchatzim/CubismUP3D/Re300_long_levels3_big_domain/forceValues_surface_0.dat", skiprows=1)
#   time = data[:,1]*speed/(2*radius)
#   totDrag  = 2*data[:,3] /(radius*radius*np.pi*speed*speed)

#   fDragTarget = sp.interpolate.interp1d(time, totDrag)

#   time = np.linspace(100,time[-1], 1001)
#   y = sp.fftpack.fft(fDragTarget(time))

#   plt.plot(sp.sp.sp.y)
#   plt.show()

def plotVortexSheet():
  data = np.loadtxt("/scratch/snx3000/mchatzim/CubismUP3D/debug6/7/gamma_integral_000000020.txt")
  theta = data[:,0]
  gamma = data[:,1]
  speed = 0.125

  tick_pos= np.linspace(0,np.pi,5)
  labels = ['0','$\pi$/8','$\pi$/4','$\pi$/2','$\pi$']
  plt.xticks(tick_pos, labels)

  plt.plot(theta, -1.5*np.sin(theta), "k--", label="analytic")
  plt.plot(theta, gamma, ".", label="present")
  plt.xlabel("Angle $\\theta$")
  plt.ylabel("Average Vortex Sheet Strength $\\frac{1}{2\pi u}\int \gamma d\phi$")
  plt.legend()
  plt.show()

def plotRefinement():
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/poissonTol7/"
  cases= [ "1000" ]
  levels = np.arange(3,7)
  levels = levels[::-1]
  timesteps = [ "1e-2.25", "1e-2.5", "1e-2.75", "1e-3", "1e-3.25", "1e-3.5"]
  refinement = "vorticity"

  if refinement != "vorticity":
    if refinement == "space":
      runnames = [ "sphereRe{}_levels{}".format(cases[0], level) for level in levels]
    elif refinement == "time":
      runnames = [ "sphereRe{}_levels5_dt{}".format(cases[0], timestep) for timestep in timesteps]
    
    # Get integration bounds for refinement study
    lowerBound = 0
    upperBound = 8

    # container for QoI
    meanGridpoints = []
    meanErrorDrag = []
    for i in range( len(runnames) ):
      # Simulation parameters
      speed = 0.125
      radius = 0.0625

      # Get drag
      data = np.loadtxt(rootSCRATCH+runnames[i]+"/forceValues_surface_0.dat", skiprows=1)
      time = data[:,1]*speed/radius
      presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
      viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
      totDrag =  2*data[:,3] /(radius*radius*np.pi*speed*speed)

      # Get number of blocks
      blockData = np.loadtxt(rootSCRATCH+runnames[i]+"/diagnostics.dat", skiprows=1)
      blockTime = blockData[:,1]
      numBlocks = blockData[:,-1]

      # get samples according to integration bounds
      indices = (time > lowerBound) & (time < upperBound)
      timeDiff = upperBound-lowerBound
      time = time[indices]

      print("Number of Samples: ", len(presDrag[indices]))

      if refinement == "space":
        if levels[i] == 6:
          dragTarget = totDrag[indices]
          fDragTarget = sp.interpolate.interp1d(time, dragTarget)
        else:
          # errorDrag = np.abs( (totDrag[indices] - dragTarget ) )
          errorDrag = np.abs( (totDrag[indices] - fDragTarget(time) ) )

          #plot instantaneous error
          plt.plot( time, errorDrag, label="{} levels".format(levels[i]) )

          #compute number of gridpoints
          meanGridpoints.append( sp.integrate.simps( numBlocks*8**3, blockTime) )

          #compute mean error drag
          meanErrorDrag.append( sp.integrate.simps( errorDrag, time) )
      elif refinement == "time":
        if timesteps[i] == "1e-4":
          dragTarget = totDrag[indices]
          fDragTarget = sp.interpolate.interp1d(time, dragTarget)
        else:
          meanErrorDrag.append( sp.integrate.simps( np.abs( totDrag[indices]  - fDragTarget(time) ), time ) )

    plt.xlabel("Time $T=tu/r$")
    plt.ylabel("Instantaneous Error $\epsilon(t;l)$")
    plt.legend()
    plt.yscale("log")
    plt.grid(b=True, which='minor', color="white", linestyle='-')
    plt.show()
  else:
    ## data for vorticity for 
    #t=0.5
    # from old runs
    # h = np.array([1658880, 1803968, 2902272, 3879360, 8062336, 28699840 ])
    # vorticityDiff = [ 0.034554, 0.0264363, 0.0156306, 0.00731397, 0.00309227, 0.00130656]
    # h = np.array([1435638, 2281472, 6152192, 29620224 ])
    # vorticityDiff = [ 0.0106758, 0.00469273, 0.00176592, 0.000703645]
    h = np.array([1349632, 2080768, 5062656, 21735424 ])
    vorticityDiff = [ 0.0095116, 0.00470031, 0.00180059, 0.000749308]

  ## Plot Order ##
  if refinement == "space":
    meanGridpoints = np.array( meanGridpoints )
    plt.plot( meanGridpoints, meanErrorDrag, "o")
    plt.plot( meanGridpoints, 2*10**0*meanGridpoints**(-1/3), label="1st order", linewidth=1, linestyle="--" )
    plt.plot( meanGridpoints, 2*10**1*meanGridpoints**(-2/3), label="2nd order", linewidth=1, linestyle="--" )
    # plt.ylim([0.01,1])
    # plt.xlim([10**3,2*10**4])
  elif refinement == "time":
    h = np.array( timesteps[1:] ).astype(np.float)
    plt.plot( h, meanErrorDrag, "o")
    plt.plot( h, 2*10**3*h, label="1st order", linewidth=1, linestyle="--" )
    plt.plot( h, 1.5*10**5*h**(2), label="2nd order", linewidth=1, linestyle="--" )
    # plt.ylim([0.02,5])
    plt.xlabel("Timestep")
  elif refinement == "vorticity":
    plt.plot( h, vorticityDiff, "o")
    plt.plot( h, 2*10**-1*h**(-1/3), label="1st order", linewidth=1, linestyle="--" )
    plt.plot( h, 6*10**1*h**(-2/3), label="2nd order", linewidth=1, linestyle="--" )
    # plt.ylim([0.001,0.02])
    # plt.xlim([10**6,5*10**7])
    plt.xlabel("Number of Gridpoints")

  plt.xlabel("Number of Gridpoints")
  plt.ylabel("Error")
  plt.xscale("log") #, base=2)
  plt.yscale("log")
  plt.legend()
  plt.grid(b=True, which='minor', color="white", linestyle='-')
  plt.show()

if __name__ == '__main__':
  # plotValidation()
  # plotVortexSheet()
  plotRefinement()
  # computeStrouhal()
