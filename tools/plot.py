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

  ### Plot components of Drag ###
  # plt.plot(time, presDrag, color=lighten_color(colors[i],1.2))# , label=runname+", $C_p(t={})=${}".format(time[-1],presDrag[-1]))
  # plt.plot(time, viscDrag, color=lighten_color(colors[i],1.4))# , label=runname+", $C_v(t={})=${}".format(time[-1],viscDrag[-1]))
  ###############################

  # plt.plot(time, totDrag, color=lighten_color(colors[i],1), label="present (Re = {})".format(label) )
  plt.plot(time, totDrag, color=lighten_color(colors[i],1), label="present ({} levels)".format(i+3))

  plt.xlabel("Time $T=tu/D$")
  # plt.xlabel("Time $T=tu/r$")
  plt.ylabel("Drag Coefficient $C_D=2|F_x|/\pi r^2u^2$")


### "The impulsive starting of a sphere", C.-Y. Wang ###
def dragWang( Re, t ):
  t = t*2 #Wang uses t*(u/r), not t*(u/D)
  Re = int(Re)/2 # Wang uses radius-based Re
  return 12*np.sqrt(1/(np.pi*t*Re))+12/Re
########################################################

def plotValidation():
  ### Experimental Data from Roos et al. (1971) ###
  targetValues = { "300": 0.629, "500": 0.547, "1000": 0.472}
  #################################################

  ### Simulation parameters ###
  speed = 0.125
  radius = 0.0625
  #############################

  ### Path to data ###
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  rootVALID = "/project/s929/pweber/sphereValidationData/"
  ####################

  ### Re to plot ###
  cases= [ "1000" ]
  ##################

  ### Legend Hack for initial times ###
  # time = np.linspace( 0, 1, 1001 )
  # plt.plot(time, dragWang( cases[0], time )+5, linestyle="--", color="black", label="Wang (1969)")
  #####################################

  for i in range( len(cases) ):
    ### Plot analytic solution by Wang at initial times ###
    # time = np.linspace( 0, 1, 1001 )
    # plt.plot(time, dragWang( cases[i], time ), linestyle="--", color=lighten_color(colors[i],0.75))
    #######################################################

    ### Plot experimental value by Ross ###
    # time = np.linspace( 0, 200, 1001 )
    # Cd   = np.full(time.shape, targetValues[cases[i]])
    # plt.plot(time, Cd, color=lighten_color(colors[i],0.5), linestyle="--", label="Ross (1971)")
    #######################################

    ### Plot simulation data by Ploumhans ###
    # ploumhansData = np.loadtxt(rootVALID+"Re"+cases[i]+".txt", delimiter=",")
    # plt.plot(ploumhansData[:,0], ploumhansData[:,1], marker="2", linestyle='None', color=lighten_color(colors[i],0.5), label="Ploumhans et al. (2002)") # , Winckelmans, Salmon, Leonard, and Warren
    #########################################################

    if cases[i] == "300":
      ###  Plot simulation data by Johnson ###
      # johnsonData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Johnson.txt", delimiter=",")
      # plt.plot(johnsonData[8:,0], johnsonData[8:,1], marker="+", linestyle='None', color=lighten_color(colors[i],0.5), label="Johnson and Patel (1998)")
      ########################################

      ### Plot simulation data by Spietz ###
      # spietzData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Spiez.txt", delimiter=",")
      # plt.plot(spietzData[8:,0], spietzData[8:,1], marker="--", linestyle='None', color=lighten_color(colors[i],0.5), label="Spietz et al. (2017)") # , Hejlesen, and Walther
      #######################################

      ### Plot present results ###
      runnames = [ "sphereRe300_levels6", "sphereRe300_levels7" ]
      for runname in runnames:
        plotDrag( rootSCRATCH, runname, speed, radius, i, cases[i] )
      ############################
    elif cases[i] == "1000":
      ### Plot simulation data by Mimeau ###
      # mimeauData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Mimeau.txt", delimiter=",")
      # plt.plot(mimeauData[:,0], mimeauData[:,1], marker="+", linestyle='None', color=lighten_color(colors[i],0.5), label="Mimeau et al. (2016)") #, Cottet, and Mortazavi
      ######################################

      ### Plot simulation data by Spietz ###
      # spietzData = np.loadtxt(rootVALID+"Re"+cases[i]+"-Spietz.txt", delimiter=",")
      # plt.plot(spietzData[:,0], spietzData[:,1], marker=".", linestyle='None', color=lighten_color(colors[i],0.5), label="Spietz et al. (2017)") # , Hejlesen, and Walther
      ######################################

      ### Plot present results ###
      runnames = ["sphereRe1000_levels3_Euler"] #+ [ "sphereRe1000_levels{}".format(level) for level in np.arange(4,8) ]
      for j, runname in enumerate(runnames):
        plotDrag( rootSCRATCH, runname, speed, radius, j, cases[i] )
      ############################
  
  # plt.legend(facecolor="white", edgecolor="white", ncol=2 )#, loc="lower center", bbox_to_anchor=(0.5, -0.3))
  plt.legend()
  plt.xlim([0,0.5])
  plt.ylim([0,50])
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

def plotVerification():
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  cases= [ "1000" ]
  levels = np.arange(3,8)
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
    upperBound = 0.5

    # container for QoI
    meanGridpoints = []
    meanErrorDrag = []
    for i in range( len(runnames) ):
      # Simulation parameters
      speed = 0.125
      radius = 0.0625

      # Get drag
      data = np.loadtxt(rootSCRATCH+runnames[i]+"/forceValues_surface_0.dat", skiprows=1)
      time = data[:,1] #*speed/radius
      # presDrag = 2*data[:,9] /(radius*radius*np.pi*speed*speed)
      # viscDrag = 2*data[:,12]/(radius*radius*np.pi*speed*speed)
      totDrag =  2*data[:,3] /(radius*radius*np.pi*speed*speed)

      # Get number of blocks
      blockData = np.loadtxt(rootSCRATCH+runnames[i]+"/diagnostics.dat", skiprows=1)
      blockTime = blockData[:,1]
      numBlocks = blockData[:,-1]

      # get samples according to integration bounds
      indices = (time > lowerBound) & (time < upperBound)
      timeDiff = upperBound-lowerBound
      time = time[indices]

      # print("Number of Samples: ", len(totDrag[indices]))

      if refinement == "space":
        if levels[i] == 7:
          dragTarget = totDrag[indices]
          # dragTarget = dragWang(1000, time)
          # fDragTarget = sp.interpolate.interp1d(time, dragTarget)
        else:
          errorDrag = np.abs( (totDrag[indices] - dragTarget ) )
          # errorDrag = np.abs( (totDrag[indices] - fDragTarget(time) ) )

          ### plot instantaneous error ###
          # plt.plot( time, errorDrag, label="{} levels".format(levels[i]) )

          #compute number of gridpoints
          meanGridpoints.append( sp.integrate.simps( numBlocks*8**3, blockTime) / timeDiff )

          #compute mean error drag
          meanErrorDrag.append( sp.integrate.simps( errorDrag, time) / timeDiff )
      elif refinement == "time":
        if timesteps[i] == "1e-4":
          dragTarget = totDrag[indices]
          fDragTarget = sp.interpolate.interp1d(time, dragTarget)
        else:
          meanErrorDrag.append( sp.integrate.simps( np.abs( totDrag[indices]  - fDragTarget(time) ), time ) )

    # plt.xlabel("Time $T=tu/D$")
    # plt.ylabel("Instantaneous Error $\epsilon(t;l)$")
    # plt.legend()
    # plt.yscale("log")
    # plt.grid(b=True, which='minor', color="white", linestyle='-')
    # plt.show()
  else:
    ## data for vorticity for 
    #t=0.5
    # from old runs
    # h = np.array([1658880, 1803968, 2902272, 3879360, 8062336, 28699840 ])
    # vorticityDiff = [ 0.034554, 0.0264363, 0.0156306, 0.00731397, 0.00309227, 0.00130656]
    # h = np.array([1435638, 2281472, 6152192, 29620224 ])
    # vorticityDiff = [ 0.0106758, 0.00469273, 0.00176592, 0.000703645]
    # run with high PT
    # h = np.array([1349632, 2080768, 5062656, 21735424 ])
    # vorticityDiff = [ 0.0095116, 0.00470031, 0.00180059, 0.000749308]
    # run with equal times
    h = np.array([1335296, 2080768, 4962304, 20860928 ])
    vorticityDiff = [ 0.0101682, 0.00494003, 0.00208448, 0.000849824]

  ## Plot Order ##
  if refinement == "space":
    meanGridpoints = np.array( meanGridpoints )
    plt.plot( meanGridpoints, meanErrorDrag, "o")
    plt.plot( meanGridpoints, 2*10**1*meanGridpoints**(-1/3), label="1st order", linewidth=1, linestyle="--" )
    plt.plot( meanGridpoints, 6.6*10**3*meanGridpoints**(-2/3), label="2nd order", linewidth=1, linestyle="--" )
    plt.plot( meanGridpoints, 2*10**6*meanGridpoints**(-3/3), label="3rd order", linewidth=1, linestyle="--" )
    plt.ylim([0.01,1])
    plt.xlim([4*10**6,10**8])
  elif refinement == "time":
    h = np.array( timesteps[1:] ).astype(np.float)
    plt.plot( h, meanErrorDrag, "o")
    plt.plot( h, 2*10**3*h, label="1st order", linewidth=1, linestyle="--" )
    plt.plot( h, 1.5*10**5*h**(2), label="2nd order", linewidth=1, linestyle="--" )
    # plt.ylim([0.02,5])
    plt.xlabel("Timestep")
  elif refinement == "vorticity":
    plt.plot( h, vorticityDiff, "o")
    plt.plot( h, 3*10**-1*h**(-1/3), label="1st order", linewidth=1, linestyle="--" )
    plt.plot( h, 6*10**1*h**(-2/3), label="2nd order", linewidth=1, linestyle="--" )
    plt.plot( h, 1.2*10**4*h**(-3/3), label="3rd order", linewidth=1, linestyle="--" )
    plt.ylim([0.0005,0.02])
    plt.xlim([10**6,3*10**7])
    plt.xlabel("Number of Gridpoints")

  plt.xlabel("Number of Gridpoints")
  plt.ylabel("Error")
  plt.xscale("log", base=2)
  plt.yscale("log")
  plt.legend()
  plt.grid(b=True, which='minor', color="white", linestyle='-')
  plt.show()

def plotNumBlocks():
  rootSCRATCH = "/scratch/snx3000/pweber/CubismUP3D/"
  cases= [ "1000" ]
  runname = [ "sphereRe1000_levels3_Rtol4" ]

  data = np.loadtxt(rootSCRATCH+runname[0]+"/diagnostics.dat", skiprows=1)

  plt.plot(data[:,1], data[:,-1])
  plt.show()

if __name__ == '__main__':
  plotValidation()
  # plotVortexSheet()
  # plotVerification()
  # plotNumBlocks()
