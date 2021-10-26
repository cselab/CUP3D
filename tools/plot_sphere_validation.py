import numpy as np
import matplotlib.pyplot as plt

def plotValidation():
  root = "/project/s929/mchatzim/SphereValidation/"
  rootVALID = root + "validation-data/"
  speed = 0.125
  radius = 0.0625
  cases= ["300","1000"] #Re to plot
  runname= ["Re300","Re1000"]

  colors = ["lightblue","maroon"]
  fig, axs = plt.subplots(2)
  msize = 6
  for i in range( len(cases) ):
    if cases[i] == "300":
      ploumhansData = np.loadtxt(rootVALID+"Re300.txt"        , delimiter=",")
      johnsonData   = np.loadtxt(rootVALID+"Re300-Johnson.txt", delimiter=",")
      spietzData    = np.loadtxt(rootVALID+"Re300-Spiez.txt"  , delimiter=",")
      axs[i].plot(ploumhansData[ :,0], ploumhansData[ :,1], marker="2", color=colors[i], linestyle='None', label="Ploumhans et al. (2002)" ,markersize=msize,markeredgecolor="black") # , Winckelmans, Salmon, Leonard, and Warren
      axs[i].plot(johnsonData  [8:,0], johnsonData  [8:,1], marker="^", color=colors[i], linestyle='None', label="Johnson and Patel (1998)",markersize=msize,markeredgecolor="black")
      axs[i].plot(spietzData   [8:,0], spietzData   [8:,1], marker="H", color=colors[i], linestyle='None', label="Spietz et al. (2017)"    ,markersize=msize,markeredgecolor="black") # , Hejlesen, and Walther
    elif cases[i] == "1000":
      ploumhansData = np.loadtxt(rootVALID+"Re1000.txt"       , delimiter=",")
      mimeauData    = np.loadtxt(rootVALID+"Re1000-Mimeau.txt", delimiter=",")
      spietzData    = np.loadtxt(rootVALID+"Re1000-Spietz.txt", delimiter=",")
      axs[i].plot(ploumhansData[:  ,0], ploumhansData[:  ,1], marker="2", color=colors[i], linestyle='None', label="Ploumhans et al. (2002)",markersize=msize,markeredgecolor="black") #, Winckelmans, Salmon, Leonard, and Warren
      axs[i].plot(mimeauData   [:45,0], mimeauData   [:45,1], marker="^", color=colors[i], linestyle='None', label="Mimeau et al. (2016)"   ,markersize=msize,markeredgecolor="black") #, Cottet, and Mortazavi
      axs[i].plot(spietzData   [:  ,0], spietzData   [:  ,1], marker="H", color=colors[i], linestyle='None', label="Spietz et al. (2017)"   ,markersize=msize,markeredgecolor="black") #, Hejlesen, and Walther
    #Present results
    data = np.loadtxt(root+runname[i]+"/forceValues_surface_0.dat", skiprows=1)
    time = data[:,1]*(speed/(2*radius))
    totDrag  = data[:, 3]*(2.0/(radius*radius*np.pi*speed*speed))
    axs[i].plot(time[::10], totDrag[::10], color=colors[i], label="Present method (Re="+cases[i]+")")
    axs[i].set_xlabel("Time")# $T=tu/D$")
    axs[i].set_ylabel("Drag Coefficient")# $C_D=2|F_x|/\pi r^2u^2$")
    axs[i].set_ylim([0,1])
    axs[i].legend(ncol=2,loc="lower right",prop={'size': 8})

  axs[0].set_xlim([6.5,50])
  axs[1].set_xlim([6.5,25])
  plt.tight_layout()
  plt.savefig("sphere-validation.eps")

if __name__ == '__main__':
  plotValidation()
