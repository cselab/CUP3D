import numpy as np
import scipy
import matplotlib.pyplot as plt
import argparse
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import glob
import pickle
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def readData (root, fish):
  v1 = []
  for i in range(fish):
      vi =           pd.read_csv(root+   "/velocity_"+str(i)+".dat",delimiter=' ',usecols=['time','CMx','CMy','CMz','vel_x','vel_y','vel_z'])
      vi = {**vi ,**(pd.read_csv(root+"/powerValues_"+str(i)+".dat",delimiter=' ',usecols=['defPower','EffPDef']))}
      vi = {**vi ,**(pd.read_csv(root+"/forceValues_"+str(i)+".dat",delimiter=' ',usecols=['drag']))}
      v1.append(vi)
  names = v1[0].keys()

  t = np.asarray(v1[0]["time"])
  v = {}
  for n in names:
      temp = np.zeros((fish,len(t)))
      for j in range(fish):
          temp[j,:] = v1[j][n]
      v[n] = temp
  v["speed"] = np.sqrt(v['vel_x']**2+v['vel_y']**2+v['vel_z']**2)

  cot   = np.zeros((fish,len(t)))
  ryz   = np.zeros((fish,len(t)))
  delta = np.zeros((fish,len(t)))
  xm0 = np.average(np.asarray(v["CMx"][:,0]))
  xmm = np.average(np.asarray(v["CMx"]),axis=0)
  ymm = np.average(np.asarray(v["CMy"]),axis=0)
  zmm = np.average(np.asarray(v["CMz"]),axis=0)
  for i in range (fish):
      ryz  [i,:] = ((v["CMy"][i]-ymm)**2+(v["CMz"][i]-zmm)**2)**0.5
      delta[i,:] = (((v["CMx"][i]-xmm) - (v["CMx"][i][0]-xm0))**2 + (v["CMy"][i]-v["CMy"][i][0])**2 + (v["CMz"][i]-v["CMz"][i][0])**2)**0.5

  u = np.asarray(v['speed'])
  d = np.abs(np.asarray(v['defPower']))
  for j in range (len(t)):
      if (t[j]+1.0 <= t[-1]):
        idx = np.logical_and(t>=t[j],t<=t[j]+1.0)
        cot[:,j] = np.trapz(np.asarray(d[:,idx]),t[idx])/np.trapz(u[:,idx],t[idx])
      else:
        cot[:,j] = cot[:,j-1]

  v["COT"] = cot
  v["defPower"] = -v["defPower"]
  v["ryz"] = ryz
  v["delta"] = delta
  print("\nReading completed.")
  return v

def correlationCoefs(data):
    t = np.asarray(data["time"][0])
    names = data.keys()
    fish  = data["time"].shape[0]
    tidx  = data["time"].shape[1]
    coefs = {}
    for n1 in ["CMx","ryz","delta"]:
        for n2 in ["drag","COT","EffPDef","defPower","speed"]:
            print(n1,n2)
            corr = np.zeros(tidx)
            for j in range(tidx):
                _,_,corr[j],_,_ = stats.linregress(data[n1][:,j],data[n2][:,j])
            coefs[n1+n2]=corr
    return coefs

def correlationCoefs1(data,tmin,tmax,n1,n2,kount):
    t = np.asarray(data["time"][0])
    indices = np.logical_and(t >= tmin,t<=tmax)
    t = t[indices]
    tmin = t[ 0]
    tmax = t[-1]
    names = [n1,n2]
    q = {}
    for n in [n1,n2]:
        e = np.zeros(fish)
        for i in range(fish):
            e[i] = np.trapz(data[n][i][indices],t)/(tmax-tmin)
        q[n]=np.abs(e)

    _,_,r_value,_,_ = stats.linregress(q[n1],q[n2])
    return r_value

#def createDragPlot(data,fish):
#    fig, axs = plt.subplots(1, 1)
#    for idx in range(fish):
#      x = data["time"][idx]
#      #y = data["defPower"][idx]
#      y = data["drag"][idx]
#      dydx = data["CMx"][idx]-np.average(np.asarray(data["CMx"]),axis=0) - (data["CMx"][idx][0]-np.average(np.asarray(data["CMx"][:,0]),axis=0))
#      points = np.array([x, y]).T.reshape(-1, 1, 2)
#      segments = np.concatenate([points[:-1], points[1:]], axis=1)
#      norm = plt.Normalize(dydx.min(), dydx.max())
#      lc = LineCollection(segments, cmap='viridis', norm=norm)
#      lc.set_array(dydx)
#      lc.set_linewidth(2)
#      line = axs.add_collection(lc)
#      if idx == 0:
#        fig.colorbar(line, ax=axs)
#        axs.set_xlim(x.min(), x.max())
#        axs.set_ylim(y.min(), y.max())
#        #axs.set_ylim(-4e-6,4e-6)
#    plt.show()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--root' ,required=True            ,type=str)
  parser.add_argument('--root2',required=False,default="",type=str)
  args = vars(parser.parse_args())
  root  = args['root' ]
  root2 = args['root2']
  fish = len(glob.glob1(root,"powerValues*.dat"))

  data = {}
  try:
      print("Reading data from file...")
      with open(root+'-data.proc', 'rb') as myfile:
           data = pickle.load(myfile)
      print("Reading data from file completed.")
  except:
      data = readData(root,fish)
      print("Writing data to file...")
      with open(root+'-data.proc', 'wb+') as myfile:
           pickle.dump(data, myfile)
      print("Writing data to file completed.")
  '''
  coefs = {} 
  try:
      print("Reading coefs from file...")
      with open(root+'-coefs.proc', 'rb') as myfile:
           coefs = pickle.load(myfile)
      print("Reading coefs completed.")
  except:
      coefs = correlationCoefs (data)
      print("Writing coefs to file...")
      with open(root+'-coefs.proc', 'wb+') as myfile:
           pickle.dump(coefs, myfile)
      print("Writing coefs to file completed.")
  '''
  t  = np.asarray(data["time"][0])
  tm = np.arange(t[0],t[-1],0.1)
  fig, axs = plt.subplots(3, 1)
  fig.suptitle(root, fontsize=16)
  j = 0
  names1 = ["x","$r_a$","$\delta$"]
  names2 = ["$r_{p,Drag}$","$r_{p,CoT}$","$r_{p,\eta}$","$r_{p,P_{def}}$","$r_{p,U}$"]
  i1 = 0
  for n1 in ["CMx","ryz","delta"]:
      i2 = 0
      for n2 in ["drag","COT","EffPDef","defPower","speed"]:
          c1 = []
          kount = 0
          for tmax in tm:
              c1.append(correlationCoefs1(data,tmax,tmax+1.0,n1,n2,kount))
              kount += 1 
          n = names2[i2]
          axs[j].set_title("$p=$"+names1[i1])
          axs[j].plot(tm,c1          ,label=n)
          axs[j].set_ylim(-1,1)
          axs[j].yaxis.set_ticks(np.linspace(-1,1,5))
          axs[j].grid()
          axs[j].set_aspect(2)
          i2 +=1
      i1 +=1
      j += 1
  axs[2].set_xlabel("time")
  axs[2].legend(ncol=5,bbox_to_anchor=(1.1,-0.55),fancybox=True, shadow=True)
  plt.tight_layout()
  plt.savefig(root+"--all.pdf")
  plt.clf()
  if (root2 != ""):
    plt.title(root)
    for i in range(fish):
      plt.plot(t,data["COT"][i],color='blue')
    data1 = readData(root2,1)
    plt.plot(data1["time"][0],data1["COT"][0],color='orange')
    plt.plot(t,np.average(np.asarray(data["COT"]),axis=0),color='red')
    plt.show()
    plt.savefig("cot-"+root+'.pdf')
