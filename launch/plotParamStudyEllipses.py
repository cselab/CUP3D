import os
import math
import numpy as np
from matplotlib.patches import Ellipse
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec

fig = figure(figsize=(16,9))
fig.canvas.set_window_title('Falling Cylinders')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

data = []
u = []
v = []
t = []
dtdt = []
x = []
y = []
a = []
n = []
ells = []

gs = gridspec.GridSpec(2,3)
p  = fig.add_subplot(gs[:,:-1],aspect='equal')
pu = fig.add_subplot(gs[0,2],aspect='equal')
pa = fig.add_subplot(gs[1,2])#,aspect='equal')
colorList = ['#d9001c','#d95100','#0065d9','#0088d9','#5aa600']
counter = 0

for dirName, subDirList, fileList in os.walk(rootDir):
	for file in fileList:
		if "diagnostics.dat" in file and "Andersen_0211_50FPS" in dirName and "_CFL0.1_bpd16" in dirName:
			fname=dirName+'/'+file
			print(fname)
			data.append(np.genfromtxt(fname))
			if "_T_" in file:
				counter=3
			if "_F_" in file:
				counter=1
			
			idx = len(data)-1
			dataset = data[idx]
			u.append(dataset[:,9])
			v.append(dataset[:,10])
			t.append(dataset[:,11])
			dtdt.append(dataset[:,12])
			pu.plot(u[idx],v[idx],label=file,color=colorList[counter],linewidth=1.0)
			pa.plot((t[idx]+math.pi)%(math.pi*2)-math.pi,dtdt[idx],label=file,color=colorList[counter],linewidth=0.0,marker='.',ms=1.)
			
			x.append(dataset[:,7])
			y.append(dataset[:,8])
			a.append(dataset[:,11])
			n.append(x[idx].size)
			ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.05, height=0.00625, angle=a[idx][i]*360/(2*math.pi))
						 for i in range(n[idx])])

			if "CFL0.1" in file:
				increment = 5
			if "CFL0.01" in file:
				increment = 50

			for i in range(1,n[idx],increment):
				e = ells[idx][i]
				p.add_artist(e)
				e.set_clip_box(p.bbox)
				e.set_alpha(float(i)/float(n[idx]))
				e.set_facecolor(colorList[counter])

pu.set_xlim(-.3,.3)
pu.set_ylim(-.4,.2)
pa.set_xlim(-math.pi,math.pi)
pa.set_ylim(-5,5)

pu.set_xlabel('u',fontsize="30")
pu.set_ylabel('v',fontsize="30")
pa.set_xlabel(r'$\theta$',fontsize="30")
pa.set_ylabel(r'$\.\theta$',fontsize="30")

p.set_xlim(0, 1)
p.set_ylim(0, 1)
p.set_xlabel('x',fontsize="30")
p.set_ylabel('y',fontsize="30")

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
fig.savefig("ParamStudyLambda.png",dpi=300,transparent=True)
show()