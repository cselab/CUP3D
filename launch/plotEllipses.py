import os
import math
import numpy as np
from matplotlib.patches import Ellipse
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 4})

fig = figure(figsize=(4,3), dpi=280)
fig.canvas.set_window_title('Falling Ellipses')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'
h = [0.0125, 0.01, 0.0125, 0.025]

data = []
x = []
y = []
a = []
n = []
ells = []
p = []

plotIdx = 1
	
for I in [1, 2, 3, 4]:
	p.append(fig.add_subplot(1, 4, plotIdx, aspect='equal'))
	fig.tight_layout()
	bIdx = 0
	for B in [16, 32, 64, 128]:
			
		dirName = 'Fields_bpd'+str(B)+'_ic'+str(I)
		fileName = dirName+'_diagnostics.dat'
		fullName = rootDir+dirName+'/'+fileName
		if os.path.isfile(fullName):
			data.append(np.genfromtxt(fname=fullName))
			idx = len(data)-1
			dataset = data[idx]
			x.append(dataset[:,7])
			y.append(dataset[:,8])
			a.append(dataset[:,11])
			n.append(x[idx].size)
			ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.05, height=h[I-1], angle=a[idx][i]*360/(2*math.pi))
						 for i in range(n[idx])])
			increment = 25
			for i in range(1,n[idx],increment):
				e = ells[idx][i]
				p[plotIdx-1].add_artist(e)
				e.set_clip_box(p[plotIdx-1].bbox)
				e.set_alpha(float(i)/float(n[idx]))
				e.set_facecolor((float(bIdx)*.33,float(3-bIdx)*.05,float(3-bIdx)*.33))
        #		else:
        #			print fullName
		bIdx = bIdx+1

	p[plotIdx-1].set_xlim(0, 1)
	p[plotIdx-1].set_ylim(0, 1)
	p[plotIdx-1]
	plotIdx = plotIdx+1

fig.subplots_adjust(hspace=0,wspace=0)
plt.setp([axis.get_xticklabels() for axis in fig.axes[:]], visible=False)
plt.setp([axis.get_yticklabels() for axis in fig.axes[:]], visible=False)
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()

