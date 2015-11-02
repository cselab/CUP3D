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

data = []
x = []
y = []
a = []
n = []
ells = []
p = []

plotIdx = 1

for H in ['.0125', '.00625', '.003125']:
	for N in ['0.0000001', '0.000001', '0.00001', '0.0001']:
		p.append(fig.add_subplot(3, 4, plotIdx))
	#, aspect='equal'))
		fig.tight_layout()
		dIdx = 0
		for D in ['1.01', '1.1', '1.5', '2.']:
			dirName = 'Fields_paramTest_bpd64_H'+H+'_Nu'+N+'_rhoS'+D
			fileName = dirName+'_diagnostics.dat'
			fullName = rootDir+dirName+'/'+fileName
			if os.path.isfile(fullName):
				data.append(np.genfromtxt(fname=fullName))
				idx = len(data)-1
				dataset = data[idx]
				x.append(dataset[:,7]-.3+.2*dIdx)
				y.append(dataset[:,8])
				a.append(dataset[:,11])
				n.append(x[idx].size)
				ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.05, height=float(H)*2, angle=a[idx][i]*360/(2*math.pi))
							 for i in range(n[idx])])
				increment = 10
				for i in range(1,n[idx],increment):
					e = ells[idx][i]
					p[plotIdx-1].add_artist(e)
					e.set_clip_box(p[plotIdx-1].bbox)
					e.set_alpha(float(i)/float(n[idx]))
					e.set_facecolor((float(dIdx)*.33,float(3-dIdx)*.05,float(3-dIdx)*.33))
			else:
				print fullName
			dIdx = dIdx+1

		p[plotIdx-1].set_xlim(-.3, 1.3)
		p[plotIdx-1].set_ylim(0, 1)
		p[plotIdx-1]
		plotIdx = plotIdx+1

fig.subplots_adjust(hspace=0,wspace=0)
plt.setp([axis.get_xticklabels() for axis in fig.axes[:]], visible=False)
plt.setp([axis.get_yticklabels() for axis in fig.axes[:]], visible=False)
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()