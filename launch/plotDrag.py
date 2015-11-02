import os
import math
import numpy as np
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib

fig = figure()
fig.canvas.set_window_title('Drag')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

data = []
x = []
d = []

pRe40 = fig.add_subplot(311)
pRe100 = fig.add_subplot(312)
pRe1000 = fig.add_subplot(313)

for dirName, subDirList, fileList in os.walk(rootDir):
	if "Drag" in dirName and "3010" in dirName:
		for file in fileList:
			if "diagnostics.dat" in file:
				data.append(np.genfromtxt(fname=dirName+'/'+file))
				idx = len(data)-1
				dataset = data[idx]
				x.append(dataset[:,1])
				d.append(dataset[:,6])
				if "Re40" in file:
					pRe40.plot(x[idx],d[idx],label=file)
				#					pRe40.legend()
				if "Re100_" in file:
					pRe100.plot(x[idx],d[idx],label=file)
				#					pRe100.legend()
				if "Re1000" in file:
					pRe1000.plot(x[idx],d[idx],label=file)
#					pRe1000.legend()
"""
handlesRe40, labelsRe40 = pRe40.get_legend_handles_labels()
pRe40.legend(handlesRe40, labelsRe40)
handlesRe100, labelsRe100 = pRe100.get_legend_handles_labels()
pRe100.legend(handlesRe100, labelsRe100)
handlesRe1000, labelsRe1000 = pRe1000.get_legend_handles_labels()
pRe1000.legend(handlesRe1000, labelsRe1000)
"""
pRe40.set_xlim(0,5)
pRe40.set_ylim(0,10)
pRe100.set_xlim(0,15)
pRe100.set_ylim(0,2)
pRe1000.set_xlim(0,5)
pRe1000.set_ylim(0,2)

pRe40.set_title("Re40")
pRe40.set_xlabel('Time')
pRe40.set_ylabel('c_D')

pRe100.set_title("Re100")
pRe100.set_xlabel('Time')
pRe100.set_ylabel('c_D')

pRe1000.set_title("Re1000")
pRe1000.set_xlabel('Time')
pRe1000.set_ylabel('c_D')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()