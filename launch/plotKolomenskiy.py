import os
import math
import numpy as np
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib

fig = figure()
fig.canvas.set_window_title('Falling Cylinders')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

data = []
t = []
v = []
xref = [0, 2]
vrefk = [-.53, -.53]
#vrefk = [-.39, -.39]

pku = fig.add_subplot(111)

for dirName, subDirList, fileList in os.walk(rootDir):
	for file in fileList:
		if "diagnostics.dat" in file and "Kolomenskiy" in file and "1e6" in file:
# and "Avgu" in file:
			data.append(np.genfromtxt(fname=dirName+'/'+file))
			idx = len(data)-1
			dataset = data[idx]
			t.append(dataset[:,1])
			v.append(dataset[:,10])
			pku.plot(t[idx],v[idx],label=file)
			pku.legend()


#handles, labels = pku.get_legend_handles_labels()
#pku.legend(handles, labels)

pku.plot(xref,vrefk)

pku.set_xlim(0,2)
pku.set_ylim(-1,0)

pku.set_title("Kolomenskiy")
pku.set_xlabel('Time')
pku.set_ylabel('velocity v')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()
