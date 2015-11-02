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
data2 = []
x = []
y = []
r = []
v = []
xa = []
a = []
aT = []
xref = [0, 20]
rref = [157, 157]
vref = [-0.0443709505194559, -0.0443709505194559]

p1 = fig.add_subplot(221)
p2 = fig.add_subplot(222)
p3 = fig.add_subplot(223)
p4 = fig.add_subplot(224)

for dirName, subDirList, fileList in os.walk(rootDir):
	if "FallingCylinder_Namkoong_2810" in dirName:
		for file in fileList:
			if "diagnostics.dat" in file and "DLM10" in file and "CFL0.1" in file:
				#and "128" in file:
				data.append(np.genfromtxt(fname=dirName+'/'+file))
				idx = len(data)-1
				dataset = data[idx]
				x.append(dataset[:,1])
				r.append(dataset[:,6])
				y.append(dataset[:,8])
				v.append(dataset[:,10])
				p1.plot(x[idx],r[idx],label=file)
				p2.plot(x[idx],v[idx],label=file)
				p3.plot(x[idx],y[idx],label=file)
				p3.legend()
                    
			if "addedmass.dat" in file:
				data2.append(np.genfromtxt(fname=dirName+'/'+file))
				idx = len(data2)-1
				dataset2 = data2[idx]
				xa.append(dataset2[:,0])
				a.append(dataset2[:,2])
				aT.append(dataset2[:,3])
				p4.plot(xa[idx],a[idx],label=file)
				p4.plot(xa[idx],aT[idx],label=file)

p1.plot(xref,rref)
p2.plot(xref,vref)

#handles, labels = p3.get_legend_handles_labels()
#p3.legend(handles, labels)

p1.set_xlim(0,20)
p1.set_ylim(0,200)
p2.set_xlim(0,20)
p2.set_ylim(-.06,0)
p3.set_xlim(0,20)
p3.set_ylim(0,1)
p4.set_xlim(0,20)
p4.set_ylim(-.1,0)

p1.set_title("Namkoong")
p1.set_xlabel('Time')
p1.set_ylabel('Re')

p2.set_title("Namkoong")
p2.set_xlabel('Time')
p2.set_ylabel('velocity v')

p3.set_title("Namkoong")
p3.set_xlabel('Time')
p3.set_ylabel('position y')

p4.set_title("Namkoong")
p4.set_xlabel('Timestep')
p4.set_ylabel('acceleration')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()