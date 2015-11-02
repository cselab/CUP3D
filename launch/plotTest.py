#import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.patches import Ellipse
from pylab import figure, show

data = np.genfromtxt(fname='/cluster/scratch_xp/public/cconti/CubismUP/test_fields_diagnostics.dat')


x = data[:,7]
y = data[:,8]
a = data[:,11]

h = 0.0125
n = x.size

ells = [Ellipse(xy=(x[i],y[i]), width=0.05, height=h, angle=a[i]*360/(2*math.pi))
		  for i in range(n)]


fig = figure()

a = fig.add_subplot(111, aspect='equal')

increment = 1
for i in range(1,n,increment):
	e = ells[i]
	a.add_artist(e)
	e.set_clip_box(a.bbox)
	e.set_alpha(float(i)/float(n))
	e.set_facecolor((.1,.45,1.))


a.set_xlim(0, 1)
a.set_ylim(0, 1)

show()