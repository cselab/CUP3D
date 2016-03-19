import numpy as np
import math

rootDir = '/scratch/daint/cconti/'
dirName = 'test'
fileName = 'Samara_512_SP_160316_Dora_Euler_DLM1_lowDump_diagnostics.dat'
fullName = rootDir+dirName+'/'+fileName

data = []
step = []
time = []
dt = []
r00 = []
r01 = []
r02 = []
r10 = []
r11 = []
r12 = []
r20 = []
r21 = []
r22 = []
u = []
v = []
w = []
n = []
x = []
y = []
z = []
qw = []
qx = []
qy = []
qz = []

data.append(np.genfromtxt(fname=fullName))
idx = len(data)-1
dataset = data[idx]

step = dataset[:,0]
time = dataset[:,1]
dt = dataset[:,2]
r00 = dataset[:,10]
r01 = dataset[:,11]
r02 = dataset[:,12]
r10 = dataset[:,13]
r11 = dataset[:,14]
r12 = dataset[:,15]
r20 = dataset[:,16]
r21 = dataset[:,17]
r22 = dataset[:,18]
u = dataset[:,19]
v = dataset[:,20]
w = dataset[:,21]

for i in range(len(step)):
	qw.append(np.sqrt(max(0,1+r00[i]+r11[i]+r22[i]))/2.)
	qx.append(np.sqrt(max(0,1+r00[i]-r11[i]-r22[i]))/2.)
	qy.append(np.sqrt(max(0,1-r00[i]+r11[i]-r22[i]))/2.)
	qz.append(np.sqrt(max(0,1-r00[i]-r11[i]+r22[i]))/2.)	

for i in range(len(step)):	
	qx[i] = math.copysign(qx[i],r21[i]-r12[i])
	qy[i] = math.copysign(qy[i],r02[i]-r20[i])
	qz[i] = math.copysign(qz[i],r10[i]-r01[i])

	#print qw[i]

x.append(0)
y.append(0)
z.append(0)
for i in range(1,len(step),1):
	x.append(x[i-1]+dt[i]*u[i])
	y.append(y[i-1]+dt[i]*v[i])
	z.append(z[i-1]+dt[i]*w[i])

	print step[i],time[i],x[i],y[i],z[i],u[i],v[i],w[i],qw[i],qx[i],qy[i],qz[i]
