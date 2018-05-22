# File:   make_MV.py
# Date:   Sun 02 Nov 2014 01:03:18 PM CET
# Author: Fabian Wermelinger
# Tag:    matrix view
#         https://solarianprogrammer.com/2013/05/22/opengl-101-matrices-projection-view-model/
# Copyright (c) 2014 Fabian Wermelinger. All Rights Reserved.
import sys
import math

def cross(x,y):
    a = x[1]*y[2] - x[2]*y[1]
    b = x[2]*y[0] - x[0]*y[2]
    c = x[0]*y[1] - x[1]*y[0]
    return (a,b,c)

def dot(x,y):
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]

def norm(x):
    return math.sqrt(dot(x,x))

def unit(x):
    invIXI = 1.0/norm(x)
    return (x[0]*invIXI, x[1]*invIXI, x[2]*invIXI)

if (len(sys.argv) != 7):
    print "USAGE: make_MV.py center_x center_y center_z eye_x eye_y eye_z"
    sys.exit(1)

UPvector = (0, 1, 0) # camera rotation relative to world space

cx = float(sys.argv[1])
cy = float(sys.argv[2])
cz = float(sys.argv[3])
px = float(sys.argv[4])
py = float(sys.argv[5])
pz = float(sys.argv[6])

Z = unit((px-cx, py-cy, pz-cz))
X = unit( cross(UPvector, Z) )
Y = unit( cross(Z, X) )
P = (-dot(X,(px,py,pz)), -dot(Y,(px,py,pz)), -dot(Z,(px,py,pz)))

# write in column major order
print "-MV=%f,%f,%f,0,%f,%f,%f,0,%f,%f,%f,0,%f,%f,%f,1" % (X[0],Y[0],Z[0],X[1],Y[1],Z[1],X[2],Y[2],Z[2],P[0],P[1],P[2])
