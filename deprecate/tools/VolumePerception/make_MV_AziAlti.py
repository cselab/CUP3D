# File:   make_MV_AziAlti.py
# Date:   Fri 19 Dec 2014 11:43:35 AM CET
# Author: Fabian Wermelinger
# Tag:    Compute MV based on azimuth, altitude, r and center
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
    print "USAGE: make_MV.py <r> <azimuth> <altitude> <center_x> <center_y> <center_z>"
    sys.exit(1)

r    = float(sys.argv[1])
azi  = float(sys.argv[2])/180.0 * math.pi
alti = float(sys.argv[3])/180.0 * math.pi
cx   = float(sys.argv[4])
cy   = float(sys.argv[5])
cz   = float(sys.argv[6])

# compute eye
ebx = r * math.cos(alti) * math.cos(azi)
eby = r * math.sin(alti)
ebz = r * math.cos(alti) * math.sin(azi)
Eb = (ebx, eby, ebz)
E  = (ebx+cx, eby+cy, ebz+cz)

# compute up vector
xi = unit( cross(Eb, unit(cross((0,1,0), Eb))) )

# compute view space
Z = unit(Eb)
X = unit( cross(xi, Z) )
Y = unit( cross(Z, X) )
P = (-dot(X, E), -dot(Y, E), -dot(Z, E))

# write in column major order
print "-MV=%f,%f,%f,0,%f,%f,%f,0,%f,%f,%f,0,%f,%f,%f,1" % (X[0],Y[0],Z[0],X[1],Y[1],Z[1],X[2],Y[2],Z[2],P[0],P[1],P[2])
