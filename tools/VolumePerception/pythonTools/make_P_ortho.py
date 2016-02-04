# File:   make_P.py
# Date:   Sun 02 Nov 2014 01:48:42 PM CET
# Author: Fabian Wermelinger
# Tag:    make projection matrix, coordinates are given in view space
#         http://www.songho.ca/opengl/gl_projectionmatrix.html
# Copyright (c) 2014 Fabian Wermelinger. All Rights Reserved.
import sys
import math

if (len(sys.argv) != 5):
    print "USAGE: make_P.py width aspect near far"
    sys.exit(1)

width  = float(sys.argv[1])
aspect = float(sys.argv[2])
n      = float(sys.argv[3])
f      = float(sys.argv[4])
r = width/2.0
l = -r
t = r/aspect
b = -t

print "-P=%f,0,0,0,0,%f,0,0,0,0,%f,0,%f,%f,%f,1" % (2.0/(r-l), 2.0/(t-b), -2.0/(f-n), -(r+l)/(r-l), -(t+b)/(t-b), -(f+n)/(f-n))
