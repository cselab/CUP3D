#!/bin/bash

# SIM SHOULD GIVE Cd ~ 1.13  (5% OFF AT CURRENT RESOLUTION. SHOULD GET BETTER WITH REFINEMENT)

# https://arc.aiaa.org/doi/pdfplus/10.2514/3.6164
# http://www.mne.psu.edu/cimbala/me325web_Spring_2012/Labs/Drag/intro.pdf
# Reference values for comparison
# Re     Cd
# 20 	 2.85
# 50 	 1.58
# 100 	 1.08
# 200    0.78
# 644    0.525
# 1330   0.436
# 2000   0.430

BASENAME=validateSphere_Re100_smallest
NNODE=16
WCLOCK=24:00:00
FFACTORY=factorySphere

OPTIONS=
OPTIONS+=" -bpdx 64 -bpdy 64 -bpdz 64"
OPTIONS+=" -2Ddump 0"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -uinfx 0.1"
OPTIONS+=" -length 0.05"
OPTIONS+=" -lambda 1e6"
OPTIONS+=" -nu 0.00005"
OPTIONS+=" -tdump 1"
OPTIONS+=" -tend 30"
