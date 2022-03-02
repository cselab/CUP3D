#!/bin/bash
NNODE=4
FACTORY='CarlingFish L=0.2 T=1.0 xpos=0.3 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
'
OPTIONS=
OPTIONS+=" -extentx 1.0"
OPTIONS+=" -bpdx 8 -bpdy 4 -bpdz 4"
OPTIONS+=" -tdump 1.0 -tend 20.0 "
OPTIONS+=" -CFL 0.3 -nu 0.00001"
OPTIONS+=" -levelMax 4 -levelStart 3 -Rtol 0.5 -Ctol 0.05"
