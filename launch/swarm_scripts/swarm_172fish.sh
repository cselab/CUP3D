#!/bin/bash
NNODE=256

BPDX=${BPDX:-8}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}

NU=${NU:-0.00004}
BC=${BC:-freespace}

L=0.2

FACTORY=
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.4 zpos=0.4 bFixToPlanar=0 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.4 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.6000000000000001 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.6000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.4 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.4 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.6000000000000001 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.6000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.4 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.4 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.6000000000000001 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.6000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.4 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.6000000000000001 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.8 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.8 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.0 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.0 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.8 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.8 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.0 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.0 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.8 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=0.8 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.0 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.0 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=0.8 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.0 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.2000000000000002 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.2000000000000002 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.4000000000000001 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.4000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.2000000000000002 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.2000000000000002 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.4000000000000001 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.4000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.2000000000000002 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.2000000000000002 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.4000000000000001 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.4000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.2000000000000002 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.4000000000000001 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.6 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.6 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.6 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.6 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.6 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4 ypos=1.6 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.2 ypos=1.6 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.4 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.4 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.6000000000000001 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.6000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.4 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.4 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.6000000000000001 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.6000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.4 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.4 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.6000000000000001 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.6000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.4 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.6000000000000001 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.8 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.8 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.0 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.0 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.8 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.8 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.0 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.0 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.8 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=0.8 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.0 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.0 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=0.8 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.0 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.2000000000000002 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.2000000000000002 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.4000000000000001 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.4000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.2000000000000002 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.2000000000000002 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.4000000000000001 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.4000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.2000000000000002 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.2000000000000002 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.4000000000000001 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.4000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.2000000000000002 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.4000000000000001 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.6 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.6 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.6 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.6 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.6 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8 ypos=1.6 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.6000000000000001 ypos=1.6 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.4 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.4 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.6000000000000001 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.6000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.4 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.4 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.6000000000000001 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.6000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.4 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.4 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.6000000000000001 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.6000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.4 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.6000000000000001 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.8 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.8 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.0 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.0 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.8 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.8 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.0 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.0 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.8 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=0.8 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.0 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.0 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=0.8 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.0 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.2000000000000002 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.2000000000000002 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.4000000000000001 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.4000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.2000000000000002 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.2000000000000002 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.4000000000000001 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.4000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.2000000000000002 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.2000000000000002 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.4000000000000001 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.4000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.2000000000000002 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.4000000000000001 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.6 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.6 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.6 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.6 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.6 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2000000000000002 ypos=1.6 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0 ypos=1.6 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.4 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.6000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.4 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.6000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.4 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.6000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.4 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.8 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.0 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.8 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.0 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.8 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.0 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=0.8 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.2000000000000002 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.4000000000000001 zpos=0.6000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.2000000000000002 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.4000000000000001 zpos=1.0 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.2000000000000002 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.4000000000000001 zpos=1.4000000000000001 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.2000000000000002 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.6 zpos=0.4 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.6 zpos=0.8 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.6 zpos=1.2000000000000002 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.4000000000000001 ypos=1.6 zpos=1.6 bFixToPlanar=0 heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 2.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 100.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.7 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 5 -levelStart 4 -Rtol 0.1 -Ctol 0.01"
OPTIONS+=" -Advection3rdOrder=true"
