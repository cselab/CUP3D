#!/bin/bash
NNODE=365

BPDX=${BPDX:-8}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}

NU=${NU:-0.00001}
BC=${BC:-freespace}

L=0.1

FACTORY=
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=1.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=1.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=2.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=2.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=3.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=3.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=3.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=3.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=3.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=3.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=3.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.25 ypos=3.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.0 ypos=3.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=1.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=1.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=2.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=2.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=3.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=3.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=3.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=3.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=3.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=3.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=3.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.75 ypos=3.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=1.5 ypos=3.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=1.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=1.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=2.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=2.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=3.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=3.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=3.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=3.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=3.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=3.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=3.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.25 ypos=3.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.0 ypos=3.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=1.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=1.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.25 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.25 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.25 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.25 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.25 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.5 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.75 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.5 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.75 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.5 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.75 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.5 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.75 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=2.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=2.75 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=3.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=3.0 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=3.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=3.0 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=3.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=3.0 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=3.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.75 ypos=3.0 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=2.5 ypos=3.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=1.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.25 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.25 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.25 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.25 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.5 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.75 zpos=1.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.5 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.75 zpos=1.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.5 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.75 zpos=2.25 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.5 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.75 zpos=2.75 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=2.5 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=3.0 zpos=1.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=3.0 zpos=1.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=3.0 zpos=2.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=3.0 zpos=2.5 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.1 T=1.0 xpos=3.0 ypos=3.0 zpos=3.0 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+="-extentx 4.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.05 -tend 50.0 "
OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.3 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 8 -levelStart 4 -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -Advection3rdOrder=true"
