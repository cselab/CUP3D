#!/bin/bash
NNODE=32
BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}
NU=${NU:-0.00004}
BC=${BC:-freespace}


FACTORY=
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.1080277540977597 ypos=0.5343802978933964 zpos=0.6246744912135267 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.7405501354398961 ypos=0.4682364004907246 zpos=0.5121285898299222 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.33069446025092 ypos=0.7447318304702608 zpos=0.5576071182846442 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.115921692189601 ypos=0.5796368928459141 zpos=0.23494788421956636 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.9987434613740864 ypos=0.3253123366614823 zpos=0.4285159914781921 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.2217581282636398 ypos=0.22072911954886365 zpos=0.5727356699443727 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.3209419620743377 ypos=0.35720737433530925 zpos=0.39081910068769266 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.3310579947272236 ypos=0.5928797987866349 zpos=0.4114116969727389 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.5014896077895912 ypos=0.49412037254334146 zpos=0.3264989170971151 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0729046019377853 ypos=0.6441139994660597 zpos=0.4092795898592447 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8382388203974418 ypos=0.7042843030115031 zpos=0.44761876441960974 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.9672721814851465 ypos=0.26503815749546733 zpos=0.6384544763886322 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.0215888253312801 ypos=0.7136342113755899 zpos=0.5918294991790969 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8713534847387574 ypos=0.6032812445540414 zpos=0.7352872624031679 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=1.3978755105542626 ypos=0.4130738613142635 zpos=0.6883847240426346 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8859337955747113 ypos=0.6150827533559149 zpos=0.2948157096762496 heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 2.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 50.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.5 -use-dlm -1 -nu ${NU}"
OPTIONS+=" -levelMax 5 -levelStart 3 -Rtol 0.1 -Ctol 0.01"
OPTIONS+=" -TimeOrder 1"
