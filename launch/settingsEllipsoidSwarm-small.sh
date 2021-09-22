#!/bin/bash
NNODE=256
BPDX=${BPDX:-4}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}
NU=${NU:-0.00001}
BC=${BC:-freespace}

FACTORY=
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.44203111357100755 ypos=0.553707267494575 zpos=0.5150228924008031 bCorrectPosition=true heightProfile=stefan widthProfile=stefan bFixFrameOfRef=1 
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.3831634682290578 ypos=0.5298091854157021 zpos=0.43056922629181105 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.5336218294486462 ypos=0.5877606002632008 zpos=0.47829225652189866 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.7875976084981433 ypos=0.5289913392250897 zpos=0.48747954242131875 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.7151360457123023 ypos=0.4683734992041756 zpos=0.5349900112799744 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.7479936903547912 ypos=0.5753697621285049 zpos=0.44655904358958987 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.9157064375391348 ypos=0.4540178145236272 zpos=0.5074154600413379 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.7305006445318668 ypos=0.5449869126610296 zpos=0.5565664875593824 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.4868532167605708 ypos=0.48157722279993964 zpos=0.40717860796179134 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.282360888743479 ypos=0.512856952308338 zpos=0.5575555216625814 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.3852349783552834 ypos=0.447380486027519 zpos=0.4417914896442429 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.8106370274160579 ypos=0.49844297201454113 zpos=0.4366296965361509 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=0.49370783500755977 ypos=0.49150427869243907 zpos=0.47127784895166175 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 2.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 30.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.3 -use-dlm -10 -nu ${NU}"
OPTIONS+=" -levelMax 7 -levelStart 3 -Rtol 1.0 -Ctol 0.05"
OPTIONS+=" -TimeOrder 2"
OPTIONS+=" -poissonTol 5e-7 -poissonTolRel 0"
OPTIONS+=" -verbose 0"
