#!/bin/bash
NNODE=${NNODE:-64}
PSOLVER="iterative"

FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.44203111357100755 ypos=0.55370726749457500 zpos=1.01502289240080310 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan bFixFrameOfRef=1
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.38316346822905780 ypos=0.52980918541570210 zpos=0.93056922629181105 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.53362182944864620 ypos=0.58776060026320080 zpos=0.97829225652189866 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.78759760849814330 ypos=0.52899133922508970 zpos=0.98747954242131875 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.71513604571230230 ypos=0.46837349920417560 zpos=1.03499001127997440 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.74799369035479120 ypos=0.57536976212850490 zpos=0.94655904358958987 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.91570643753913480 ypos=0.45401781452362720 zpos=1.00741546004133790 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.73050064453186680 ypos=0.54498691266102960 zpos=1.05656648755938240 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.48685321676057080 ypos=0.48157722279993964 zpos=0.90717860796179134 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.28236088874347900 ypos=0.51285695230833800 zpos=1.05755552166258140 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.38523497835528340 ypos=0.44738048602751900 zpos=0.94179148964424290 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.81063702741605790 ypos=0.49844297201454113 zpos=0.93662969653615090 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=0.49370783500755977 ypos=0.49150427869243907 zpos=0.97127784895166175 Correct=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -poissonTol 1e-4 -poissonTolRel 1e-1"
OPTIONS+=" -bMeanConstraint 2"
OPTIONS+=" -extent 2.0 -bpdx 4 -bpdy 2 -bpdz 4"
OPTIONS+=" -tdump 0.0 -tend 100.0 -fsave 100"
OPTIONS+=" -CFL 0.7 -nu 0.00001 -verbose 1 -restart 0 -checkpointsteps 200 "
OPTIONS+=" -levelMax 6 -levelStart 4 -Rtol 20.0 -Ctol 1.0"