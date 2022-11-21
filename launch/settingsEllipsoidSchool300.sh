#!/bin/bash
NNODE=512
BPDX=${BPDX:-32}
BPDY=${BPDY:-16}
BPDZ=${BPDZ:-16}
NU=${NU:-0.000004}

FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.826625385337721 ypos=2.390220494302034 zpos=2.1516843428127133 bFixFrameOfRef=1
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.10987654101314 ypos=2.2192573492996117 zpos=1.9839290098545428  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.910572588760578 ypos=2.6619609394253336 zpos=1.253944500103106  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.537178575923566 ypos=1.8344906415104214 zpos=1.0689277451789008  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.957942924798192 ypos=1.122476365052767 zpos=2.1997295693261742  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.649708279279376 ypos=1.718398770515929 zpos=1.6335251439394554  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8009297749800233 ypos=1.285411907063072 zpos=2.882373465075418  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.4761232545631517 ypos=0.9365910628229492 zpos=2.080407813070982  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8137159602596418 ypos=2.4141846792102495 zpos=2.5491073157808533  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.723191456859159 ypos=1.8515971301941874 zpos=1.4738016107141187  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.422616532910009 ypos=2.9074262013553716 zpos=2.6247435902076512  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.500222267594691 ypos=2.2376765471241433 zpos=0.9547620231129295  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8700637179718558 ypos=1.1253692546474459 zpos=2.3789252773148157  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.327279107954431 ypos=1.7058576759008193 zpos=2.028729554003969  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.280289740313054 ypos=2.3210788239800446 zpos=1.7121959799183548  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.9847973298904065 ypos=2.286520951910949 zpos=1.3531458122157645  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.7554666452086396 ypos=1.9010210808761219 zpos=2.8530745878030572  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.619917958644496 ypos=3.110251160590134 zpos=1.7749483875143766  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.9157912997701536 ypos=2.256349498943436 zpos=2.6139306965373845  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.086634702196664 ypos=2.769014669380851 zpos=1.5920413508770992  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.650237905461229 ypos=1.303152065082993 zpos=1.570239611470268  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.2861791205564685 ypos=1.9397613640823164 zpos=1.25230000525303  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.187127147419799 ypos=2.815251796369413 zpos=2.7222561208605756  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.135270081718224 ypos=2.6599011790750144 zpos=2.0832397141259125  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.5999512274596306 ypos=2.4598472348182696 zpos=3.0234520298624754  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.789085737699066 ypos=2.89168032725832 zpos=1.8040496670343453  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.863896130740082 ypos=1.9380197146065483 zpos=2.459731627908184  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.569486758072301 ypos=2.0191636125039683 zpos=1.9201382289757398  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.255145795285411 ypos=1.4629978835509316 zpos=2.3689226027999557  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.42622087754799 ypos=1.1069896961178451 zpos=2.284164124670256  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.5822721556844175 ypos=1.7325964062074557 zpos=1.8817509336872495  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.493396102422196 ypos=2.1946652592703937 zpos=2.377575573088295  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.244226241622911 ypos=1.714171880498814 zpos=1.0535388144750553  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.25694564333502 ypos=1.6341095430740848 zpos=1.483928037799144  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.675607570269249 ypos=1.3121525891624273 zpos=1.7144653304704065  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.193296333605712 ypos=1.2105735054922317 zpos=2.3421433849690225  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.209793129579657 ypos=1.315739675840303 zpos=1.5274880644800801  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.3900083988328276 ypos=1.2391229659661365 zpos=1.5331637471604598  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.1631381727289507 ypos=2.4756038860328315 zpos=2.7858005562417256  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.7011570701346215 ypos=1.5506392547040488 zpos=1.6549789459797142  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.92575017594742 ypos=2.9986240474595176 zpos=2.640608968159101  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.0586307710712366 ypos=2.032295197692446 zpos=1.6324553249249272  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.216363314393598 ypos=1.6297114516761502 zpos=2.420402681691047  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.576698130342359 ypos=2.111949102059379 zpos=1.6630242565526416  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.1210052706126215 ypos=2.1327270769140347 zpos=0.9631688595430772  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.975210973658013 ypos=2.5341307199574947 zpos=1.5188530666708975  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.544382440047817 ypos=2.460820846027824 zpos=2.2297147031101225  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.010631626916099 ypos=1.342406940814439 zpos=1.710469967202576  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8914116610123144 ypos=2.172862823138429 zpos=2.372166743317259  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.0469383445276845 ypos=1.6765152190370247 zpos=2.9537944768463404  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.0401458229506435 ypos=2.0006972497473208 zpos=1.0599389406271686  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.456033866136989 ypos=2.159256139575216 zpos=1.3985377440336202  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.905979722858654 ypos=1.1823444663202185 zpos=2.680721457595382  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.6292358253507753 ypos=2.803631706278818 zpos=1.3006314807852628  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.644139507046281 ypos=2.3577846117272725 zpos=1.6824517679991096  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.110040444190085 ypos=2.9986491705988247 zpos=2.2418126494627555  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.23448823409265 ypos=2.449411848023566 zpos=2.03162557557613  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.19195246365844 ypos=1.7821933772634115 zpos=2.7864311544347156  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.88447597815053 ypos=2.320482990280406 zpos=2.3247808784512687  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.221163974340506 ypos=2.225843328832797 zpos=1.4674212735585264  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.707500222261185 ypos=2.178325180857592 zpos=3.12927237928102  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.9995910369786523 ypos=1.8731682736932846 zpos=2.2216309739239226  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.9636344830727332 ypos=1.9915954581780437 zpos=1.9769570098899536  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.290556704480799 ypos=2.597199550038673 zpos=2.3738707684045828  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.2215352084974245 ypos=2.210109765109996 zpos=3.061370727727385  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.527567069150037 ypos=1.940015939261044 zpos=2.558117623296189  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.5758142169977547 ypos=1.785932987871667 zpos=1.6541678663669692  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8787767593214943 ypos=1.5727782186047403 zpos=2.7454637527460757  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.121120063427125 ypos=2.08629104540738 zpos=1.9574224733413446  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.343555661786026 ypos=2.576000396962454 zpos=1.4694510484594654  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.754903083840798 ypos=2.4603637627968062 zpos=1.7836158269108715  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.047364620268201 ypos=1.9874567396153513 zpos=1.919720645539796  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.202732062397465 ypos=1.069693056278242 zpos=2.1420209342859837  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.937288557887308 ypos=1.6997226097102172 zpos=1.6633232787470922  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.3784701220571565 ypos=1.2432903657327226 zpos=2.720529242731512  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.250650846650401 ypos=1.8375019419094247 zpos=2.7274482563920075  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.704765183478773 ypos=2.252741667303469 zpos=2.0214517440149065  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.783786215304633 ypos=2.020002555295442 zpos=1.8948852679435204  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.343426188418915 ypos=2.4045620950523063 zpos=2.7830477661713457  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.663087064388479 ypos=2.3129033556501066 zpos=0.8931909105961533  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.170170094702329 ypos=1.0454140144603004 zpos=1.5664786799312802  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.521086665552558 ypos=2.1730554683674472 zpos=2.7303084382642506  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.3315424832882123 ypos=2.1176376520094276 zpos=1.397793833894517  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.450405311892256 ypos=2.9283180748367146 zpos=1.5411499983890389  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.6380913185637835 ypos=1.381529296691605 zpos=2.334329935437805  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.836310658648454 ypos=2.516441979866064 zpos=2.9684158101482367  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.870374368323752 ypos=1.5963286475624503 zpos=2.482202000587651  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.8869801644989685 ypos=1.9897706404972602 zpos=2.9180428311261055  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.395482074724716 ypos=2.163971552820543 zpos=2.4587580965317564  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.0874690594913416 ypos=1.5411742689494305 zpos=1.237108430740463  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.268737676834232 ypos=2.2253004554294233 zpos=2.74802837915772  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.743289527629953 ypos=2.30848505587732 zpos=2.78681699494281  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.1196390198032224 ypos=1.7753813573338952 zpos=2.2052533919293813  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.569884572310116 ypos=2.63709778323092 zpos=2.5431116896247534  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.384800651482703 ypos=2.2428804432569125 zpos=2.230931326300617  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.2472538227569947 ypos=2.565271584601553 zpos=2.0542006766636214  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.547028971667811 ypos=1.5528499361767292 zpos=2.1805638124059397  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.3053330622574344 ypos=2.3128969487930453 zpos=1.999404269637664  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.0257768701274292 ypos=2.613139705273868 zpos=1.5662677044470246  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.862245320318919 ypos=2.254999146919136 zpos=1.7786718620638984  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.4307436528275788 ypos=0.8880782452825473 zpos=1.869244065146307  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.2860580851452452 ypos=2.1925660293423745 zpos=2.204350296012617  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.324588270405455 ypos=1.056507379001324 zpos=1.9561429272781434  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.1281100035589064 ypos=2.284918344031408 zpos=3.0078612445887276  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.1167806592667517 ypos=1.6911765542377004 zpos=1.817100616822448  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8164477987415637 ypos=2.8136730734465414 zpos=2.2576787060610393  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.172146334947643 ypos=1.8917973685232023 zpos=2.1863826752163704  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.6015245029922776 ypos=2.047504549453203 zpos=2.1726292833988046  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8136234699094405 ypos=2.602870949193777 zpos=2.3104874443277748  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.3542274766688633 ypos=2.189452586332749 zpos=2.6581805986098237  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.614078023700281 ypos=1.3055661999253185 zpos=2.0153391428406318  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.695578753149553 ypos=1.3085817967445235 zpos=1.3685486296433058  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.7060751403877945 ypos=1.6723378198645242 zpos=1.3354085871818486  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.881552467489518 ypos=0.8617572638795037 zpos=1.849730777946267  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.0679782254526735 ypos=2.06815754699197 zpos=3.08913139462277  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.710943944466782 ypos=2.7418546151997765 zpos=1.705214513536532  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.873439816088888 ypos=1.889379461888417 zpos=1.9464182859973542  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.201181946865965 ypos=1.3900050944677567 zpos=2.631284244606303  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.022270869933541 ypos=2.9222395977995475 zpos=2.3372297692450945  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.248744280807596 ypos=1.0381155975711867 zpos=2.511054733642725  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.7512589798920906 ypos=2.1828607152971937 zpos=2.7587515245177503  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.686845134615665 ypos=1.6256662844397862 zpos=2.1834195457914998  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.2968207224475723 ypos=2.9910601602719686 zpos=1.8658523734815478  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.426286600011104 ypos=2.1280015611879324 zpos=2.549174093842988  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.306488084559608 ypos=1.530496372577068 zpos=0.97160287701062  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.177876525003819 ypos=1.7357066310851468 zpos=2.952769380828144  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.140997849248572 ypos=1.9921152894723444 zpos=1.033531196815059  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.830401012979681 ypos=1.5240129697245797 zpos=1.9306492788493206  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.5711931917629256 ypos=3.0438900637076474 zpos=1.6431248129789728  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.842042134400636 ypos=1.07032509219919 zpos=2.1526601042906246  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.5285374196671824 ypos=2.627535490667178 zpos=2.493291529297369  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.6766557377518945 ypos=1.8446869531430266 zpos=2.122621890287007  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.666892516217399 ypos=1.4087312101247278 zpos=1.396830333834675  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.784464071941893 ypos=1.0261356478091097 zpos=1.373261682642829  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.25492920615533 ypos=1.9243029467732773 zpos=2.760741778655712  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.541266213744273 ypos=1.462274119241051 zpos=2.5914211991684897  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.299669725663721 ypos=2.596218495381723 zpos=2.2696047122279785  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.3089209314381605 ypos=1.8255183329430977 zpos=1.6836706984816887  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.5971341268652575 ypos=1.6068572769428824 zpos=2.6503566742345517  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8665973805341474 ypos=2.268235433735847 zpos=1.107085999470538  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.861431752625092 ypos=1.7084109271468115 zpos=2.2969220888762796  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.816779340838648 ypos=1.4889871208465117 zpos=1.141806548135904  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.576332808072711 ypos=2.2363012348303317 zpos=1.1729758110372805  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.299647783546181 ypos=2.0941173548139846 zpos=1.1286846561821418  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.710544832312326 ypos=2.2008290074365413 zpos=1.4783467597528146  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.572916604091095 ypos=1.6614954866079927 zpos=1.4095013902732132  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.6446122755524164 ypos=2.66042956786555 zpos=2.557919327531874  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.4105866099813116 ypos=1.310594895920119 zpos=1.7421494910398623  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.163404380774067 ypos=0.8693725229975287 zpos=2.0684556998116106  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.212540697532866 ypos=2.0977341774757163 zpos=1.8135086057572911  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8466174657511596 ypos=2.4876565992126776 zpos=1.4971870289662152  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.81820612423697 ypos=2.8687488856359984 zpos=1.379652154346607  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.33403729570245 ypos=1.513173305476039 zpos=1.3092324761636283  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.3374093275894823 ypos=1.9100284775376644 zpos=2.4068454578746823  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.833525692180634 ypos=1.4775036130444497 zpos=1.9180975112125314  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.5089483713146095 ypos=2.3948887122083335 zpos=1.9960612381154579  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.3154339344376393 ypos=1.4061213975336262 zpos=2.7873371721135682  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.77063541841046 ypos=1.8671228340137165 zpos=2.7670454591893194  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=6.062251270029096 ypos=1.8019012906593275 zpos=2.034593670122009  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8484503308862985 ypos=2.1337635538773174 zpos=1.1816157659770354  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.700986360832679 ypos=2.556624731392625 zpos=0.9936907133139552  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.942208192401235 ypos=1.982513897315447 zpos=1.422193093380211  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.520852453550851 ypos=1.921054297000304 zpos=1.1030632775470313  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.422260744629582 ypos=1.539928101717327 zpos=2.1453720398630827  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.2907355795623023 ypos=2.483763961610259 zpos=1.6673599092417324  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.464607129502138 ypos=2.4887852644626185 zpos=2.1596395820374066  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8889487702983474 ypos=2.0661523246390248 zpos=1.4265795064745037  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.9786992294845045 ypos=2.5157493536312323 zpos=2.240038451647209  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.211074662323714 ypos=2.258757336805813 zpos=2.3764791618730796  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.755185146158804 ypos=2.860487017516198 zpos=2.315681003874283  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.3949213819917023 ypos=1.994077872628971 zpos=1.662114406061852  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=6.018072740280319 ypos=1.7110120717394328 zpos=1.7331285695178615  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.182363627406833 ypos=1.7873443418295343 zpos=2.9310642972402015  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.161690076876647 ypos=2.526457690745352 zpos=1.8907339778365995  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.969857510858098 ypos=3.0174779808786143 zpos=1.9127775272019623  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=6.225736499319481 ypos=2.277532978963602 zpos=2.1041346515877786  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.556788696838127 ypos=1.6839324168134386 zpos=2.296691300659947  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.883223330288485 ypos=2.1312349042817305 zpos=1.739286078540992  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.290940355121009 ypos=1.4658063209009145 zpos=2.1666566390777766  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.127425596521303 ypos=2.5005644336681843 zpos=1.4002189936226839  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.316980563462807 ypos=2.194035114495772 zpos=1.415998088662912  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.792885810384369 ypos=2.4109833212648297 zpos=2.370902927519197  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.048845015251885 ypos=2.324516759841655 zpos=1.7516296940619982  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.861131004008836 ypos=2.841033916767567 zpos=2.0485222114312163  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.181733240513844 ypos=2.1465767292444027 zpos=1.1906875911342025  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.1885999277715618 ypos=1.257046571923508 zpos=1.7289664257660373  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.969087592394466 ypos=1.2952923560102392 zpos=1.0932178783389315  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.324861770625796 ypos=2.478771568500011 zpos=2.5665382046540612  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.765138250759631 ypos=2.4640390949967164 zpos=1.748538221943286  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.922804705969495 ypos=2.615929893736468 zpos=2.2415003429816864  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.674401586276683 ypos=2.03191021616978 zpos=1.9033992365084955  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.464706656716035 ypos=1.6328105563727133 zpos=2.38802066757315  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.020604244836614 ypos=2.661494896082859 zpos=2.269293670882019  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.06252253741642 ypos=1.6392382834408301 zpos=0.9340078929306892  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.344944277589923 ypos=1.8038700278944924 zpos=2.574808130277022  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.325021044467363 ypos=1.2649148261863783 zpos=2.0760464113024204  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.678766982266688 ypos=2.187811715029956 zpos=2.459403010343421  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.892744932617542 ypos=2.135530711892178 zpos=1.3534483567133815  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=6.12180904281635 ypos=2.093900339094304 zpos=2.3988682922518225  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.9724690524143793 ypos=2.97003287272141 zpos=1.3525094699792128  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.670009186041846 ypos=2.7504540672733215 zpos=1.892964631838177  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.6575203367645 ypos=3.0708863188340487 zpos=2.2483707539522193  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.00484184671194 ypos=1.6252244812527374 zpos=2.209662779178264  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8105308193004657 ypos=2.09050631295278 zpos=2.410631695621158  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.739932895447939 ypos=1.5989449030774534 zpos=1.3518498143559263  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.347623650720436 ypos=1.1229714251512086 zpos=1.38147979030594  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.605991637475174 ypos=1.6570333580969998 zpos=1.7805637307081363  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.389254878283272 ypos=1.2822765851129168 zpos=2.291921423893032  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.685905369691525 ypos=2.2463961009072118 zpos=1.5478477017608174  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.90819164011983 ypos=1.3249509309101066 zpos=2.4264545219898  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.3594408639039575 ypos=1.8075885238145406 zpos=1.8599552522153764  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.2772527363290815 ypos=1.619160660029316 zpos=1.9780842124465963  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.226086467595115 ypos=2.7475230484686928 zpos=1.3310155781948405  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.674356570897164 ypos=2.1321387348358667 zpos=2.202263616167982  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.336905072054367 ypos=1.3865036050600748 zpos=2.9338349085256348  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.7627891805248135 ypos=2.077467307538448 zpos=1.6739276132237615  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.763426348508801 ypos=1.5600976453335074 zpos=2.4364164569865214  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8855611475785037 ypos=1.397877547294342 zpos=1.340454968143734  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.5841715384733255 ypos=2.8508590320020835 zpos=2.0620054200972144  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.63781467192016 ypos=1.2191283022485422 zpos=2.1392415598200505  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.383301651891491 ypos=2.062123114936283 zpos=2.836002961182359  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.132477437685713 ypos=2.379494737469335 zpos=1.1605141915233694  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.171837323549092 ypos=1.151006523288098 zpos=1.8775036310410276  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.878767284572982 ypos=1.6467712824905598 zpos=1.973948970814434  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.277839980928463 ypos=1.5912198987483412 zpos=2.205669780828739  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.9653675589257125 ypos=1.096121210811099 zpos=1.4619832953666967  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.682512356725769 ypos=1.7567936391734276 zpos=3.0842122361051962  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=6.09298904211237 ypos=1.5074454521521168 zpos=1.9866868279044445  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.704713589975209 ypos=1.9113492413271842 zpos=1.279717510926599  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.423660515245783 ypos=0.893660774409363 zpos=1.9221073313878385  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.836575398642489 ypos=2.454489706623689 zpos=1.2361572514255605  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.065464161868458 ypos=1.8352171471266598 zpos=1.3924227815163983  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.8957916774452706 ypos=2.567246207930307 zpos=1.8662065130946321  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.049128827742861 ypos=2.5267989547640197 zpos=1.11181044968676  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.025275184707638 ypos=0.9518470528141945 zpos=1.9003713385952228  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=6.176561920707558 ypos=2.1262382891071168 zpos=1.758635815354298  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.6902971183802404 ypos=1.2080366036696737 zpos=1.8828559502625644  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.949304559351503 ypos=1.2334261053361886 zpos=1.5973179802873607  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.72036533245263 ypos=0.972779556745279 zpos=2.2142666819479926  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.48991480249386 ypos=2.9076284761296196 zpos=2.0684896964544777  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.859271649462656 ypos=1.1553511599333204 zpos=2.6250937498319296  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.539117145809641 ypos=2.6808557182232624 zpos=2.4649746402139012  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.866968831579412 ypos=1.9635077571628827 zpos=1.0938679957595059  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.0754307662264693 ypos=1.9896559590929173 zpos=1.5900061108754993  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.3474752684682074 ypos=1.8908234222572164 zpos=2.2915282911209203  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.17822454203136 ypos=2.9388788154459924 zpos=2.0403968827507515  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.358002113708979 ypos=2.3839422819481775 zpos=1.808816466127225  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.103620836148759 ypos=2.9148774773136523 zpos=1.9753964829213881  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.860629619093536 ypos=2.4245920791864375 zpos=1.5157980686412633  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.404884383434384 ypos=1.678219057987996 zpos=1.3182612493598453  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.895752857844418 ypos=1.6550467604944707 zpos=1.3505776542115984  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.77866846495513 ypos=1.0342730359018464 zpos=2.530408287952227  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.795524375830919 ypos=2.613057430958264 zpos=2.7384796953560797  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.799078549409183 ypos=2.319492019054804 zpos=2.3058308393743894  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.212236592566112 ypos=2.662011355312552 zpos=1.2320301560958984  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.246893915120177 ypos=1.5802731025194932 zpos=1.1154073619314522  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.711053071311555 ypos=2.564869017280511 zpos=2.084367391956131  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.637230139015965 ypos=1.5790612043324506 zpos=2.777749217133041  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.107597959398744 ypos=1.9595867258550572 zpos=2.98018409607322  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.7807686705354957 ypos=2.79787986197306 zpos=2.461649884112123  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.7265611935070697 ypos=2.7862397836391097 zpos=1.6523805756020233  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.9107393065641785 ypos=1.598991161137372 zpos=1.585220552617553  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.207906706802245 ypos=1.5652799774595898 zpos=2.545970898796443  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.6885922639920268 ypos=1.8998297667085484 zpos=1.694503278654507  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.446054584495829 ypos=1.6988946031985002 zpos=2.8211912621494823  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.14941102177273 ypos=1.9992196753634648 zpos=1.5916453213821649  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.595721515317553 ypos=2.1097433529747414 zpos=1.297127330943545  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.102300939203735 ypos=1.2777332235601644 zpos=2.449163180584696  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.1422626346632376 ypos=2.7129585936739815 zpos=1.3391928574016152  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.302259950970698 ypos=1.851347657555421 zpos=3.040419359353366  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.242288859513058 ypos=1.3849298099364717 zpos=1.5727538729343533  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.29159105712147 ypos=1.4335092404425844 zpos=1.8437214716665768  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.392943788288364 ypos=1.431030362892777 zpos=2.5167597344839274  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.1310777294950025 ypos=2.4845368786678046 zpos=2.482195778376659  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.486122413747852 ypos=2.6481877874134443 zpos=1.1063023047045564  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.248766064657867 ypos=2.841785909580018 zpos=1.5787433616562825  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.419790425601705 ypos=2.373451321523938 zpos=2.3993723605683885  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.462898200314324 ypos=2.4515600385722673 zpos=1.2568153559177708  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.8613515826289695 ypos=2.279961561491114 zpos=1.796726116955066  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.628265906340495 ypos=1.9030368855694169 zpos=2.2750309929760126  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.444093771848673 ypos=1.9453111517743786 zpos=1.641893133683533  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.559042213532442 ypos=2.761309642631419 zpos=1.3318021929439983  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.6401334980465276 ypos=0.9308802675501757 zpos=2.405014642186297  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.614332018749826 ypos=1.0572876755027352 zpos=1.4062651238191834  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.7965709264264045 ypos=1.1567141715727756 zpos=1.9564421903440006  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.437782718987762 ypos=1.596531990118397 zpos=1.8952891912139045  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.221844520413788 ypos=1.020836434018639 zpos=2.479090052706518  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.334429399176257 ypos=1.1346801993706592 zpos=2.280191716864198  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.372682708450375 ypos=2.0645075843189775 zpos=3.178062480863285  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.3747733863658467 ypos=1.4186736534976947 zpos=2.078697313925324  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.21575982649086 ypos=2.2758104325987767 zpos=2.930288552003408  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.143874237638417 ypos=1.6640369994403863 zpos=1.6330029950394913  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.345430506275063 ypos=1.6342003809780794 zpos=2.32511237633942  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.623368897925456 ypos=2.95886669369072 zpos=2.0594235873410587  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.031013846182351 ypos=3.0490545137547658 zpos=1.658014968381332  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=4.349287570704228 ypos=1.8551386231840958 zpos=1.4784848259116599  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.4585542457560354 ypos=2.2113669207009603 zpos=1.4271452490982037  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=3.230953617460462 ypos=2.380741027115603 zpos=1.4192230789751217  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=5.059189974453234 ypos=2.5064175734036596 zpos=1.9975085293581287  
"
FACTORY+="StefanFish L=0.2 T=1.0 CorrectPosition=1 CorrectPositionZ=1 CorrectRoll=1 xpos=2.6591624178615763 ypos=1.2877449276734139 zpos=1.5292894800762127  
"

OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -extent 8.0 -bpdx 16 -bpdy 8 -bpdz 8"
OPTIONS+=" -tdump 0.1 -tend 100.0"
OPTIONS+=" -CFL 0.3 -nu 0.000008 -lambda 1e10 "
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -levelMax 7 -levelStart 3 -levelMaxVorticity 6 -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -restart 0 -checkpointsteps 1000 "
