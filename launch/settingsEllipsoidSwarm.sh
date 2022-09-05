#!/bin/bash
NNODE=512
BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}
NU=${NU:-0.00001}
PSOLVER="cuda_iterative"


FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.603979659359861 ypos=1.6706107417150406 zpos=1.814174499548682 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan bFixFrameOfRef=1 
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6002862640035835 ypos=1.786393273465106 zpos=1.7794975413301144 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.144476090525486 ypos=2.004416112229979 zpos=1.8657457133401358 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3679997769548873 ypos=2.1032678974006465 zpos=1.4911510243720048 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.620137215110823 ypos=1.6867329489587157 zpos=1.9662275878522404 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7962677588361116 ypos=2.0217818673409225 zpos=2.533352841796526 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6381831935918507 ypos=1.8805659605815281 zpos=2.158903051847623 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7656412870687594 ypos=2.169434682145135 zpos=2.3429462242505443 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6377946318627457 ypos=2.2459046031490253 zpos=2.392954217974509 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.32807246113325 ypos=2.0786877675420166 zpos=1.8427261794816754 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9447866007483356 ypos=2.2088407529633733 zpos=2.3867083875967507 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.784538762006627 ypos=2.1634525122055783 zpos=2.3240236011962216 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.125904712693283 ypos=2.119418285844678 zpos=2.4039235029414914 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3200097835336275 ypos=1.4840279653973685 zpos=1.8285706657148355 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6642041998381654 ypos=2.247118096894978 zpos=2.048320371882238 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.249546913886253 ypos=2.0729991696890453 zpos=2.0796860643401733 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.540116270445905 ypos=1.822894357686476 zpos=1.6028358760890926 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0201724516872024 ypos=1.9220192330405896 zpos=2.4472273756798115 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.903368417683583 ypos=1.8704118001669379 zpos=2.3620840218709738 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.874988072682961 ypos=2.293326119518834 zpos=2.0962626971656757 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9580837209292037 ypos=2.1601365713695233 zpos=1.8466029708319296 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.500452640359144 ypos=1.561816451054336 zpos=2.036891190594705 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5776821574339337 ypos=2.290331551579325 zpos=1.730919928932506 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9729218145171368 ypos=2.331557581088527 zpos=1.6620949706027985 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1982442507261606 ypos=1.4767081545884408 zpos=1.9229177749716282 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3323885420163175 ypos=2.170032800261188 zpos=1.586399576467766 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.45437427315514 ypos=2.19149403799901 zpos=1.590353653401807 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8957393802845743 ypos=1.7893536406763872 zpos=1.8686721442902248 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6199192661210784 ypos=2.252007291449397 zpos=2.1918973113917137 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.667378077554017 ypos=1.597348489789476 zpos=2.013879249604721 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.835499642560351 ypos=1.8115513991766758 zpos=1.6474843628901996 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.823923088658541 ypos=1.8637220713188476 zpos=2.29181738641595 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.12230625165122 ypos=1.940278335770679 zpos=2.1958714590203545 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4719942498850873 ypos=2.2387333105717664 zpos=2.1990412572420652 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8044909625591656 ypos=1.9076250155747796 zpos=1.8724476507603087 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.041601614584345 ypos=2.304326305142255 zpos=2.1313470518356104 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8742207069490826 ypos=2.079749115920102 zpos=1.6877440134368817 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7717461489184028 ypos=2.4257257796938294 zpos=2.251232763933939 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.403515606328694 ypos=2.3774544595693197 zpos=2.000140037443242 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4262062119891827 ypos=1.6721929493735308 zpos=1.6570840793885466 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.535955644208952 ypos=2.3926709153734014 zpos=2.207996264289361 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.91650641108246 ypos=1.9889435574268863 zpos=2.2870277397600076 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6559734183536308 ypos=2.005898616871936 zpos=1.8945885904357336 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1235430980872145 ypos=2.358901615615801 zpos=1.7546342666546335 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7588731161508977 ypos=2.230318082545237 zpos=1.677911239040071 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1460602568316682 ypos=1.958561853838407 zpos=2.2813085671370517 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.192514256057233 ypos=2.274118312184007 zpos=1.82366613735588 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5930872264900335 ypos=2.179587185076904 zpos=2.1580277383858397 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.627914020769219 ypos=2.3225121189167424 zpos=1.776632717201636 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.83314486241288 ypos=2.3041827889950586 zpos=1.7143134116476664 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0338792313296414 ypos=2.0535541131078023 zpos=2.3590020817855404 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.333126512425667 ypos=1.9857392219928118 zpos=1.6129177127873169 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2338883693179494 ypos=1.8271054297082305 zpos=1.7684598588175375 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0518850695517283 ypos=2.2085036364886457 zpos=2.004372324683112 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8017155992083613 ypos=1.766995694412346 zpos=2.0235599598274443 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.805125840647167 ypos=1.7485182159873156 zpos=2.5001220997263793 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8725090788953263 ypos=1.975006087502846 zpos=2.079776630962869 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9124296140080954 ypos=1.8577595515119403 zpos=2.085411060743192 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.460763947672957 ypos=2.293850810095587 zpos=2.24605684717897 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2119966340426385 ypos=2.247405050467795 zpos=2.0926668246430387 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.274654303923397 ypos=2.1448219972356157 zpos=1.9353074511508148 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.741462662323762 ypos=1.7036777719799974 zpos=1.9358529706128034 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.935215179288954 ypos=2.0148981751406563 zpos=1.5103086685722418 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.99168412504101 ypos=2.3353386416795794 zpos=1.813809057732875 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.752562552957797 ypos=1.7865765290471665 zpos=1.8718505533811127 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.605194297931106 ypos=1.7805801052131693 zpos=2.04469742639549 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8508307970095363 ypos=2.0237208781171203 zpos=2.090153475923813 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.321450285464072 ypos=1.8789548514359982 zpos=2.543924570384117 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6424970475513985 ypos=1.6367568915069577 zpos=2.3543828560226228 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2495968925393313 ypos=2.4332208220522547 zpos=2.0972738352698403 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.266141938630236 ypos=1.8173215510162473 zpos=2.060918277884109 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5829210784363554 ypos=2.2469297994172193 zpos=2.125454744019284 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.961004592900915 ypos=2.0772235987008023 zpos=1.7298152909303666 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4715370321663843 ypos=1.7922145830432528 zpos=2.3801366480374946 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.91281794590688 ypos=2.146775532729253 zpos=2.347854315514853 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1530632585705924 ypos=2.341753927520849 zpos=1.985687135811768 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0991073983649158 ypos=1.6697518359842904 zpos=1.7291809570272074 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.591381023455834 ypos=1.6751822327888513 zpos=2.3297478044970474 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0515836283302766 ypos=2.4855866850237507 zpos=1.9520947608043329 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.198175814306102 ypos=1.7999353696604945 zpos=2.215979298523834 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3618533797399515 ypos=1.6439867813047544 zpos=1.814999395466307 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9372001174331093 ypos=2.2502709087815678 zpos=2.0617501628498363 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.978704929443134 ypos=1.8181953822804726 zpos=1.9245464780007626 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.573618159517034 ypos=1.865223058306612 zpos=2.3268975652142467 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.022024382116765 ypos=1.7781337470693004 zpos=1.5538861093537202 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.109778033114075 ypos=2.0181106857297446 zpos=2.0896546669249454 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7860178126631268 ypos=2.256364845735358 zpos=1.7796349191765815 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.428623023789527 ypos=2.4436462708548206 zpos=2.146828489499692 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5432056648264654 ypos=2.1348330037876235 zpos=2.4535041666637434 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4146014619619827 ypos=1.8783788454092247 zpos=2.406723327267026 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6482488520550507 ypos=2.0956527803386433 zpos=2.4883067037940982 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.868890805853909 ypos=1.9106707850212554 zpos=2.191954091472016 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8742014521343324 ypos=2.0210034048653376 zpos=2.319860615334799 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3636061594582194 ypos=1.6330932323785612 zpos=2.0474821177562372 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6679403226734446 ypos=1.5634486085974222 zpos=2.1825003782821586 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.447409373981588 ypos=1.8330932105881212 zpos=2.0481016486641397 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.574668715617991 ypos=2.214591199159488 zpos=1.7649616174558862 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.569168696925256 ypos=1.9251143814013207 zpos=1.5523441367055455 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.505678964437473 ypos=1.836603158897988 zpos=2.4048270603648634 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3000252982690137 ypos=2.2413699835000878 zpos=2.2977190793890605 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0753966975538 ypos=1.5572129568205768 zpos=1.7913501408842096 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.684853548815897 ypos=2.0319481306940825 zpos=1.5350119959723745 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9149612957341504 ypos=2.5452073351755655 zpos=2.1565713997952938 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9632226930058034 ypos=1.5149395034433166 zpos=2.3002639459917176 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5078030649145817 ypos=1.6204797361370302 zpos=1.9463590631728476 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3027309955000037 ypos=1.432883653063651 zpos=1.9334360142215294 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5953707336040064 ypos=1.5166333589969776 zpos=1.8374223714805307 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8387629376226227 ypos=2.1507320750510246 zpos=2.119354290417249 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6840056490151554 ypos=2.491996666812211 zpos=1.8727251479676625 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.862345330350746 ypos=2.369187982355746 zpos=1.9939209968907954 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1901114590085915 ypos=1.4538549071128655 zpos=1.8460414488737609 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2403385320452305 ypos=2.043344873402622 zpos=2.1656277839105935 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.068498483943145 ypos=2.031713844867048 zpos=2.435890313649928 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6457562941039803 ypos=2.0794833178719228 zpos=1.5675209888588495 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3269191339667885 ypos=1.81062011452989 zpos=2.4005482065026333 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4652351569843733 ypos=1.7377552078815857 zpos=1.9048490134415883 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.989817005748063 ypos=1.626393339965571 zpos=1.6854171607977126 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.456231108273125 ypos=1.9004767130061095 zpos=1.588936593363784 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.95040169740054 ypos=2.3548265081081983 zpos=1.5519002312584593 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0724690026308568 ypos=2.2770864581049537 zpos=1.762534987324013 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2318538379704207 ypos=1.7499347904468407 zpos=2.1490688496539603 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5227072130734007 ypos=1.8921516855817602 zpos=2.103691385190239 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7519635636888893 ypos=2.398775453918438 zpos=1.8040799673885433 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.72807297525419 ypos=1.703785518915378 zpos=2.0783634844830017 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.995473186256607 ypos=2.2848322488799737 zpos=2.32701401489934 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.526811909274684 ypos=1.7609978635392456 zpos=2.1608257076112665 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4379784843398573 ypos=1.910484976546073 zpos=1.8263774307494314 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1913434606591076 ypos=1.52231259840208 zpos=1.8670673778060651 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.575269976630144 ypos=1.7239918727421686 zpos=2.375960642305814 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6923746377952393 ypos=1.5861486766601018 zpos=1.8367391323070419 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.557264663198907 ypos=2.4421637062114705 zpos=1.895127753241104 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.701775501765404 ypos=1.9926108660469883 zpos=2.123080529169905 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6273960610175906 ypos=1.7718696376262109 zpos=1.9537070782633132 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.247460573952476 ypos=1.5775264286839927 zpos=1.7897775645453164 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3816559084502966 ypos=1.704792135477143 zpos=2.1978976072791787 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6461314312545023 ypos=1.9398417697262187 zpos=2.205403377258981 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3177930387648824 ypos=2.312025355345101 zpos=1.9414219591236486 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3536546325463354 ypos=1.662524565192367 zpos=1.8551801220730268 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8234385608214287 ypos=1.6580124465346517 zpos=1.6256711561395996 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0380354594857004 ypos=1.977133498220015 zpos=1.8112078483775327 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0030100660622074 ypos=1.8642205074516152 zpos=1.8004046541186198 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7321154401189314 ypos=1.8118907047792254 zpos=1.9230006872945073 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0186837207194053 ypos=1.807581958505752 zpos=1.679166545141544 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1469542289978945 ypos=1.4470763774393376 zpos=2.214441390856593 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9326261133325264 ypos=2.067611148935324 zpos=1.8307963996921135 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.278302516234581 ypos=1.6741219245830115 zpos=1.9788343345512405 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.577924250725218 ypos=1.9887788385684522 zpos=2.398886649269307 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.570234309285346 ypos=1.821026404336715 zpos=2.1162446937572823 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1400550500026045 ypos=1.976864558373098 zpos=1.9651277047270914 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.343083845230457 ypos=1.615629939640998 zpos=2.1950018477770845 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0282634071586836 ypos=1.7069048956900166 zpos=2.084812135858525 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2460534911432197 ypos=1.9638458977342887 zpos=1.4309513669367013 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7138416285422022 ypos=2.2635999948086756 zpos=1.826970765046047 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4828613782449325 ypos=2.1588525179880493 zpos=1.691432053290799 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.363189609411299 ypos=2.01942046647435 zpos=2.1514024631754696 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7569104792495547 ypos=2.078604430708655 zpos=1.6466254772964417 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2695000835418537 ypos=1.8975386424638236 zpos=2.3433021231668154 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.453072990866635 ypos=2.1639555214620434 zpos=2.283482547255086 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.177852863358917 ypos=1.7722103970859897 zpos=2.5304628183733042 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8787103129113407 ypos=1.6770893538929164 zpos=1.5729940646641871 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.388547766477169 ypos=1.9311847646184779 zpos=2.5086040983229756 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2802871294999103 ypos=2.040569605030244 zpos=1.6896936936548534 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3436104306564407 ypos=2.015671780282148 zpos=2.507133357346846 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0848056115568716 ypos=1.9415307109907651 zpos=1.4033545964138385 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8819256167564014 ypos=2.1238177775422264 zpos=1.4500771553516532 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5232627784920276 ypos=1.9768100527446906 zpos=2.1741296728805186 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.719201335866507 ypos=1.8099355451537595 zpos=1.7177229916135914 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5811070936648814 ypos=1.9869477431089135 zpos=1.9337636682374244 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4550911141032787 ypos=1.7213020621307433 zpos=2.239111447086235 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.47773077168303 ypos=2.3114728103653817 zpos=2.1067328772021683 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9100119864104395 ypos=1.9353711420139723 zpos=2.0554432750263163 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5575253106038853 ypos=1.982113238508278 zpos=2.464721863314112 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.504139831284035 ypos=2.0492731120205154 zpos=2.2243863650435514 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.871811373875312 ypos=1.8646021110568625 zpos=1.4694444663939197 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.080349972949727 ypos=2.0066774678677253 zpos=2.022206632653204 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.511368854917211 ypos=2.0551483744392733 zpos=1.8588111085090218 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0793834480534903 ypos=2.5117751746799533 zpos=1.8013665674495383 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1413524713661882 ypos=1.9607023401654806 zpos=1.782484805711214 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.810061410129931 ypos=2.028456032860857 zpos=2.2676854645246087 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7028425591898797 ypos=2.419662790977525 zpos=2.084385356966528 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.729472281523705 ypos=2.3632428820927323 zpos=2.2933377565519617 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.781051399611906 ypos=1.6627950293979505 zpos=1.7792698348659515 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9849073562816866 ypos=1.7119934293994226 zpos=2.277465727585359 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2418803615259586 ypos=1.775951812085201 zpos=2.0111250827496727 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.002025294187382 ypos=1.6563253679704435 zpos=2.0886741803577054 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7652890462432063 ypos=1.595591968460302 zpos=1.8906380770154159 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.546772610414455 ypos=2.331186914016555 zpos=1.6398125082474933 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2330410796335856 ypos=1.6096904406986599 zpos=1.685271607510452 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6317743443569324 ypos=1.8343403983901345 zpos=1.7913848312251226 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.761960780404861 ypos=1.491825674071735 zpos=2.019373605824588 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2371253697171447 ypos=2.464601639832649 zpos=1.8645777126175331 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.782726270461738 ypos=1.9366659554925774 zpos=2.0274060550674857 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6609583852634002 ypos=2.3272802593179605 zpos=1.8609159132915898 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.658035811293833 ypos=2.3073776620738213 zpos=2.2753240967122688 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.802885338937078 ypos=2.3486634955091747 zpos=1.6780455714522728 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4085732456133555 ypos=2.038332665661879 zpos=2.273224356882515 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.052166590014998 ypos=1.993952160417249 zpos=2.161295008554018 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2142587530336213 ypos=1.865981029824316 zpos=2.424955852410192 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3182901752472422 ypos=2.2661006268725252 zpos=2.0735390188820335 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.942111451103644 ypos=1.5532240016331296 zpos=2.2364361366297496 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5312658001939523 ypos=1.6219860339257943 zpos=1.6625037416319353 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2536455618442193 ypos=2.454426235739236 zpos=2.029589266096258 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.749468385405877 ypos=2.27899555484792 zpos=1.971951184757882 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.687092849929532 ypos=1.8840900166064447 zpos=2.4533760632538866 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3345638126362016 ypos=2.444141279437927 zpos=1.7925759230750293 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.638840549829944 ypos=2.1560450408552327 zpos=2.2894783361677633 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.964970152462628 ypos=2.062354140654507 zpos=1.921577234478 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4184307345477976 ypos=1.857789001237218 zpos=2.339856342454379 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.571541658168017 ypos=2.0324887535703495 zpos=2.051910998850682 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9977737835628524 ypos=1.8480206088297495 zpos=1.5143359748841452 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5612761913248567 ypos=2.076683477489744 zpos=2.492065362088076 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8287863419352526 ypos=1.6230169963024967 zpos=1.6746568268780315 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6700416619796883 ypos=1.9027197051451956 zpos=1.7000020723089633 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5732980943187393 ypos=2.382730645838801 zpos=1.8568642045280979 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7118743123708318 ypos=2.0748684115013987 zpos=2.1163589549426276 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.185666545389431 ypos=1.803520433597373 zpos=2.098761487721649 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.92513735218642 ypos=1.7457917908263814 zpos=1.823702326866011 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.789034327202528 ypos=1.7867117965628556 zpos=2.2880991549810448 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5356965907322024 ypos=2.2239167407078506 zpos=1.9179837295764062 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9353702668607085 ypos=1.7015537151946143 zpos=2.009069775720321 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.947648195656936 ypos=2.059499951102868 zpos=2.030126135842581 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1461222184568394 ypos=2.443945056478634 zpos=1.6118001589971767 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7557585723022346 ypos=1.955214438055112 zpos=1.8532645570811943 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.975925038388761 ypos=2.0971359530528138 zpos=1.5220799574924935 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5289236093162826 ypos=2.0445848407570457 zpos=2.44651972476142 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.999454432675008 ypos=2.2726838140004375 zpos=1.665642934630747 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5436989836212986 ypos=1.9256482662679888 zpos=1.8703826558103986 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8778007906480525 ypos=2.354077075761952 zpos=2.3333736498270876 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.286284738744367 ypos=1.9760526916207335 zpos=2.2825755500309235 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.729141403596389 ypos=2.516342514924321 zpos=2.1651844904304443 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.635984708268197 ypos=2.1643326213480827 zpos=2.353678348484098 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4251173957302488 ypos=1.9374576834469643 zpos=1.9739561636426388 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0587241494866895 ypos=2.1751875373504763 zpos=1.9351203479983192 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4827730628828992 ypos=2.5005762358635217 zpos=2.222684234793453 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3287701215975143 ypos=2.3551200691762895 zpos=2.0329847792256683 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.746557969249568 ypos=1.9680809223063687 zpos=1.7701789046474767 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0601964769974215 ypos=1.729487795191757 zpos=2.0108926771982394 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4598595428018784 ypos=2.1467597110235466 zpos=1.6351876001610295 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0816619544275508 ypos=2.0860512649856373 zpos=1.901893584056424 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.381431177991221 ypos=1.935179435725983 zpos=2.2648910135308324 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.898549840818433 ypos=2.2129866857210576 zpos=2.00327512611691 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.794850811274682 ypos=1.8841962556064829 zpos=2.068384783837569 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9147002839649114 ypos=2.2786688218676865 zpos=2.162071590069783 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.65730130662959 ypos=1.6881960438350516 zpos=1.6896577557767132 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.445007966859906 ypos=1.8870955936245986 zpos=1.690060987714388 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1332404847426316 ypos=1.6389227990358228 zpos=2.012868386158401 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4178190628101013 ypos=2.0259755353566984 zpos=1.9033481298762414 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3027012503230195 ypos=2.2845635876170887 zpos=1.77999941238443 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.947830046251011 ypos=2.0046154468720196 zpos=2.209272720227334 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.612258725340929 ypos=2.2565806650951807 zpos=2.441821047195037 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0583611152131036 ypos=2.16950900634581 zpos=2.4059607089007624 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4975538182616566 ypos=1.8168778879660457 zpos=2.1643167791716493 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1763695354118253 ypos=1.9625297431199562 zpos=1.9973506930157647 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4823635851751047 ypos=2.3636336229164665 zpos=1.709000780187062 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.937044179907859 ypos=2.23839845383406 zpos=2.4704151884192047 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4257839947286057 ypos=2.2164923493735427 zpos=2.1911729168866705 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6628987840343417 ypos=1.6534587138516668 zpos=1.9138327922579974 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.264198943677058 ypos=1.8831042353547112 zpos=2.2709234525743325 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.32324920035156 ypos=2.2933695862172487 zpos=2.0149399623758417 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.284665057881159 ypos=2.1370285191620866 zpos=2.1336695915186894 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.704336851733735 ypos=1.8841798802390786 zpos=2.2941885060338 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1605298843893097 ypos=2.2890124703787693 zpos=2.1736719732900287 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.785256734873514 ypos=1.9399597670951416 zpos=2.101891595361183 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4128847930697397 ypos=2.29517304699537 zpos=1.9991797086119554 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5374532167762176 ypos=1.6555876688660756 zpos=2.0337560933114247 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.596146115448355 ypos=1.9963733962510548 zpos=2.0967562657765098 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.389931447215756 ypos=2.4021843036785056 zpos=2.1713428882578674 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2006762688127135 ypos=2.539699145322418 zpos=1.7603731887457734 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.906340328457982 ypos=2.061348670552664 zpos=1.4453041526614534 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.815316766699014 ypos=2.1689455359484118 zpos=1.9639109193434872 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.910327620419786 ypos=1.7973201190590844 zpos=1.967873452728302 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.816680018655689 ypos=2.4587245408855027 zpos=1.8200333951490246 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5483622125297716 ypos=1.5718527777788633 zpos=1.8283886514774574 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7354914418828433 ypos=1.6991216822319721 zpos=2.039417266814222 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.645993509465907 ypos=2.132777244905151 zpos=2.2143128861223595 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6629523937215267 ypos=2.3415793825815343 zpos=1.6564771494416148 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9580656010976045 ypos=1.8669337024423531 zpos=2.5525560644918377 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5589256271451655 ypos=1.5921636622140147 zpos=2.2142119446780932 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8075538506281292 ypos=2.1056839889699552 zpos=2.4499897389411416 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.805114454851673 ypos=2.256628917330278 zpos=2.312216004043854 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8180213175820694 ypos=2.2790289276736373 zpos=2.303721326035744 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2032494157954288 ypos=1.7505290556928708 zpos=1.731135866575603 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.758288435884486 ypos=2.413936896093176 zpos=1.9548793032166691 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.115707153344611 ypos=2.212777928419044 zpos=2.5405427548516704 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.185114636300431 ypos=2.1363284604364092 zpos=2.0793957584501204 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3947491232302682 ypos=2.071163630986688 zpos=2.445278879020255 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.682688185054429 ypos=1.8350882391054255 zpos=1.8634448272954072 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2425699526478264 ypos=2.1169176264527376 zpos=1.9886001918082965 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.174791429618321 ypos=2.1710660348380877 zpos=1.4541847625354536 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5436317250019442 ypos=2.3971429553458306 zpos=2.360874202872735 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.774506426684647 ypos=1.8035123329223834 zpos=1.7451359347962283 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.468945681971481 ypos=1.784144150859404 zpos=2.4388842442364767 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.870549816271378 ypos=2.074379002882372 zpos=1.580745429705986 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7093777634772453 ypos=1.9708265714144273 zpos=2.3421735578713068 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4160949283266353 ypos=1.6616520632368315 zpos=2.1302523955329598 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6052159624913283 ypos=1.7038763417016594 zpos=2.3739860105319073 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3401392211086978 ypos=1.5067256367994877 zpos=1.9828865704999288 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2198571705348917 ypos=1.786008933176965 zpos=1.9020763608046443 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2869133427690707 ypos=2.1168860458785637 zpos=1.7663897670250988 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4980260621719026 ypos=2.0972323603748193 zpos=1.5977743690557857 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6891646609782507 ypos=1.967829171292102 zpos=1.937482389304612 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.773429278347769 ypos=2.4861248435337435 zpos=1.7764533153286923 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3674933228375994 ypos=2.477662934749087 zpos=2.244545101905839 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6761338851303136 ypos=1.731514161608402 zpos=1.8425988998387677 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3071807106986872 ypos=1.8720379638933866 zpos=2.1525761164730484 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7232879041324245 ypos=2.426766206968883 zpos=1.8779295235207465 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1319479679123075 ypos=2.0828973423590402 zpos=2.091035495338387 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5403514837067913 ypos=1.87466116002252 zpos=2.2845216593642714 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.834767721869182 ypos=2.0838635191832986 zpos=2.0531839955643716 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.618676650942933 ypos=2.242790853090058 zpos=2.436747194961091 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.574806310497988 ypos=2.2379452391047843 zpos=2.1100077128920023 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0627382210946608 ypos=1.8795423794364896 zpos=2.1046066002807673 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.523650590358928 ypos=2.3146334524898458 zpos=1.980791962309997 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.729721629019696 ypos=2.517161265755226 zpos=2.0321579349791823 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5751456129487607 ypos=1.6752485219136408 zpos=2.219874514884328 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.495904866428039 ypos=2.2215823615863117 zpos=1.643825118482618 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.496511192217962 ypos=1.537625775934109 zpos=1.743411009404624 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8511799157354742 ypos=2.0899952721755897 zpos=2.0348518151612813 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4586500165327987 ypos=1.5956793776870497 zpos=1.8950425870024588 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9109928940779413 ypos=2.2037241007647936 zpos=1.4779281531151076 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.474864487449702 ypos=1.5083464878282686 zpos=1.8464188779820458 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6985787053577046 ypos=1.8309602188920115 zpos=2.337565811049097 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2416326449288797 ypos=2.204514449381621 zpos=2.0808742287720183 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.853204863906901 ypos=2.289166779702446 zpos=1.8549466600028413 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.524842406107695 ypos=2.2836110720973397 zpos=2.3365372996550526 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9060615187382632 ypos=1.9469080852919922 zpos=1.9809958727908477 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.087685397481004 ypos=1.5742256989905108 zpos=1.9188506387562414 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3354258485013824 ypos=2.5633983282522044 zpos=1.8855592298067045 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1198033984778926 ypos=2.0940530348848045 zpos=2.2990363400387577 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.951002229877209 ypos=1.5644460992836624 zpos=1.7925861532192437 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0807054128624864 ypos=1.4461458361469903 zpos=2.0703352070958654 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5895506868713993 ypos=1.8916749911853752 zpos=2.239662016392667 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4470473089129805 ypos=1.92568460493752 zpos=2.166000547828523 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4976687270117224 ypos=1.7758638294895357 zpos=1.8692283806838368 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.613421486495739 ypos=2.0596141516154947 zpos=2.243101739197878 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9325342754257533 ypos=2.4756044744438976 zpos=2.0476865240540554 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0001560311913025 ypos=2.478264041968354 zpos=1.6872634173189256 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3923962344080856 ypos=2.170048995938494 zpos=1.7551643465005777 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1151335406695173 ypos=2.1098430216559176 zpos=2.131726392736395 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.419296096179701 ypos=1.998141647672068 zpos=1.9498461274898666 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8474568320419453 ypos=2.0882558014561345 zpos=1.6217417234539995 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8772606259455897 ypos=2.2009181456350224 zpos=2.2882923332842404 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.919503062958913 ypos=2.3134906969098776 zpos=2.3962066783956772 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.688230855597616 ypos=2.2289217327321893 zpos=1.9622107501062847 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5122938716543506 ypos=2.1964732798002697 zpos=1.7380931158139536 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.878196949578736 ypos=1.8134280322310685 zpos=2.4507207020493675 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.770537063025514 ypos=2.437416861899557 zpos=2.1392832161528377 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.571519236096497 ypos=2.276762100608991 zpos=1.6349074045874934 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5772719528329597 ypos=2.004470157568356 zpos=2.1919315201205776 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5787168252399613 ypos=1.9871013371390989 zpos=1.7697103544069448 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.096789562809425 ypos=2.4311865863220405 zpos=2.2712558330985075 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2253683839196166 ypos=2.4210698335773535 zpos=2.3155168889593893 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7968856571218264 ypos=1.9111808569995308 zpos=1.9542220950973066 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.401460461802147 ypos=1.9126626071456432 zpos=1.4626693502629848 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.852581314734015 ypos=1.9296058559029103 zpos=2.4481743968114626 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7797696743135707 ypos=1.7616806850440174 zpos=1.9419529070910018 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4177368572827396 ypos=2.337388709737172 zpos=2.0715037386798594 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7337066362648756 ypos=1.9692565723792959 zpos=2.0623444140129914 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2144615470883857 ypos=1.9828616920297801 zpos=1.719981446459016 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8639714473232973 ypos=2.059362366073499 zpos=2.360064954111132 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5498430659054283 ypos=2.44448540604399 zpos=2.0420238012178795 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6905105675152043 ypos=2.2041331214834545 zpos=1.5793073991583173 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6927327399004133 ypos=2.2604228960073085 zpos=2.354763362299863 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.440477101672508 ypos=2.0224016377389917 zpos=1.684452456205931 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.75161373280804 ypos=1.8243957800227053 zpos=1.660761252651946 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6659953265979364 ypos=2.145883726327188 zpos=1.54857323526449 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4076862168368125 ypos=2.190538143159624 zpos=2.4665923788317894 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9529630737730126 ypos=1.7638918810660913 zpos=2.498407769934807 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.886426343813653 ypos=1.9259349399610652 zpos=1.7763219763324387 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.520454207327296 ypos=2.1792780453650993 zpos=2.48123203373527 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6010126246495635 ypos=2.206491464743417 zpos=2.2832453057751163 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4465946069885747 ypos=1.9782657797115102 zpos=1.7113368158596938 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5214407469931954 ypos=1.8692817304522087 zpos=1.8301101351283873 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9499550928422864 ypos=2.031509012716368 zpos=2.084114035534638 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.386978952905903 ypos=1.7249539437506427 zpos=2.153801798803329 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.515545494341637 ypos=1.609943638297 zpos=1.7198971182890403 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2933340832355147 ypos=2.510626021433931 zpos=2.164525333639968 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.039521025886921 ypos=1.853980491514862 zpos=1.9871651129150318 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.803203543849928 ypos=1.5625533815385684 zpos=1.6867499863332724 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9235688564909 ypos=1.6240247330375572 zpos=2.141826652440627 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2337436964810204 ypos=2.052843815481473 zpos=2.5100572607161773 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7966713588969783 ypos=2.242627140946209 zpos=2.217855330020393 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.800711731189168 ypos=2.3216688952253812 zpos=2.280983451012459 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.537264285066502 ypos=2.02217794179982 zpos=1.642321686457143 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1931193799004602 ypos=1.7083910185975748 zpos=2.4007405001088737 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9513020299676183 ypos=1.8721824407101963 zpos=1.8511817104519923 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0536562425641653 ypos=2.145810348457213 zpos=1.5686404553354292 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.615944161361161 ypos=2.1971441104024736 zpos=1.8499511515155267 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.413268750268642 ypos=1.7013562007584762 zpos=2.4460392427855044 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4331010229054595 ypos=2.351220012883897 zpos=1.5718092382642999 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.249869107321645 ypos=2.4330602969182245 zpos=1.9733478397843918 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3330904287516434 ypos=2.380120551902928 zpos=1.8620536182172125 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.557692929092426 ypos=1.7876718305292394 zpos=1.73743735582288 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.347939658039227 ypos=1.5627914360123873 zpos=2.1655858638674736 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.912453465899682 ypos=2.1897209259507067 zpos=1.8120418020924984 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1411070931707985 ypos=2.014582370709871 zpos=2.2328849494807517 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6695987770708314 ypos=1.6134006795561173 zpos=1.75886984622948 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.349667674583877 ypos=1.565211076860635 zpos=1.8346528436780527 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1729932055324275 ypos=1.8985903506652118 zpos=1.9831759153445256 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8312709131375593 ypos=1.988685612985064 zpos=1.9226818937572947 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4749850167525698 ypos=1.6812635918190413 zpos=2.401402148978715 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5046856855536475 ypos=2.051165714233103 zpos=2.003906079076129 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.207122868642169 ypos=1.6503079380032413 zpos=1.6488116187686246 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.037081660641918 ypos=2.3626222856278374 zpos=2.468089908326238 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8366557773833523 ypos=2.24840750395168 zpos=2.231768469423581 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.319020585547797 ypos=1.8475841441978118 zpos=1.9097373623878697 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.887547869121095 ypos=1.9606239460410293 zpos=1.622506019019643 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.856316383207799 ypos=1.8510496662532114 zpos=1.6712492966615846 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8019157025229586 ypos=2.283734476614419 zpos=2.371705364966454 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3726492766206206 ypos=1.9036409148427482 zpos=1.8076133968372483 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.488873969728596 ypos=2.050958058622079 zpos=1.6312863500433306 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9515688385139423 ypos=1.6723212244166614 zpos=2.0410205456215262 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3343676551473194 ypos=2.1125366809667523 zpos=2.3432498154016814 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7397010007051175 ypos=1.659432293567418 zpos=1.8315547093887785 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2570068426607293 ypos=2.157855024113334 zpos=2.5665766144499678 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.940894164671667 ypos=1.8061987775230954 zpos=1.9195199397632614 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2754511391132217 ypos=2.344614508840441 zpos=1.813962886001845 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3050828286474125 ypos=2.290115412253204 zpos=1.6411176928961124 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5276713497321524 ypos=2.001462057165038 zpos=1.8633359789362056 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.109860662243 ypos=1.471361219212437 zpos=1.7998795553231401 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1229489996500077 ypos=1.8757045393424314 zpos=2.149955269484944 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2940736857255035 ypos=1.7065930440951698 zpos=2.2881555691521385 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.954530360240802 ypos=1.7466984685342017 zpos=2.157770681019109 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2255630515097975 ypos=2.4301690129327804 zpos=1.7314400892267892 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0406437182632255 ypos=2.4517645967602215 zpos=2.0473594106076876 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8379207599857943 ypos=2.2458772734098194 zpos=2.1646064474569267 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3876964180710454 ypos=1.9462371063346258 zpos=2.086825516053725 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.694994980974574 ypos=1.7228020962265747 zpos=1.8150837141963594 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9069690535758332 ypos=2.1430204250857305 zpos=2.058168673657419 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5969496193321824 ypos=2.1268536625950762 zpos=2.032510780784302 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.341336413021827 ypos=2.0754682398190636 zpos=2.1741723541815463 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0456972514156035 ypos=1.99923484645897 zpos=1.9192960492067446 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9582958525194734 ypos=2.2493332044049392 zpos=2.0738769015679757 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8517407537152044 ypos=1.9342065566682083 zpos=1.5462478067465422 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.382882520426454 ypos=2.1979001988958555 zpos=1.6871015167802939 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.476926527558424 ypos=2.2682247215987843 zpos=1.9221476597064953 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6833057694659193 ypos=2.2269826813466103 zpos=2.0051361749656884 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2201418253264773 ypos=2.4937861466992284 zpos=2.3096255236963557 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8030351563895204 ypos=1.8361045720184461 zpos=1.8191978138016924 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4528280215277856 ypos=1.9852709900941066 zpos=1.896714685603454 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.886818994188437 ypos=2.4559150768434166 zpos=1.7045704049315598 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2458627878938247 ypos=2.3012421018196973 zpos=2.2469913925775504 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.610265135424111 ypos=1.9378066584057572 zpos=1.9448856422690675 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.608205827014801 ypos=2.4697906800665237 zpos=2.102471897358915 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4930041735650215 ypos=2.3327744159725494 zpos=2.346865268576551 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1786684243486265 ypos=1.5315512116432897 zpos=2.303373861414658 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6067738661490414 ypos=2.192389194256829 zpos=1.866156669057471 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.083158071400519 ypos=1.7048250177421165 zpos=1.487968840480474 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6769358916533372 ypos=1.749938534882916 zpos=2.2151327769847966 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0097075815296273 ypos=2.402927930912443 zpos=1.5617394003044602 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4984786808280406 ypos=2.118270750282995 zpos=1.7535640217379675 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8783419775347943 ypos=1.9575184947054525 zpos=1.8300127632697427 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3207412115556822 ypos=2.1542746460169977 zpos=2.439245082849833 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.85343842596512 ypos=2.1788887074277548 zpos=2.195874415378216 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.456460095334337 ypos=2.1165240163412435 zpos=1.5714491846107463 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.992100968054512 ypos=1.8761292028017151 zpos=1.439878705133661 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4383215534958014 ypos=1.4672363133391797 zpos=1.9810751908083721 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8594495291183524 ypos=2.4151163358881806 zpos=2.052436917084886 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.450930522300941 ypos=2.3452425437905453 zpos=1.6779050422686355 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.067245113327741 ypos=2.1022416716762087 zpos=2.1929954560233416 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.865168728442368 ypos=1.6069223324129402 zpos=1.835335072183216 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.906149370430683 ypos=1.5064149194826102 zpos=1.7021535873127518 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9495702690076606 ypos=2.055601652000586 zpos=1.8089504634084852 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.704504199676881 ypos=1.7695881798742739 zpos=1.6747684890881616 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.668726790083547 ypos=1.5037246672987958 zpos=1.9702530519066361 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.74429250088694 ypos=2.002379460761161 zpos=2.2230222340035666 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.438387213412512 ypos=2.24929063173187 zpos=1.8239965417370427 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.945286757934672 ypos=1.797589576560148 zpos=1.7586526509792855 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.121166006779364 ypos=2.062321583757711 zpos=1.406244223948315 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.655991898774088 ypos=2.203075718885966 zpos=1.813155030175248 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3288688265484736 ypos=2.0846195924512894 zpos=2.554751368900243 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.055524177361906 ypos=1.8839758317766757 zpos=1.6774576840856064 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6830975422054824 ypos=1.7287680167980752 zpos=2.130705206079431 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.554695953250502 ypos=2.070808090020666 zpos=2.3554049068806915 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.719268997781064 ypos=1.9949426073170582 zpos=1.9565858961270748 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4245345848583737 ypos=2.0667616314986095 zpos=1.5687469969413659 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.997025725163094 ypos=2.081768529063928 zpos=2.1290568009988053 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.184145592234329 ypos=2.1448586744502283 zpos=2.226869069449034 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8694100375508556 ypos=1.8157315434193666 zpos=1.805744503579737 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0655852218667956 ypos=2.1856531294694173 zpos=2.152567410180823 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.374694492073793 ypos=2.360751764656036 zpos=2.4109739266945462 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2364774789743387 ypos=2.264816107644088 zpos=2.434630971535128 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8399155831200207 ypos=1.8057721719256445 zpos=2.1941703238648693 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.364405658756169 ypos=1.5627509105845254 zpos=2.2295160412485866 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.707836061705091 ypos=2.4480379440072464 zpos=1.706432101636037 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0336286295068327 ypos=2.034378882573436 zpos=1.7600183418225992 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.245886079649499 ypos=1.8801781177402594 zpos=2.065056436547157 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5910497703212387 ypos=1.8626893662904822 zpos=2.0385957886761292 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.813950440187617 ypos=2.0439607207007087 zpos=2.4440742073188946 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.231324393313657 ypos=1.8620360398616627 zpos=1.6307134447445966 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2545475044440924 ypos=1.4807642411060122 zpos=2.18672245716019 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2752346456332395 ypos=2.2017102613008688 zpos=1.4744746983832666 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7655160428131933 ypos=1.8369607865115356 zpos=2.4102591487010034 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.589769770252501 ypos=1.5032267654950093 zpos=2.1257730234289864 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7290780738702742 ypos=2.160981117869355 zpos=1.7034237271747792 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.990812908018505 ypos=1.9763502994652882 zpos=2.1880239957565863 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5626597458785114 ypos=1.8802448079832952 zpos=1.6053085128628135 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6855940833973437 ypos=2.120047351829517 zpos=1.671979675101135 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.429459342473609 ypos=1.7103715367678327 zpos=2.05398149565776 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3905062060583724 ypos=1.6708966028904353 zpos=2.3596940992138324 CorrectPosition=true CorrectZ=true CorrectRoll=true heightProfile=danio widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 8.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0.1 -tend 100.0 "
OPTIONS+=" -CFL 0.4 -nu ${NU}"
OPTIONS+=" -levelMax 7 -levelStart 4 -Rtol 5.0 -Ctol 0.1"
OPTIONS+=" -bMeanConstraint 2 "
OPTIONS+=" -poissonSolver ${PSOLVER}"
