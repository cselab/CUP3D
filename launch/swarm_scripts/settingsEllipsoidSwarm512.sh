#!/bin/bash
NNODE=256
BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}
NU=${NU:-0.00004}
BC=${BC:-freespace}


FACTORY=
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8079147869767676 ypos=2.428248697006022 zpos=1.8825573529659696 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.184358180722814 ypos=2.16614301661999 zpos=2.354049812634989 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.59688787289342 ypos=2.285189348474806 zpos=1.8620724715196435 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.842012127476236 ypos=1.7179262148913235 zpos=2.1472479053339457 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.903404748073509 ypos=1.7175417289643082 zpos=2.3598712774142627 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.052821627079624 ypos=2.0577361529041758 zpos=2.413036731746632 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.059520283124013 ypos=1.998564915242297 zpos=2.3435957957675924 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.274928731247296 ypos=2.0166499242526856 zpos=1.9469921790727271 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.349341663742045 ypos=1.979966918936756 zpos=2.059097047765718 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.913363701098108 ypos=2.296170568922477 zpos=2.2576980358978718 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.803133419615935 ypos=2.1173668905309464 zpos=2.4138045345736026 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.911891449698668 ypos=1.8726538335621146 zpos=1.6888484601715927 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.019626774721702 ypos=2.379910276052397 zpos=2.0182497545186324 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8223724491261435 ypos=1.993358531285214 zpos=1.7388740491914365 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.811350809550636 ypos=2.0259059779034856 zpos=2.07505530690124 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.657467307459874 ypos=2.0936992333209665 zpos=2.4669336608138672 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.050540064221057 ypos=1.7841777122760607 zpos=1.9103620954890042 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.468159119595299 ypos=2.2573461655802474 zpos=2.061027109348931 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.372343594905772 ypos=2.048094061175167 zpos=1.6638115106836646 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.10993316691059 ypos=2.032235186308089 zpos=2.230328826234416 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.749600976284597 ypos=1.7438365705448868 zpos=2.3143944256522326 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.415456030201676 ypos=1.990377543351621 zpos=2.2970589647125457 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9215621503309954 ypos=2.078475945254817 zpos=2.065772133048021 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.6865195326244464 ypos=1.899186967882114 zpos=1.9145584893536185 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.403834549164142 ypos=1.8740283359637795 zpos=2.0762933206690475 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4263667607049686 ypos=2.1435358251022363 zpos=1.5919891099274401 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9549888058826683 ypos=2.333553500221468 zpos=2.237235984925145 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.521357649864497 ypos=2.141548950953678 zpos=2.1110816049837156 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.758423236956139 ypos=1.9628255450719994 zpos=1.6045689074645257 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.600275274415542 ypos=1.5832037544411937 zpos=2.0660326230612567 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.573618206586102 ypos=2.2564155122677403 zpos=2.1895702215687067 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.16023611049437 ypos=1.6584009410426785 zpos=1.922486400614563 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8228055360295317 ypos=1.6507543334940387 zpos=2.060957526852899 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.491137604915578 ypos=1.8910236171627393 zpos=1.9813336530861265 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.127630516103028 ypos=1.917222942124279 zpos=2.297755348080485 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.298833595318516 ypos=2.0758576915264313 zpos=2.083100729307597 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2002697921119654 ypos=1.900074938787899 zpos=2.208164424465154 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9920249313881526 ypos=1.854135907993944 zpos=2.231295581931857 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2182247040342173 ypos=2.2272311134302587 zpos=2.032677965751056 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.14271638980883 ypos=2.002899838080577 zpos=1.9200015163168882 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.919809273565975 ypos=1.8461891914921955 zpos=1.6309170453452122 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3642174190093304 ypos=1.6926645232311284 zpos=2.0288835940382324 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.945645918902422 ypos=2.1427524637570947 zpos=2.0512232632093745 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.495662347790434 ypos=2.184289723508462 zpos=1.6142040245806217 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.195761560626964 ypos=2.023562696019655 zpos=1.835175098006595 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.767709868739606 ypos=2.337991948561369 zpos=1.848476931179033 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.041641101846054 ypos=2.339667079705495 zpos=1.7204795613182196 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.722986169146105 ypos=2.308629641060814 zpos=1.9140509179549139 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1692240692511198 ypos=1.9954761431595354 zpos=1.6841044150748714 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.057704496419532 ypos=1.8490406364603014 zpos=2.1290597524194745 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.537482616722383 ypos=2.284553925172821 zpos=1.6724587136445586 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.550801989722714 ypos=1.8320647010403348 zpos=1.6966877021976916 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.1904734657734215 ypos=1.9609401420792663 zpos=1.9365975330926362 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.565421940264866 ypos=2.429464129171681 zpos=1.866094818821195 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.223382628194689 ypos=2.2181222750490495 zpos=1.857448605806057 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.090335328581515 ypos=2.093029202297687 zpos=1.7649103274804392 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.103342075857179 ypos=1.748448063517058 zpos=2.0042648321993846 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5432831475680233 ypos=1.6430084978448414 zpos=2.016793966305439 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.185611408457489 ypos=1.7694733581655655 zpos=2.2465947113212184 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.127522206421871 ypos=1.7989218074939062 zpos=2.097718612742245 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.7984590344425495 ypos=2.3054646705693567 zpos=1.6846448911072611 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7559649168141447 ypos=2.137715625920036 zpos=2.28700955503016 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.988368367620419 ypos=1.6232171850258674 zpos=1.9160644561986757 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.971002765269604 ypos=1.728246375533815 zpos=2.2060260592625323 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.156870197707276 ypos=1.8978234634581623 zpos=1.5652239028790964 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.782957266457263 ypos=1.6149728835857773 zpos=1.8788644918436301 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6600266781997908 ypos=2.0291677277079225 zpos=2.1123846012220953 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3596904708376636 ypos=2.1567646869904915 zpos=2.2460388448012747 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.014910047410173 ypos=1.7149931885723346 zpos=1.8518067018057662 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.290353817614597 ypos=2.0757865441955965 zpos=1.6895024296192669 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.269634502082997 ypos=1.6537235230077816 zpos=1.9150331129484293 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.8078173345335316 ypos=1.9823129890881797 zpos=1.800816109110322 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7302176122512622 ypos=1.9814244062201098 zpos=2.346474889761163 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.198974145010384 ypos=2.0472332384811875 zpos=2.1692605810089716 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.340936183514966 ypos=2.176819052053817 zpos=2.0426763551736653 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6148085824092058 ypos=2.166992038563033 zpos=2.3410703272412596 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4022014975567356 ypos=2.154358594138242 zpos=1.6899621473975972 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.165721095798961 ypos=1.8392023966311992 zpos=1.7566943318994355 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.387570236323515 ypos=2.221676653856228 zpos=1.6749991933913206 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.2740735852638405 ypos=2.083120484116949 zpos=2.0477958288821037 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3154765322901443 ypos=2.099623825397835 zpos=2.0261189862464977 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.308109919815035 ypos=2.1095684280926332 zpos=1.9736196583061114 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.77374079362241 ypos=1.828203080615267 zpos=2.410983823690278 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.565208146823208 ypos=1.9248236679071997 zpos=1.9989815992757922 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.577764840599423 ypos=1.7722462511629145 zpos=2.059669414599621 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.62116454091616 ypos=1.9498884335348783 zpos=2.34511042405938 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5570921358881145 ypos=2.365210710165168 zpos=1.8672887218886598 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.703496761830449 ypos=1.9779219669308907 zpos=1.973141038385451 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.5070105483626435 ypos=2.019556465305541 zpos=2.022240636603784 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.627648102800681 ypos=2.0857560442396905 zpos=2.25182653614543 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3052221047305514 ypos=2.2150446325425763 zpos=1.9116723818862302 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8954114875564154 ypos=2.2653227936878824 zpos=2.230195804753024 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.174395184963317 ypos=1.9712224775594314 zpos=2.390912852151006 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2025194277359437 ypos=1.9624148578584828 zpos=2.0932807370152293 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2663208537654724 ypos=1.599062614229148 zpos=1.9878181817512914 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.502375308098268 ypos=1.8335042145857434 zpos=2.145291593272142 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5755521890974493 ypos=2.0286677632130563 zpos=1.5803537326297759 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.079597185334881 ypos=1.895534853081218 zpos=1.9822091732207856 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.8695157539380687 ypos=2.232991707928447 zpos=1.889642438144818 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.518737144730129 ypos=1.9361250370355467 zpos=2.1020953813428775 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.106114234566276 ypos=2.202504757849789 zpos=1.9497240423104036 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.276341432365332 ypos=2.321460242282993 zpos=1.9440346452522277 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.753350996858531 ypos=1.7057040334302198 zpos=2.2124145360693386 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.191181922722685 ypos=2.013098920951855 zpos=2.212665092462594 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.152982483288519 ypos=1.9938218101649514 zpos=2.2553591667420636 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8365469901569003 ypos=1.7886159351808206 zpos=1.6404554702859708 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.244276309774592 ypos=1.6407563461795112 zpos=2.063567993694416 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.724108236855024 ypos=2.033674460555426 zpos=2.177033251039868 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.682069913473789 ypos=1.6855904663370156 zpos=2.3090304278582012 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.105789084350042 ypos=2.2559832081279176 zpos=1.7930480713402352 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.3998540389998135 ypos=1.7362646536926825 zpos=2.3216347795505206 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.409606568093418 ypos=1.6574336968960561 zpos=2.1573175868088907 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.6480602887667373 ypos=2.0177557799996655 zpos=2.0314003841275845 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.111359553158718 ypos=1.567716532953709 zpos=2.0840602528950773 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.374983553099897 ypos=1.5528954141364397 zpos=2.1382586860810884 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.544597616244379 ypos=1.9129035549771958 zpos=1.6804935105244836 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.872992067097497 ypos=2.219344352075884 zpos=2.1495194624432177 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.402948949221425 ypos=2.121956127729501 zpos=1.8497777064155825 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.2614245048303285 ypos=2.144213586919012 zpos=1.7384325927636017 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6917028919902815 ypos=2.046351149186428 zpos=2.1902152573413995 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.882907261204931 ypos=1.9137684117866365 zpos=1.992845740068377 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7224748704027366 ypos=1.9647178278743445 zpos=2.306918329685248 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.266052059211592 ypos=2.0827500936015504 zpos=1.907619534494065 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.0853659672411995 ypos=1.6285604509878036 zpos=1.8206704445489676 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.604922476607994 ypos=2.12577551780146 zpos=1.669485355694509 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.898891060531133 ypos=2.1966826933511165 zpos=1.850088166563831 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.70693244067608 ypos=1.778030548084158 zpos=2.1052021235630995 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.121346283856853 ypos=2.133552799384634 zpos=2.465532261819641 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.440961665677719 ypos=2.0044690492731365 zpos=2.2073517501371342 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.045689468995085 ypos=1.9322529077661235 zpos=2.132026806031629 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.664252440928976 ypos=1.5735087061096238 zpos=2.171582721393073 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1058391285382814 ypos=2.2736042698136676 zpos=2.1478697425790862 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.249576643380834 ypos=1.7829684821096934 zpos=1.7033549462409228 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8489809334512675 ypos=1.8092554892073798 zpos=2.42237654161665 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.370123697592832 ypos=1.591397611332539 zpos=1.9361768399908286 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.282823119689238 ypos=1.6229942552726244 zpos=1.7730005206630979 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.416727277164141 ypos=2.1228968077343153 zpos=2.1394982968563623 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.682624687681023 ypos=2.0863381130628875 zpos=1.7258988128440427 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.102622815858624 ypos=1.6784801051680287 zpos=1.8441021668557918 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.090474920853054 ypos=2.2791835032798944 zpos=2.376549019347066 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8723127384925737 ypos=1.865409367921688 zpos=1.8287369341536868 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.437847729307342 ypos=2.2561626145474887 zpos=1.7434490684431094 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.371568947659955 ypos=1.9830838436525755 zpos=1.9298148947924074 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5010540793641023 ypos=2.0192467902858326 zpos=2.0813372749636323 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9292598497108115 ypos=1.9247768208924299 zpos=1.960048978441515 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.617113122962163 ypos=1.8839082993132321 zpos=1.7703487625163754 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.097254404203105 ypos=1.8906900548385124 zpos=1.8723130990691046 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.858410194152713 ypos=1.9689800219353695 zpos=1.9420700128478283 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.233102814291892 ypos=2.224901791686111 zpos=2.3417719716732313 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.646184849841421 ypos=2.0730627579128176 zpos=1.8155060303436392 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2792951723762482 ypos=2.3038279364971532 zpos=2.0810558002443447 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.037400422338206 ypos=2.1504840556467353 zpos=2.0135961226405956 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.148040981112208 ypos=1.740401911679812 zpos=2.0884093078003185 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.673843839704001 ypos=1.6325067547975443 zpos=2.049263082165515 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7304299744305554 ypos=2.143864599315122 zpos=2.078830815092053 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.703211087319836 ypos=2.2698379609296464 zpos=2.2846501493478715 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.92890803902779 ypos=2.166926134652347 zpos=1.5804959941676586 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0281118486066894 ypos=2.167830811768882 zpos=1.7067510687408658 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.365605832530008 ypos=2.052727333355389 zpos=1.9655875789524555 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.907077136874524 ypos=2.316046027043784 zpos=2.3590489449955383 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.199952611695709 ypos=2.2448520527115905 zpos=2.0246918805153067 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.855714825957479 ypos=2.072798459062482 zpos=1.7218974661634014 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.761448218052623 ypos=1.8275906415055676 zpos=2.2545729017245613 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8050490044903476 ypos=2.012864587732484 zpos=1.5169663698475166 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0009562672170937 ypos=1.781367707056862 zpos=1.7825648857946506 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.0917295418309525 ypos=2.218876604137272 zpos=2.249973619028204 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.832280534459388 ypos=2.3834946326422015 zpos=2.045592115585231 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.243809646311519 ypos=1.777445098235471 zpos=2.2918056207078727 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.057638757622132 ypos=1.8289326430108757 zpos=1.9224348752418576 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.887500545262121 ypos=2.033597020752486 zpos=1.6846067811139216 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.44588152584305 ypos=1.7141118393427375 zpos=1.925084563187215 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.128692663975647 ypos=2.0031558289202236 zpos=1.9016671150201745 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.803415326333346 ypos=2.1922226158119558 zpos=2.2495948165333086 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.06822712795803 ypos=1.651282694538173 zpos=2.137738793557452 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.951894380957896 ypos=2.3022577680999814 zpos=1.6204269295232476 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.658519476013837 ypos=2.0194664906523183 zpos=1.8915996666943726 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.506084745533578 ypos=2.2383076539502924 zpos=2.1055012213105058 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.50658336826157 ypos=2.110430877518574 zpos=1.9980870168763776 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.237191294561343 ypos=1.6683434031128455 zpos=2.168289009565067 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3355728253295513 ypos=1.7067559299540904 zpos=1.773636618291747 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9360105247884616 ypos=1.8260132781660006 zpos=1.6962058066691066 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.651600181426676 ypos=1.6476145775573419 zpos=2.186115023573062 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.158811085584262 ypos=1.7491648344738522 zpos=2.0260264962356 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.255676851389218 ypos=1.5069646418117466 zpos=1.9759602725688121 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.621663346343116 ypos=2.4009526232608946 zpos=2.1900131898128765 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.513972934934065 ypos=1.9195829986043464 zpos=1.8729722979582386 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.143222344034789 ypos=1.8691402078389219 zpos=2.2035193258344368 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.460065569095672 ypos=1.6978325404956633 zpos=2.0439863105874143 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.628778234666262 ypos=2.212189560343795 zpos=1.643920087578263 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.921613422704164 ypos=1.9204532030048636 zpos=1.9908654679762723 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7842833736029298 ypos=1.8515530548480825 zpos=1.571998610340798 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.634821207223654 ypos=2.2314043478172607 zpos=1.784217807853387 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.0688895886411185 ypos=1.901621740031666 zpos=1.9172870539653801 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.456349757939348 ypos=2.349610705726241 zpos=2.3211337373443457 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.984112515280339 ypos=2.3499096658429552 zpos=2.1838415931655217 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.678451854165452 ypos=2.3080669769805273 zpos=2.0574276507782674 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.323863224876307 ypos=2.0281793220082687 zpos=1.999523303242728 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.9838573712289 ypos=2.3325123238996857 zpos=1.940236060154651 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2722985185551 ypos=2.3358152489958077 zpos=2.024367925557096 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.188052460376784 ypos=1.7262149784091623 zpos=2.1889071745526905 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8729330104342523 ypos=2.018268712167271 zpos=2.1474422337347074 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.921076970328804 ypos=1.844921550526153 zpos=2.3454013944285212 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.581034197305937 ypos=2.284830095051859 zpos=1.6825868251215734 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.484265317516221 ypos=2.178354242594007 zpos=1.9090623562659073 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6833296670382176 ypos=1.6062570504925817 zpos=2.2530193275594446 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.213913999409867 ypos=1.9143016932643928 zpos=1.965346706554471 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.631998686293709 ypos=2.3479648294986104 zpos=1.7475968528868895 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.807017373578872 ypos=1.810530754327861 zpos=1.9767994858474052 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4199051986551763 ypos=1.7501700103483213 zpos=1.7674996362711297 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.7754844303644304 ypos=1.8173146159016986 zpos=1.8610827886532046 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.063108899639007 ypos=1.9750450010988443 zpos=2.0886238816247795 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.188971182097101 ypos=1.952742046673319 zpos=1.6886374134559017 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.442642238322839 ypos=1.6875420182588703 zpos=1.825273962383358 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.1839927635003775 ypos=2.418013470669761 zpos=2.187854086783929 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.1693007824180714 ypos=1.7258443128982337 zpos=1.8977216134699557 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.105055626307838 ypos=2.191149007990258 zpos=2.258922739676788 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.051713017326805 ypos=1.5479313940994777 zpos=1.8365243864643848 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9544179403279838 ypos=1.9023552061044575 zpos=2.466790486067235 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5832851861081942 ypos=1.730538162942937 zpos=1.7267851096814428 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.234142471796526 ypos=1.9750964339276174 zpos=2.240081788508115 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.0052040304011545 ypos=1.944151729106908 zpos=1.8440568547051932 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.483121238998948 ypos=2.1672394915789273 zpos=2.0677670142810474 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9065185468936576 ypos=2.365748441251794 zpos=2.188781202218809 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.4428552730237385 ypos=2.0572213414967626 zpos=1.8742097378676412 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.127109249057125 ypos=1.7465413949690607 zpos=2.10104063470567 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.472182232398391 ypos=2.185292211128731 zpos=1.64089184757148 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.319012205069447 ypos=1.7316529074930072 zpos=1.9448903290077215 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4142910111641496 ypos=2.2065695431937056 zpos=1.7359151430590067 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7777074068286782 ypos=1.95804039083966 zpos=2.1308514072029214 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.072249030302925 ypos=1.7432932063471478 zpos=1.7609831489762398 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.086611524537103 ypos=1.9010057846455286 zpos=1.7923038308331884 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.592623379292726 ypos=1.9657895910157117 zpos=1.809248412740943 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.693644838771433 ypos=2.206885172146239 zpos=2.413947060745335 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.743794192128163 ypos=2.35043881032061 zpos=2.253574511686992 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.075121923703805 ypos=2.0733980280661566 zpos=2.3743293921584234 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6057059599574384 ypos=1.7320761002130536 zpos=1.6437642478077057 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.961788750610329 ypos=2.0359756752452633 zpos=1.764791336990276 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3948044956624983 ypos=1.6614752097138843 zpos=1.7341608826383614 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.170270787626045 ypos=1.9970583372395896 zpos=2.0010354235058636 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.3093931301633175 ypos=1.7920057662754312 zpos=1.9031282497630597 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3118175033354924 ypos=1.8607643014766828 zpos=2.3037328152368226 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.291256357309889 ypos=2.213776475147093 zpos=1.8439519338464216 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.066902452493886 ypos=1.7570480689466903 zpos=1.9329699012732682 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.798690557580705 ypos=2.138387962144691 zpos=1.8192195602528372 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.429469732745844 ypos=2.109639504768614 zpos=2.2120785476184532 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.067405830306876 ypos=1.6801314052069043 zpos=2.0104037604114877 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.410177318412309 ypos=1.707481623416196 zpos=2.116078762715288 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.969104696403375 ypos=2.4118496231426647 zpos=1.812322822406602 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3180168338114253 ypos=2.0722716030700843 zpos=1.9905885385736428 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6598680257688243 ypos=1.8720381715025032 zpos=2.358502811317285 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.585440670212036 ypos=1.9365874692478378 zpos=1.7939846697060045 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.957830404040312 ypos=2.163105232091584 zpos=1.7932811856335407 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.308606658517798 ypos=1.8874565162930825 zpos=1.8653901196179583 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.7717655379127555 ypos=1.7030035487252861 zpos=1.678582229588903 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5845589243604428 ypos=2.0843772374925416 zpos=1.6460746373835407 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5964895950439764 ypos=2.1753293563312774 zpos=2.1832861016706526 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5160252198337383 ypos=2.364413479783954 zpos=2.174241556337682 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.010093812460374 ypos=1.6573646675213283 zpos=2.009151556238767 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.34829913354955 ypos=2.1609866910651667 zpos=1.899983800688571 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.8647679373635935 ypos=1.9485042632896554 zpos=1.8826885855712685 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.078167055026483 ypos=2.000242624918598 zpos=2.1342629411862384 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.764469884022493 ypos=2.3851725912233075 zpos=1.971194022930737 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.1971525521818664 ypos=2.1074359216806173 zpos=1.7852245930061552 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8316689178358 ypos=2.132881959376942 zpos=2.419439972638721 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.278399366194198 ypos=2.057143555770488 zpos=1.8746314757304305 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.1696160084195055 ypos=1.8547364361693999 zpos=2.0050535837365304 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.571229504716797 ypos=1.9408696810960002 zpos=1.6226019781334722 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.918785461568834 ypos=2.0661932088176136 zpos=1.7382885442215377 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.166489055900776 ypos=1.879129966515898 zpos=2.457661733773544 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.923935668359014 ypos=2.256338432233455 zpos=2.202082066571826 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.872487820248042 ypos=2.2243233615528046 zpos=1.8386760736179113 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.150694525186368 ypos=2.171046479678502 zpos=2.021869645252172 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.09812283524116 ypos=1.6265286873262386 zpos=2.2888405635423292 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9014622548176017 ypos=1.8265933297313575 zpos=2.1789024824748817 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9424513083153916 ypos=1.7502224390503283 zpos=2.1596583830459335 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.8924488671480155 ypos=1.6988749445493436 zpos=1.7182502728032125 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.6921185429145735 ypos=2.137711191181036 zpos=2.349216672236426 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8075667402559357 ypos=1.8083294589961363 zpos=1.791864519444346 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7935613927755094 ypos=2.253790855454257 zpos=1.9678207628632647 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.7872583763053083 ypos=2.1149335117687054 zpos=1.9822178540718147 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.287865989467146 ypos=2.472806343170607 zpos=2.1145217916787264 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.973612221303782 ypos=1.780864105695732 zpos=2.007021861693806 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.399890425486514 ypos=1.8767733857443003 zpos=1.8228740216465273 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.014610303842755 ypos=2.377760110033627 zpos=1.8842850415423082 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.4718916925150785 ypos=2.095120182702228 zpos=1.658903784412473 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.499315614757714 ypos=2.4090161903181637 zpos=1.9867259135880906 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.890585565930525 ypos=1.9840032304685926 zpos=2.319483827230744 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.93960672024851 ypos=1.7605871665926136 zpos=2.2901331011725863 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0738978036029976 ypos=2.2112386749015935 zpos=2.1345156180678937 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.4379810821665755 ypos=1.553478392518841 zpos=1.9244065280204516 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.34509876923662 ypos=2.252251115319709 zpos=1.8452555363940126 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6325757464350725 ypos=1.882716021812483 zpos=2.2834520756323795 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.941185555152003 ypos=2.0791819096813526 zpos=2.120831903838515 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.560318453075789 ypos=2.416191520726226 zpos=1.9507669781910104 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.775608180912262 ypos=2.2513230443873655 zpos=2.0740069797012777 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5306876712796753 ypos=2.015359722323327 zpos=2.237245493959878 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7040471787946 ypos=2.3964332580377716 zpos=1.8230821450897863 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.328245489315053 ypos=2.2817982275057362 zpos=2.1419292267740158 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.710321761491394 ypos=2.091927838330933 zpos=2.3631973494476712 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.667995311548616 ypos=2.1531021854454253 zpos=2.219562423892178 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.066445281181103 ypos=1.7005160278682725 zpos=1.64035079441582 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.235056287759935 ypos=1.8634995692075134 zpos=2.36278424290411 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.407484210030664 ypos=1.7639234290108292 zpos=2.089807728724348 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0449198673974784 ypos=1.8951685042487836 zpos=2.0363725976533207 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.21057772874801 ypos=1.8612047158720997 zpos=2.049826395239831 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.354362582010239 ypos=2.437551512843211 zpos=1.8381116979548513 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.304854892625314 ypos=1.974398086732561 zpos=1.8446844444796888 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.575522090933245 ypos=2.1096432878263296 zpos=2.2068497944198326 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.089109541294202 ypos=2.103290750560899 zpos=1.7384379977885942 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.131065242596902 ypos=1.6932266931896915 zpos=1.988948046199326 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.151180939223433 ypos=2.072759722217613 zpos=1.9489393325721556 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7420879867278054 ypos=1.5500214341612415 zpos=2.0265624017572463 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6065438463020443 ypos=1.8180378087779496 zpos=2.0882848136278582 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.054647737458673 ypos=1.9867025715464748 zpos=2.1636052823306535 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.482397857075971 ypos=1.968102050542925 zpos=2.0558164721709407 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.618842603367814 ypos=1.7973685909912909 zpos=2.277715407633062 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.120611612116852 ypos=2.097009753047722 zpos=2.0859589986020826 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.748619612746271 ypos=2.0248906678217495 zpos=2.219235920752679 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.963534402306426 ypos=1.8529872179430182 zpos=2.3696328872325805 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.7226361181510783 ypos=2.0189624394615286 zpos=2.171193324396583 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4787257886796543 ypos=1.8335852444037015 zpos=2.3734174646444006 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.815212469606152 ypos=2.2723531660823544 zpos=2.320972089333003 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.708245378002439 ypos=2.183192757760626 zpos=1.958381407539779 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.343829051172421 ypos=1.7573651593322035 zpos=2.181548610562234 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2490333816749644 ypos=1.8450320297192655 zpos=1.7147414357510775 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.484365066185977 ypos=1.6244832285526818 zpos=2.22303211592858 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8941301098199785 ypos=1.9664385215806837 zpos=1.706058602886058 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.865430596906009 ypos=2.15456674873299 zpos=2.141427593722304 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.796828786074208 ypos=2.199807754710874 zpos=1.6301921558030377 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.853754784074624 ypos=2.0392326432445036 zpos=1.554711530483071 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.288492706339182 ypos=2.1591246142782294 zpos=2.1816381438590122 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2922066206565073 ypos=1.6522966935998367 zpos=1.9696777697598875 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6446191620850477 ypos=1.9119780756077733 zpos=1.5900393385818572 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4529787221608803 ypos=1.900970370641656 zpos=2.312198826644575 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.133281646584865 ypos=2.3903992892624344 zpos=1.836442183075253 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.63484152416954 ypos=2.4542381256049177 zpos=1.9540450177314135 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.300180583320981 ypos=2.202791453843588 zpos=2.300393368104041 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.2141526493596135 ypos=1.7703614760770185 zpos=2.0574917462254922 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.936970851610363 ypos=2.052142243053458 zpos=1.9212465404696406 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.812056891486599 ypos=2.425162014669382 zpos=2.22563729564167 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9983845212786266 ypos=1.9522668644090513 zpos=1.8896844954430765 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.807012683807615 ypos=2.089668357667464 zpos=1.981528227444622 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.151763461150486 ypos=2.2799394133591244 zpos=2.2247292558230267 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.517111453113103 ypos=2.0503870142523914 zpos=1.9479211050341214 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.81637375175831 ypos=2.0322097178602383 zpos=2.318339931353334 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.600194680678254 ypos=2.1892349933146145 zpos=1.9913994028267554 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.127497337620703 ypos=2.359915517179978 zpos=2.2151800494307285 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9215562323638684 ypos=2.200209569502138 zpos=1.7938686956090042 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.0448986015983275 ypos=1.5761941517842677 zpos=2.243482750769074 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6446520050280866 ypos=1.9908044243859966 zpos=1.5707671670995402 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.599291795611827 ypos=2.0350413910164447 zpos=1.7436869428164559 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.632100179963594 ypos=1.9258042223662255 zpos=2.0380521788221957 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.103826502004343 ypos=1.6085509903959363 zpos=2.1516381485808025 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.044786092189666 ypos=1.7926428164802863 zpos=2.141132210648026 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.529631974193037 ypos=2.3854200587100256 zpos=1.7461762479826222 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3501285011852606 ypos=1.8465165426197518 zpos=1.9162042114563675 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.362721741526364 ypos=2.0070062091177725 zpos=2.405291527238342 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.548101207209289 ypos=2.3077308172225015 zpos=2.034850681238871 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.567038191790325 ypos=1.62629608409275 zpos=2.2855854607697768 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.367905554397409 ypos=2.0048707255058433 zpos=2.3716661654586333 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8540521200223203 ypos=1.693046529788251 zpos=2.241515065802544 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.858198023955059 ypos=1.7705887567732206 zpos=1.8724568363328788 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5331707795068206 ypos=1.7638510282447435 zpos=1.8480929374996846 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.095325571018664 ypos=1.7779535983043653 zpos=1.7150553175074 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.0885501817922005 ypos=2.41997606685519 zpos=2.1338433318067525 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.400907032638651 ypos=2.010666689061971 zpos=1.7632712842298999 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6673805689125727 ypos=1.6399680092306723 zpos=1.7965386954945928 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.38202014595454 ypos=2.020124924247803 zpos=1.7739743917377233 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.7359443344366 ypos=1.9052698419072387 zpos=2.357744821651849 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8436669253886246 ypos=2.3909609222805988 zpos=1.703711233601687 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.069093696165941 ypos=1.8053880062594145 zpos=2.2453736053633513 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.432719669718462 ypos=2.318332227342496 zpos=2.2440365532455284 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9488344510098026 ypos=2.4021563711723735 zpos=2.1063223662757697 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.256110982460637 ypos=1.8881053262181824 zpos=1.6559250521832842 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.263199801425705 ypos=1.8604620563618535 zpos=1.852679054559836 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.622410484554499 ypos=1.7881339311672024 zpos=1.7448226891709266 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.449041726375554 ypos=1.6853371783858317 zpos=1.9995369529289888 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.047102310147016 ypos=2.283083821714642 zpos=1.9781847645060604 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4806912021155783 ypos=2.4059347739786445 zpos=2.1604579027862014 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6707091054160257 ypos=1.5920039865501312 zpos=1.994995915702949 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.05018187182393 ypos=2.4477679687958607 zpos=2.0676794630198443 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.938262653083992 ypos=2.116296607871137 zpos=2.24443088791484 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.671025017052637 ypos=1.8924944003766433 zpos=2.0665193752567594 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.298831026513219 ypos=2.4689286856299884 zpos=2.0273472277971356 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.4453840907716735 ypos=2.3341279994895183 zpos=1.6547953605183114 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3753046794213275 ypos=2.3186517130070365 zpos=2.154824093816325 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.610144635495939 ypos=2.02028588520478 zpos=1.9825776891860214 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.156011729106034 ypos=1.5252811774277375 zpos=2.016689599530662 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.813231240625563 ypos=2.0509211427437224 zpos=1.7556226626724638 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.391646715327312 ypos=1.7228468061590259 zpos=1.9267471347605287 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.050612674670957 ypos=2.094107136825689 zpos=1.8831267579422335 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.404898224087502 ypos=2.2111537882379726 zpos=2.1737168360651506 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.802863826515042 ypos=2.04207582583334 zpos=2.4095972376418886 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5712271056422646 ypos=2.11796786902836 zpos=1.894471170746664 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4782509824441856 ypos=2.39737458738762 zpos=2.0652243464839617 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.915611211133532 ypos=2.180451763965199 zpos=2.314841546898027 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.997551532027662 ypos=1.9690949203762569 zpos=1.638555377540116 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.249016401300312 ypos=2.388874015805935 zpos=2.093889029362683 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5235148636127778 ypos=2.0788574267199444 zpos=1.92815099609989 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.449305834419103 ypos=2.0769245388462223 zpos=1.8054032434933946 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.210774250361027 ypos=2.0446565804730725 zpos=1.8885448228670287 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.932697400014808 ypos=2.1645484034779505 zpos=1.684266831930766 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.094893845600577 ypos=2.0410068305106353 zpos=2.339885020013467 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1359705898125054 ypos=2.2408249230402344 zpos=2.074392166054812 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.868354822642015 ypos=1.819233376986871 zpos=1.8942770425571256 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.756119904132552 ypos=2.304388406910897 zpos=1.7446413516607082 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.62270778979846 ypos=2.0326082525650024 zpos=2.097199361305641 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.2390380523564 ypos=1.8844069804975008 zpos=1.9135194940982891 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.15413148689852 ypos=2.1714343430992913 zpos=1.960003432523069 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.2050791364040565 ypos=2.1232974847467747 zpos=2.0850607224661095 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.783025118224804 ypos=1.696083471007532 zpos=2.1520660079676337 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.420861323694002 ypos=1.969288016233871 zpos=1.5884973339881567 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.660362665550193 ypos=1.811890508003362 zpos=1.710009210065418 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.252540806315363 ypos=1.880348395627245 zpos=2.052899695805263 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.468639719833071 ypos=2.2143516093212208 zpos=2.3748304174253088 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7289429342670273 ypos=2.3908908201016708 zpos=2.005040656260524 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.575371642281873 ypos=2.0953769477767983 zpos=1.7630334120193287 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.893743128759353 ypos=1.683437928647264 zpos=1.8922102107697234 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.944514296831181 ypos=1.930320317495542 zpos=1.5313052759342105 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.393263878288356 ypos=1.593083424592996 zpos=2.2503355632903452 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0696238941711096 ypos=2.107793939563603 zpos=2.040762078392164 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.970643983419746 ypos=2.0704386198581792 zpos=1.7918581731447927 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9580192705548773 ypos=1.7854824740656077 zpos=2.199765810206114 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.8204198025176233 ypos=2.2590465907158728 zpos=2.0158952376958372 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.324751792909821 ypos=2.4015128848565936 zpos=1.8694040220188628 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.547960301809379 ypos=1.89699817208728 zpos=2.01032090809922 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.5252829054180435 ypos=2.27526777387086 zpos=1.9390130176811984 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.623378567305347 ypos=1.8160551728237855 zpos=1.577756064046985 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.817885740120024 ypos=1.8431130294887654 zpos=1.8624037099109583 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.9904278027391396 ypos=2.1579427098722785 zpos=1.8889261089659095 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.075591119994897 ypos=2.245042103822839 zpos=1.9565540902566467 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.992534984831783 ypos=1.9591914928959302 zpos=1.657125188445648 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.852093693769379 ypos=2.0466444802401114 zpos=2.0904356010581884 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.818966330554833 ypos=2.2440659722609064 zpos=1.8956748015806044 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.8213322187225285 ypos=2.1458667731659506 zpos=2.14164829656266 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.601866932532872 ypos=1.9247190864116897 zpos=1.7745963224980235 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.963994343433053 ypos=1.9889036994336964 zpos=2.2028142175247987 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.1676112581538485 ypos=1.781564161515893 zpos=1.8656680807010886 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.315971437413576 ypos=2.3444580224499014 zpos=1.7462141554335335 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5939317285295576 ypos=1.9610495050402152 zpos=1.7147221163998085 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.538848415668515 ypos=1.894730025038557 zpos=2.1789367685381977 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.143077106860899 ypos=2.0453896994529757 zpos=1.757379743769858 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2400225336070183 ypos=1.785758997947905 zpos=2.16654313793834 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.360152261215717 ypos=2.1287962089162775 zpos=1.9463535638451699 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.872063642387014 ypos=2.3089214672290783 zpos=2.029401316019215 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6675063496346705 ypos=1.9054235784275453 zpos=2.421490532515595 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.47879047814666 ypos=1.74808748395474 zpos=1.714480152669244 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9229485444125896 ypos=1.6844758121846966 zpos=2.0822351537970767 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.756532406206893 ypos=2.2695271830383597 zpos=1.6090511756203394 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.60976030631716 ypos=2.1511842989917955 zpos=2.4424682191140272 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.438656779309757 ypos=2.268838347147016 zpos=2.0064663502531133 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.969946263477596 ypos=1.7940565248306035 zpos=1.58517654736365 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.064687170663524 ypos=1.8504889415250967 zpos=1.7805184290202867 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4642230895096713 ypos=2.0103943654745104 zpos=1.8881776046304215 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.306842262437479 ypos=1.8486349527757189 zpos=1.9117742008841743 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.129740908575929 ypos=1.7059733048911925 zpos=1.9567725828922167 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6486463900553123 ypos=1.654553626299871 zpos=2.1789130120478335 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.7828089474069304 ypos=2.0701421515179055 zpos=2.017879658228792 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.028095571797951 ypos=2.2554614192765783 zpos=2.0650539148548863 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.811814242354607 ypos=1.9452008305726691 zpos=2.2751612100462206 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.382061777178275 ypos=2.379021149207178 zpos=2.138181692696223 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.130155112888911 ypos=1.9998370372872285 zpos=1.8199731515778113 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.153202861566664 ypos=1.6853224278882242 zpos=2.2131511200628196 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.145358988378645 ypos=1.947538709486455 zpos=2.4950189286277014 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.003861760636752 ypos=1.8967636496264249 zpos=1.8364916049004234 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.105100673246796 ypos=2.1767327204798512 zpos=2.2753587530141743 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9306828203997575 ypos=2.3544222822378105 zpos=2.307359648216111 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.7659669324591714 ypos=2.0806053091434933 zpos=1.8987658152763904 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7497468500939903 ypos=2.0070588680215895 zpos=2.0526704020662647 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.682120992060419 ypos=2.3446176943375816 zpos=1.9134229970616499 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.175809952027656 ypos=2.030826691880815 zpos=1.5480022810482943 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.65794908243758 ypos=2.0590700577285395 zpos=2.3119891037566496 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.395059071956029 ypos=2.0420825880514664 zpos=2.298017244093943 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6700169080846394 ypos=1.593828622968918 zpos=1.8357608793070803 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.874709308676443 ypos=2.109275204583076 zpos=1.6334702731367208 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.379663110461986 ypos=1.9786305705496054 zpos=2.4170401982653678 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.121519141314822 ypos=1.7261134492392718 zpos=2.40569292255546 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2436210174121043 ypos=2.157856083646852 zpos=1.6381007811271946 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8957743043140804 ypos=2.0227178023882706 zpos=1.8315915912383574 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9489200704393435 ypos=2.0478947796034794 zpos=2.010044410077505 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.010211381493621 ypos=1.7493127057255682 zpos=2.0287768331892346 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7444589709586222 ypos=2.0287839591107146 zpos=1.9135900318492869 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.985592948631207 ypos=2.2116437336733012 zpos=1.7976753378059458 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.113290073934148 ypos=1.9390444535837208 zpos=1.9855297109666166 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.592498056892302 ypos=1.741349988604352 zpos=1.858577068757762 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5437423389103415 ypos=1.6230307138521856 zpos=1.8894659814879802 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2556268035113387 ypos=1.8735769433525804 zpos=2.166303876716856 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.493454271761698 ypos=2.3871947731709264 zpos=2.11792317891826 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.505201800689818 ypos=2.2116539661105468 zpos=2.0715679788634582 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.329245514693419 ypos=1.8877535051754133 zpos=1.933334346716398 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.764130698783227 ypos=1.8075974292100445 zpos=1.6388428739800869 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.652036924851455 ypos=1.701082257756756 zpos=1.9008273004492164 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.374278910876386 ypos=1.854219145260695 zpos=2.1219018560063034 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.187900252882684 ypos=1.8338092890692517 zpos=1.8437970784487692 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.055040840257201 ypos=2.075239454252134 zpos=2.11109884537197 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.283018030722261 ypos=1.7643626614240207 zpos=1.9942279012680608 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9085227349197864 ypos=2.1252649452803305 zpos=1.5538894464792747 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.549128271878145 ypos=2.447946638021276 zpos=1.9938142511053287 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.2350474974644055 ypos=1.9273370919931474 zpos=1.6831045066682178 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0544918587226513 ypos=1.9238229892494738 zpos=2.1514148768374297 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.508707065180996 ypos=1.64736757546929 zpos=1.9262567451358077 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.321956208443871 ypos=2.0419062284816705 zpos=1.592793464870727 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7621108697660572 ypos=2.0947471362643832 zpos=1.829987816265009 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.063750420646761 ypos=2.1761079838418986 zpos=2.2202731286608484 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.959579248790489 ypos=1.936673342767929 zpos=2.0556650910848586 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.29406708931697 ypos=2.4174208064983187 zpos=1.9592830533299943 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1602688609850595 ypos=1.702372310094683 zpos=2.0702935704688663 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0607341265058152 ypos=2.068376068281556 zpos=2.193801562631977 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.009372763953987 ypos=1.5405457008431833 zpos=2.183507551365103 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4516933702656645 ypos=1.881747898499154 zpos=2.376279756460066 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2321681322530167 ypos=2.2264164106711264 zpos=1.845148291025294 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.314263867079226 ypos=2.248630628344882 zpos=2.3138804475317096 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 8.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.2 -tend 100.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.80 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 7 -levelStart 3 -Rtol 0.1 -Ctol 0.01"
OPTIONS+=" -Advection3rdOrder=true"
