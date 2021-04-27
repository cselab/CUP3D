#!/bin/bash
NNODE=256
BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}
NU=${NU:-0.00001}
BC=${BC:-freespace}


FACTORY=
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4998336806804398 ypos=2.278319526974906 zpos=1.6032208812671251 bCorrectPosition=true heightProfile=stefan widthProfile=stefan bFixFrameOfRef=1 
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.668371900907971 ypos=1.9159914036420154 zpos=2.3425391289234825 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.274684940836759 ypos=1.7242899937143252 zpos=1.941901373433094 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4758280036509253 ypos=1.7305563175481435 zpos=1.0561693380675288 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.842883615123978 ypos=2.2012198155292206 zpos=1.9988185793484625 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0423669883868887 ypos=2.5296336984141465 zpos=2.5233931875863425 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.051754656451916 ypos=1.0544221257088142 zpos=2.2953799130173667 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3812515925734568 ypos=2.1694907871537628 zpos=2.7350313357086646 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.394460042785757 ypos=1.5063422941544693 zpos=2.1723824136789815 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.292828223760962 ypos=2.237053751975343 zpos=2.592918843323386 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0619756585176026 ypos=2.5078504033481224 zpos=1.543979941638127 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.268492924230659 ypos=2.2162912160472272 zpos=2.2869739264619273 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8526808824118755 ypos=2.2180794732616347 zpos=2.499268257998895 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8568721855067907 ypos=1.8112254474118643 zpos=2.0395643341463114 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.728289091743981 ypos=2.0541137212756686 zpos=1.4954455956573263 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6233713303812314 ypos=2.0895614010248806 zpos=1.6950448322010663 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.702469589804149 ypos=1.226795131616901 zpos=1.5969335793154966 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1017790302388004 ypos=2.2785030313138144 zpos=1.343042270749866 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.230848759306387 ypos=2.0638581172463586 zpos=1.7750061617864588 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6556508540221575 ypos=1.7758189746047082 zpos=2.2334820634050514 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.217256444959015 ypos=2.663330053772636 zpos=2.2694168954477902 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2999344406216338 ypos=1.8882406893124422 zpos=1.750042914990401 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.6983349014767983 ypos=2.171685353907782 zpos=1.727953670741852 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.116819670474937 ypos=1.4754446717861316 zpos=2.2472053222821238 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.817353205101265 ypos=1.4095710779772923 zpos=2.6963117381331347 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1653363564021975 ypos=2.677858241323456 zpos=1.7858039517584694 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.753724622396249 ypos=1.8598796404437792 zpos=1.172362978650692 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.21033120315554 ypos=2.373577834429035 zpos=1.7373577483284552 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.180362196532435 ypos=2.1179180911741042 zpos=2.77623324425373 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5610640133197013 ypos=1.7567501287262481 zpos=1.319908397194802 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.178029543363098 ypos=2.251402301003468 zpos=2.467277323121907 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0969870362970973 ypos=1.7849258406729689 zpos=2.337085750179013 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5352860574177307 ypos=2.2629160599675666 zpos=2.272204824230524 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.253474169571236 ypos=1.5307215746314728 zpos=2.498557533455527 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.2589260114608174 ypos=2.0357360012291803 zpos=1.7782017347997803 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.552845858032502 ypos=1.8042786759142668 zpos=1.4404862624046544 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8751588386508762 ypos=2.135753909541873 zpos=2.075721555359119 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.471106136024425 ypos=2.583333337345614 zpos=1.7398333066505565 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.715252780929776 ypos=2.5242819355900736 zpos=2.153875476870075 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.599582577353193 ypos=1.7910321848587498 zpos=1.7824856109892915 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8522272643728024 ypos=2.5950398829980514 zpos=2.6528723518368347 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.57493534544426 ypos=1.541137102225867 zpos=2.075051623883174 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.325368663740601 ypos=2.6611849788524076 zpos=2.428173980834463 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.023903294039185 ypos=2.080293526230208 zpos=1.5205582037718957 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.2345785732400332 ypos=2.2069261445018284 zpos=2.2437682647156247 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.319968565231447 ypos=1.7687384436374174 zpos=2.2714018114151084 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1436948488302856 ypos=2.834898669257405 zpos=2.207494949619402 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.244086077110725 ypos=2.335731873025119 zpos=2.38423658815441 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.058921695082271 ypos=2.3803373891194757 zpos=2.5713764970925657 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.21927710737027 ypos=1.992164877009577 zpos=1.5598853455374635 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6821131634676805 ypos=1.8104327131546056 zpos=1.168004467026242 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.299662381876884 ypos=2.5848168241096436 zpos=1.6238837549667997 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4357258515108096 ypos=2.1729892428808073 zpos=2.58154273045852 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.523735988567563 ypos=2.2015172140429535 zpos=2.433817385226157 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.6415628616805114 ypos=1.6634690772424054 zpos=1.6782196129615103 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0510834538974367 ypos=1.4247049230110904 zpos=2.042417130498818 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.273365322067716 ypos=1.3424492319738592 zpos=1.530975399560912 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.668699981030996 ypos=2.5784272758445566 zpos=2.1619419804414397 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.600913598370062 ypos=1.1434399632279695 zpos=2.3632864112402125 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.792729267497731 ypos=1.1355726741844512 zpos=2.3404570159502387 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.449883509580639 ypos=1.7097881329628155 zpos=1.9172834370102294 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8340269413389643 ypos=1.4707905885453965 zpos=1.789017241844856 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3123230999692725 ypos=1.7093289477219735 zpos=2.864721260121077 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.519525371337749 ypos=2.0236383002110925 zpos=2.847936597926207 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9313186496232824 ypos=2.2822100985969875 zpos=1.07709810722086 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.430604478477822 ypos=2.65513376046135 zpos=1.6126669858647495 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.8496283342372735 ypos=1.6852285061625918 zpos=2.1057773668613673 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5356802642225893 ypos=2.016446378456024 zpos=1.024425771906046 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.950961412175716 ypos=2.5212921711611473 zpos=1.750742062917155 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7159536361871996 ypos=1.2731143827859455 zpos=2.147861216222803 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.752263049419677 ypos=2.2057942248814024 zpos=2.8850653748142827 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.06106715596064 ypos=2.109276505059264 zpos=1.707162322513741 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2777703891065686 ypos=2.1605892299448346 zpos=1.1741304810739008 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9255185185375634 ypos=2.723057940176084 zpos=1.408708609841815 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1462686883780986 ypos=1.4170713258784846 zpos=1.8158491380039659 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.605100735658626 ypos=2.068238878674115 zpos=2.0408222973014087 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.5782273892460512 ypos=2.3877324523031214 zpos=1.792666406496551 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.280955478348991 ypos=1.737748967640201 zpos=1.2174694354203077 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1388258645714058 ypos=1.426833364401113 zpos=1.1878007444129364 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.517554513386898 ypos=2.5710091916377005 zpos=1.3231879424548754 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.967454044857935 ypos=2.3019829467784954 zpos=1.7729519253207175 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.492004483106189 ypos=1.6549392463732207 zpos=1.5993183220464748 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4212072340556725 ypos=2.6943325587747573 zpos=1.8184309502251714 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4193458604701994 ypos=2.1736800735765156 zpos=1.7528510888388653 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.110086981378488 ypos=1.2557366367445026 zpos=2.646722524153145 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.119778262094086 ypos=2.3115622702279657 zpos=2.3007683904051808 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1872429046906743 ypos=1.99688192157366 zpos=1.750635197126085 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8183630237931678 ypos=1.9716355043713778 zpos=1.5469278291247548 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1863329693630926 ypos=1.5366671523664697 zpos=2.730246917535324 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7056960649008495 ypos=2.0736239025245404 zpos=1.9464741242831414 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.443809896677793 ypos=2.430684862647614 zpos=1.4379400276383731 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0913215308159634 ypos=1.9678702795592524 zpos=2.4783247086364137 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5983749997045633 ypos=2.444179803319517 zpos=1.964206061048419 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8886757632737003 ypos=2.8089762149847024 zpos=1.8097606465434215 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.272773223451885 ypos=2.072840495125703 zpos=1.636496972346396 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.44771209194998 ypos=1.8780847194458494 zpos=1.9023387688817464 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.039708295905398 ypos=1.8877334781049124 zpos=2.3655387633358487 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.443963319483444 ypos=2.167626969067576 zpos=2.125570902439166 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3907986824850793 ypos=1.1990856208753198 zpos=2.0046324539431044 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.30400208355617 ypos=1.7129799075329089 zpos=2.5082529316843876 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9054334689277828 ypos=2.357848869401563 zpos=1.7999208528126291 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4624471654001643 ypos=1.7115104074107572 zpos=1.360395124710009 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.202284134818546 ypos=1.6192262855150363 zpos=2.3763980890165994 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.9885297544564855 ypos=2.461077725549391 zpos=2.395200400746148 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.864253551391505 ypos=2.8666344539640427 zpos=2.3207386814068904 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.3430569462068231 ypos=1.9074677482063076 zpos=2.2971971314419766 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6721546885873746 ypos=1.661165068379948 zpos=2.895062527494007 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.01567614572176 ypos=2.453567174453241 zpos=1.8202926746704264 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.756092339049796 ypos=1.4284116951318988 zpos=1.896812070292552 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.276874996497717 ypos=2.815010851306898 zpos=2.196065562948648 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.885826087021983 ypos=2.492703540213375 zpos=1.553389651230802 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.007242215793533 ypos=2.77532343665123 zpos=2.0226918665096703 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.699243028722093 ypos=1.3655851240894623 zpos=1.8531648636080302 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.519338786870281 ypos=1.8918204737253488 zpos=1.8230703222245765 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.9260847693569387 ypos=2.686687892619769 zpos=1.8290866811340303 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.166753755327241 ypos=2.292189269284017 zpos=1.9814456771950155 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.368511958778617 ypos=2.3318852042888625 zpos=1.9396332027714003 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7523364193527469 ypos=2.5480774495845155 zpos=1.9809487627481877 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5680553451821773 ypos=1.5760286956535765 zpos=2.854145792741229 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.98820366458499 ypos=2.048481595107448 zpos=2.459106209601421 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.601862475332659 ypos=2.047740925222921 zpos=2.629632308458789 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2128716562870325 ypos=1.4843404554185342 zpos=2.2070958841194623 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.228878859974632 ypos=2.3868319201146075 zpos=2.577480803405243 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.2448041073287355 ypos=1.8086047263491052 zpos=1.6930579272117234 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.86446937257542 ypos=2.7235173316446035 zpos=1.8278770901006829 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3745659007463367 ypos=2.208004553238035 zpos=2.7720443564566444 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.13386467323837 ypos=2.1551719797465227 zpos=1.6407530682233384 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.763533193511905 ypos=1.5488420754956573 zpos=1.6643526346405815 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.41107408878664 ypos=2.7855361227890794 zpos=2.346388730109319 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2608118043349226 ypos=2.165329790666642 zpos=1.624872966566412 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0230182303241646 ypos=1.8170775359279479 zpos=1.9365251626077602 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.9380233478094495 ypos=2.5162075355147535 zpos=2.264323311416328 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.821148951584648 ypos=2.6122684323815895 zpos=1.5403009463875437 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.844743160970736 ypos=2.2338966033261136 zpos=1.4650701529382566 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3885126191180377 ypos=2.398024371969836 zpos=2.7262590933119712 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.949046992503095 ypos=2.5058137454557587 zpos=2.4165724781983484 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.274756078725904 ypos=1.426238368838932 zpos=1.5899327179883465 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1892586045437397 ypos=1.662850871484415 zpos=1.8845620489569566 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.522511722703456 ypos=2.1462138509680733 zpos=2.7154649916821056 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4155514407406797 ypos=1.2992080015755807 zpos=2.1208527602947775 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.396472450275808 ypos=1.9162345682856252 zpos=2.081622805405541 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.434876889679504 ypos=1.2075125787024672 zpos=1.5226099446285408 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3836471936716306 ypos=2.3004015634739226 zpos=1.260187192867617 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2964260963179823 ypos=1.8319122680211388 zpos=1.360357166330092 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.457093074735073 ypos=2.574375377300925 zpos=2.2535380046886266 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8047620896895222 ypos=1.0819027749183556 zpos=1.7428970122630436 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.246784766963414 ypos=1.6304904066359254 zpos=1.225251048867007 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.333873766663658 ypos=1.9516854653056788 zpos=1.1241902272793078 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0072190937274796 ypos=1.6665319842131785 zpos=2.271057532923071 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0121647383148336 ypos=2.399845747730954 zpos=1.569741477362387 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1499203321353084 ypos=2.0553730828310557 zpos=2.5144229651052257 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.279837881551856 ypos=1.857349592167485 zpos=1.6578261436545862 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9333011753575025 ypos=1.2122520793145681 zpos=2.490549471652132 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3549525814716397 ypos=2.4533685115116883 zpos=2.2822177531535512 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2737830919006004 ypos=1.4313192898062388 zpos=2.770741311654134 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.317114595602602 ypos=1.6961923997737334 zpos=1.292364220877406 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.4981130725042187 ypos=1.7534517856920642 zpos=2.4955147622634684 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.329809056480803 ypos=2.4562348957829174 zpos=1.7983046221698615 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0020163157456357 ypos=1.9356007853785986 zpos=2.2766412582647737 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8674569491654975 ypos=1.7521797727926753 zpos=1.4484448159581502 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.152959810798527 ypos=2.5777255069372282 zpos=2.5546669318031263 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.972811249771958 ypos=1.4541002192917922 zpos=1.5228885006429755 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.236343458052543 ypos=2.1165038720405587 zpos=2.6846919237453477 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.810802464394659 ypos=2.0778986862651134 zpos=2.190675215202257 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.115749721058285 ypos=2.59434346449319 zpos=2.414757388366997 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7026638586126617 ypos=2.1135284254296507 zpos=1.1880336842459038 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9786624950009997 ypos=1.908223891250455 zpos=1.3214835036130053 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.8006341002551105 ypos=1.6008949972718831 zpos=1.8282912402434697 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2375920452147127 ypos=2.6824828556974314 zpos=2.130584403756862 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.269400967867728 ypos=1.3996266935497599 zpos=1.884595102761665 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8853415187978497 ypos=1.7863057980939518 zpos=1.7144421097427096 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0494124253972683 ypos=1.9506308343368577 zpos=1.6387800204443872 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7778626610290729 ypos=1.3270571621581972 zpos=1.7264334305821014 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.568457469920716 ypos=2.2824840228770293 zpos=2.4465162605684503 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3960544984851966 ypos=2.4191454159633827 zpos=2.5387388957351873 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5729965576005087 ypos=1.4441780488146305 zpos=1.9296308273419098 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.081651266316588 ypos=1.9012625967180516 zpos=2.0586659962246534 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.997400313748872 ypos=2.7014142658916764 zpos=2.2581298383856376 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8091186873884255 ypos=2.110521777586561 zpos=1.0419148785309127 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4131540496798927 ypos=2.3783594926176757 zpos=2.392029329061381 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8252782002956893 ypos=1.937461505748657 zpos=1.3887909690184954 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.551856654244043 ypos=1.4545155219824841 zpos=1.7272337593772238 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8086442012494754 ypos=2.4984230719435376 zpos=1.1924666358756535 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.167392772094953 ypos=1.3068196754968988 zpos=2.37848626265094 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.086151822964976 ypos=2.4539621245567234 zpos=2.6611916561601148 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.657620237631664 ypos=2.438303182678486 zpos=2.0955081420889523 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5673893743684295 ypos=2.5501359351013395 zpos=2.4014043438762513 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.142484519772489 ypos=1.6957694013603062 zpos=1.5129241205019008 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.609888816044984 ypos=2.0725474055407203 zpos=2.8087714393263057 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.716799734770775 ypos=2.4640747949127295 zpos=2.64634034455356 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.675697947101896 ypos=2.8347449665695397 zpos=2.3535382817643034 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4642852796093804 ypos=1.4117862566964392 zpos=2.094486548253918 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.219901769861451 ypos=1.380670031901504 zpos=1.5552438574678311 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7862617812141335 ypos=1.8576629889836302 zpos=2.2931769545903604 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.31722859447314 ypos=2.3802367102553403 zpos=2.218149878629514 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.831431418981922 ypos=2.339350555881618 zpos=2.4521259870478778 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.322233144654412 ypos=1.1169787331440402 zpos=2.1166001351863732 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8933560924922643 ypos=1.748142948860211 zpos=2.1041874916462584 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.354153121534256 ypos=2.082723702947697 zpos=1.4585113390719888 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.259267297791583 ypos=1.6329639989056033 zpos=1.6709374674095319 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.575434357932995 ypos=2.4451450616109085 zpos=2.125042501968908 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.845114503502358 ypos=1.7804160429316 zpos=2.329014474886266 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6387074640461385 ypos=2.183495443160228 zpos=1.393700297608666 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.761767943488328 ypos=2.2588189881872016 zpos=1.6141584839017256 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7422087182430832 ypos=1.5512581845968036 zpos=1.6113701480185338 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9928642455008414 ypos=2.6204263296759622 zpos=1.4183576585301751 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5262383926381666 ypos=1.3471242506968282 zpos=1.5616965026192309 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9359108096620528 ypos=2.4211850271249267 zpos=1.3823682853765247 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3065519342442795 ypos=1.762321780605792 zpos=2.232949212282957 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.4990896312573212 ypos=2.181596200294715 zpos=2.263982478022876 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4728776757255795 ypos=2.5766269334420637 zpos=2.29502634457073 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.900704905435064 ypos=2.0452457447860892 zpos=1.8886686849686454 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5965157980989004 ypos=1.2935988741706712 zpos=1.6244317768916634 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3077751928358445 ypos=1.946195943011821 zpos=2.118954908505841 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2509912735720303 ypos=1.5051948207827333 zpos=1.4560725398474519 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9325616046251906 ypos=1.4608916883142302 zpos=1.8617350098222745 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2033382929328638 ypos=2.6013625822791506 zpos=2.776412637227035 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.248523973690452 ypos=1.9490579681417877 zpos=1.7328847859426384 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9596387053057014 ypos=2.061309808629678 zpos=1.7518182158718603 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.894311372347302 ypos=1.3970178619241334 zpos=1.9826098953440134 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8193549820048984 ypos=2.1966151420196076 zpos=2.1164292602463335 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9807821826842344 ypos=2.871981082519123 zpos=1.5227901740512995 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.376766571707408 ypos=2.0910254762480793 zpos=1.6610484177843803 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9116692167621756 ypos=1.3012667225924432 zpos=1.5456705317334338 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.36775462472975 ypos=1.844586191075277 zpos=1.5548166664717256 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.772393332555151 ypos=1.8534789435030756 zpos=1.7918088311284512 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.58443454279454 ypos=1.521917042054429 zpos=1.9292400004070454 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.14237959089388 ypos=1.5953043865334542 zpos=2.0990174081954645 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.721652213544145 ypos=2.2280529913451748 zpos=2.2221120460344745 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5894188729452132 ypos=2.762129648763022 zpos=1.5401734728581844 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4735977641958207 ypos=1.5868342705663006 zpos=1.1932591778983128 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6844350136828083 ypos=2.29710431533032 zpos=1.1002359002100406 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.4684245338206379 ypos=2.444321683473783 zpos=2.3665748648753224 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.28842884985545 ypos=1.9082860883855557 zpos=2.0151162380557577 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8293039303364016 ypos=1.7282639323681133 zpos=1.7674062761746683 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7232289445440334 ypos=1.8471293527648955 zpos=1.8814185425284 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.5824465990972705 ypos=2.4234629669607557 zpos=2.3275501125586375 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.748970263239857 ypos=2.1472236844593424 zpos=1.8577826643853121 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.888387969111337 ypos=2.200130212430669 zpos=2.083277573464229 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7311864007833853 ypos=2.17438282023931 zpos=2.4530233832450317 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.5203608828484825 ypos=1.6389364076283686 zpos=2.1588812250093694 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.726104577353406 ypos=2.7944133741236943 zpos=1.6549424010871803 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.3901601895714695 ypos=1.7282680002836248 zpos=1.7820430028401764 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8202613194995765 ypos=1.771556287987078 zpos=2.038566338816312 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.762381450301381 ypos=2.583061219650894 zpos=2.008411557086495 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.630739790872195 ypos=2.8850136483425435 zpos=2.014610653167594 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.408242446977917 ypos=2.2242555329497287 zpos=1.8225852884159177 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.4470681427472973 ypos=1.1764672162637821 zpos=2.1463523233657362 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.908876249455508 ypos=2.5071942296231593 zpos=2.4417300574342606 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.401478437326839 ypos=2.151926019341845 zpos=2.5446398609565413 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2445083860952066 ypos=1.4215915058694382 zpos=2.623522078573861 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.884118781946436 ypos=1.8866627283939743 zpos=2.656978524578058 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.579966830445521 ypos=2.4184890760207196 zpos=2.0399466899141574 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.0414587136959055 ypos=1.7571858725569467 zpos=2.6616942174844485 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8351352617327652 ypos=1.2967692566135396 zpos=2.18277632351082 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.795554561732585 ypos=2.1491142499034965 zpos=1.1920342219968802 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.436212147048454 ypos=2.0503841900638315 zpos=2.8059573106307147 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2417302814770443 ypos=1.5390609696511226 zpos=2.510942306775058 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.5934381910961317 ypos=1.7560722660065742 zpos=1.9632446749473826 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9728049069415947 ypos=1.5443904937387352 zpos=1.1369318231002836 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.139811770129563 ypos=2.831526638632641 zpos=2.065647709359567 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.905404863892174 ypos=2.01168271132488 zpos=1.2086382491128889 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.935463532638556 ypos=1.8053102821034916 zpos=1.9880094177650314 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.568110399849287 ypos=1.1481679575851524 zpos=1.9798427834446668 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.978251910836785 ypos=2.0938364607259636 zpos=1.889934268812372 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.204963839031143 ypos=1.2698945861548894 zpos=2.1412056518238853 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.060150153332719 ypos=1.9103864668853108 zpos=2.175125881401475 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.058426449517197 ypos=2.0650507848459845 zpos=2.520466587935968 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.4296863301901608 ypos=2.149840909268376 zpos=1.9304172468583516 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.765044724018578 ypos=1.3645792149194254 zpos=2.52319555055059 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.787128080844882 ypos=2.3002403672032536 zpos=1.377270128535784 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.625527640986748 ypos=2.165115749245423 zpos=2.5638265347544795 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.2549736540541576 ypos=2.7428563849529937 zpos=1.4749715509657504 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.8990766158686294 ypos=2.4729045019283076 zpos=1.3328712601108716 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.109018347626194 ypos=1.651317758238394 zpos=2.4486387587073977 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5738615359811767 ypos=1.3802507702787052 zpos=2.399795685281064 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.456381255757271 ypos=2.7008429367279114 zpos=2.116139563508183 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.517000595914172 ypos=2.1729504170933143 zpos=1.2922285457690532 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7258060343588866 ypos=2.423192933088254 zpos=1.5324290817617734 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.470550611662295 ypos=2.3582743464876583 zpos=2.1569970816197226 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4116769008065346 ypos=1.3404793186286297 zpos=1.7153378277142748 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.133425903366353 ypos=1.9599464643037228 zpos=1.394083600678978 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.567555351006176 ypos=1.7171053912596053 zpos=1.6140196352322622 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.68785156997754 ypos=1.8755313034913919 zpos=2.440934738773504 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8967269935419733 ypos=1.6632225813349282 zpos=1.9870035998920494 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.118755834728141 ypos=2.7249335732742086 zpos=1.9239539789251499 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7853325548950667 ypos=1.7179265573193285 zpos=1.8619694587668445 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.024676427601245 ypos=1.303223580873707 zpos=1.7130001200641027 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7465805059874504 ypos=2.34096165482984 zpos=2.069173892203801 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.0354595068707075 ypos=1.1620998613202724 zpos=2.052402634411606 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.285964353432576 ypos=2.042550207364309 zpos=2.8441314063502032 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.7371200303579393 ypos=1.5361476168434791 zpos=1.5241901972806784 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.852909250212452 ypos=2.3084531752647304 zpos=1.7224340138487182 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.6629603042257977 ypos=2.5742806839276104 zpos=1.7950833779785569 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.448457618196081 ypos=2.311118165937984 zpos=1.497875460805914 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4735971619378767 ypos=2.0978304611842344 zpos=1.0603474463159521 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9910051824252317 ypos=1.1753071368937666 zpos=2.088199692214976 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.577258733086007 ypos=2.603537621352059 zpos=2.6236632115869143 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.709563676174177 ypos=2.7978845575110727 zpos=1.4665921025375486 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9939183511912266 ypos=1.6502327008084523 zpos=1.0662711450472555 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.4637159029752116 ypos=1.4535040249652522 zpos=1.3219019598009574 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.452317335869621 ypos=1.6068499067993054 zpos=1.4620549543545622 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.9616904206970625 ypos=2.3187023823813258 zpos=1.6403175418357288 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.247701291243266 ypos=2.1659331938961603 zpos=1.7283281876825924 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.723934479543843 ypos=1.594889504378826 zpos=1.157093213498735 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.559734462139835 ypos=2.2999427913573736 zpos=2.3239912265996794 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2412296079589313 ypos=2.8792994135124768 zpos=2.108605951863569 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7649596163436847 ypos=2.0234589024464578 zpos=2.417811073976173 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.128296418347692 ypos=2.151379891065581 zpos=1.2223925794022494 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.91384486293397 ypos=1.4853551982905266 zpos=1.4363600101343637 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8653953486280508 ypos=2.6608062229360683 zpos=1.9869321512908416 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.7027774869587726 ypos=2.7474207047899384 zpos=1.6110369839819336 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3611374491193553 ypos=1.6270784789543744 zpos=2.796746538766002 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.3066185665683 ypos=2.3188523509472008 zpos=2.3919372022282097 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.80731588808621 ypos=1.9332474106429056 zpos=2.5066880583601234 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.227982184600409 ypos=1.7636280586182564 zpos=2.7701518101553186 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7441900198054325 ypos=2.0729484879162383 zpos=1.8170052420727831 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.095944440597186 ypos=2.1118773879962776 zpos=1.374999710645183 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.040552469049749 ypos=2.8805406097946253 zpos=2.1467837920789545 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.127203291037703 ypos=1.4505728854927644 zpos=2.1775450884036163 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.656426708739772 ypos=2.146788861666193 zpos=1.9301048511747227 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.089436398776035 ypos=2.193733104722713 zpos=1.073742731172829 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.7590455546866437 ypos=2.7978036648287996 zpos=1.9698773265158465 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.73554605136651 ypos=1.99904104326972 zpos=1.0169156090420897 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.343346079135967 ypos=1.2643807675397558 zpos=1.7909204512420327 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.543959692936298 ypos=2.1683182695198275 zpos=2.3054082795201416 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.2895266297484134 ypos=2.0403049686511987 zpos=1.5538247372623617 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6813420301606445 ypos=1.1860031362078791 zpos=2.3916575962813287 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.1005146664665153 ypos=2.150585037879641 zpos=1.51109207483716 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.0825030926412103 ypos=1.4556593086721568 zpos=2.1169008801979654 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.1883837044791115 ypos=1.7595473198831466 zpos=1.1653395122298456 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.5127046748236466 ypos=2.2331156656861113 zpos=1.8558844127628706 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.371453568562759 ypos=1.3712993208714317 zpos=2.377255179926777 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.844563167183272 ypos=1.4160211660075093 zpos=2.2203897245826827 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.9482125971859166 ypos=2.0721504902633527 zpos=2.3565146444749563 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.3749864983748976 ypos=2.5178651445430527 zpos=1.6208243016158885 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.6509248349954926 ypos=1.894499625785988 zpos=1.3785571288876226 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.806211157848615 ypos=1.9309900990794038 zpos=1.948422207580536 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.417553873816779 ypos=2.0297140282147224 zpos=2.0693934033368815 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.6839322944463158 ypos=2.4663195765178743 zpos=2.298253902200292 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.8440873897396584 ypos=1.9711382405228497 zpos=2.9307299826387796 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.683619361315044 ypos=2.55699215959252 zpos=1.926143475036214 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=3.9010924531978692 ypos=1.5614181710620894 zpos=2.7602312761977674 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.1976682597405406 ypos=2.1154005910420874 zpos=1.916281861111909 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=4.474486991005364 ypos=2.5793764099545418 zpos=1.5981973166975936 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.9250171531565754 ypos=1.8349679349770303 zpos=2.487023377474567 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.601888691884691 ypos=1.705259748541374 zpos=1.4037071055725956 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=2.737777541961435 ypos=1.3600238547590662 zpos=2.192669136438587 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.8376365328519868 ypos=1.9090568564795876 zpos=1.7500211764557412 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=1.0 xpos=1.946100300634783 ypos=2.037665271174951 zpos=1.6589579051652152 bCorrectPosition=true heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 8.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 100.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.4 -use-dlm -1 -nu ${NU}"
OPTIONS+=" -levelMax 7 -levelStart 4 -Rtol 0.1 -Ctol 0.01"
OPTIONS+=" -TimeOrder 2"
OPTIONS+=" -poissonTol 1e-5 "
OPTIONS+=" -poissonTolRel 1e-3 "
