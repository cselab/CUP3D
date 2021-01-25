#!/bin/bash
NNODE=64
BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}
NU=${NU:-0.00004}
BC=${BC:-freespace}


FACTORY=
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.241725980659904 ypos=1.725974108505949 zpos=1.69867559906042 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.635359094863067 ypos=1.7872626631602615 zpos=2.350613155963831 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.108823252494633 ypos=2.1024812015377936 zpos=2.2089761844101945 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.46126867026646 ypos=2.028978533516651 zpos=1.8579230219281577 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.1274473360417945 ypos=2.1665093825302284 zpos=2.2052466551307766 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8737390991441143 ypos=2.3392330301528967 zpos=1.9887876722457019 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5313296639802205 ypos=2.0860241263082693 zpos=2.0581088975251136 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.5743166791278504 ypos=2.1798273908851944 zpos=1.7132105435485518 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.723592388026626 ypos=1.8826875483940841 zpos=2.025166557343849 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.802121910016976 ypos=2.3710499755539742 zpos=1.9650638711110846 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.2961030341187625 ypos=2.190982579001435 zpos=2.2425369555773296 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.957781603049489 ypos=1.8512158103681153 zpos=2.2599435699254466 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.911091531291747 ypos=1.9696356850456853 zpos=1.5601849557454468 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6940194518317435 ypos=1.6352974499793969 zpos=1.9311262097279103 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.557066900611574 ypos=2.1178582205698566 zpos=1.5658243932230627 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.307543199247751 ypos=1.6651552527427183 zpos=2.2522272616756713 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.28826980776311 ypos=2.1160777364976266 zpos=2.39465615985996 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6677954226396894 ypos=2.3167814759222494 zpos=2.0733321866959185 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.980583919138153 ypos=1.608901884088163 zpos=1.9647361548095377 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.062594011020508 ypos=2.1682143609206945 zpos=1.7015605253553507 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.493434171932855 ypos=2.1052380726402906 zpos=2.053514562518789 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9017095399110655 ypos=2.0818428884902 zpos=2.1229276773514787 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.0082336392841915 ypos=2.1609635898042807 zpos=1.6189240350497789 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.543528701468778 ypos=1.8231264903510875 zpos=2.435280162161764 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9018974289705337 ypos=1.988224714148739 zpos=1.9283585326446742 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6553646620857356 ypos=1.9023442629520404 zpos=1.7645221189351517 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.454088875997865 ypos=2.423453366530218 zpos=1.960995379952867 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.723068536235346 ypos=1.8141247906764129 zpos=1.837680640225766 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.924063991610564 ypos=2.292315732999405 zpos=2.1316215921557657 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3781718577207003 ypos=1.990675209671464 zpos=2.0826494526548602 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.351124486853481 ypos=2.4060377866475715 zpos=1.9811090041305166 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.490353032509347 ypos=2.239699338044213 zpos=1.8257321806472189 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0185015307086926 ypos=2.319981828542893 zpos=1.911965924965441 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.413652522081371 ypos=1.99457966892407 zpos=1.7038773830645395 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.385355019200582 ypos=2.0614160360355642 zpos=2.084811118464794 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.478828819812437 ypos=1.813500792109723 zpos=1.9506687694025033 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.000340029923915 ypos=1.968635308098656 zpos=1.7471193072887385 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.598333971901903 ypos=1.9414443294580896 zpos=1.6213925949800072 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.088038496348445 ypos=1.999906600800906 zpos=1.8565354835415557 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.375116887803763 ypos=1.9996221544305683 zpos=2.4234930616060826 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.611272404115417 ypos=2.0699027778617967 zpos=1.67512890551641 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.63135729632774 ypos=1.7876279507373836 zpos=2.2138639286187485 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.519769349607002 ypos=1.6760406807387067 zpos=1.7345305401831652 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.520942183583816 ypos=2.2079379449623304 zpos=2.0349300439411695 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.547787292835658 ypos=2.1233990209154165 zpos=2.3221121030696237 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.954882194153366 ypos=1.9943313657084318 zpos=1.9029425884339357 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1546967812416624 ypos=2.0829748640733188 zpos=1.7354935798444844 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.145785628155004 ypos=1.9344744161804386 zpos=2.083874155629174 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.166029491914951 ypos=2.0035234520011955 zpos=1.781033640273639 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.557801216949854 ypos=1.8333375768681848 zpos=2.2025077217939804 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.545052046180199 ypos=1.6438149480312654 zpos=2.250241997925284 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.365313403588843 ypos=1.9529983222109042 zpos=1.9247008059569726 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1134939492563944 ypos=2.0783936674621817 zpos=2.302160567463001 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.726448371819282 ypos=2.341513341668841 zpos=1.7238188354435242 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9903246480373045 ypos=2.165436575214651 zpos=2.2630137453476986 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.468185315232396 ypos=2.021908246059358 zpos=2.065880429452986 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.081773314549496 ypos=2.0836962701589834 zpos=1.7504574219331361 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.490722158486208 ypos=2.0868455645863624 zpos=2.25934392074831 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5801036676570086 ypos=1.7974198043540872 zpos=2.0966609730715784 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.6773626444235905 ypos=2.1018291784463523 zpos=1.6462076616469627 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.845677033742439 ypos=2.074526329325121 zpos=2.2209031376987842 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.029868058255481 ypos=2.1188769035171005 zpos=1.9053181683019944 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.290717107688694 ypos=1.8825247443958233 zpos=1.884090183483202 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.035015597945393 ypos=1.7415436065418062 zpos=1.7652672204344666 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.22318633350373 ypos=2.2681074312182865 zpos=1.9027914945339066 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.1656436023733106 ypos=1.983445079693648 zpos=1.5425139702712158 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.765636540567593 ypos=2.3025624125275135 zpos=1.9738539947993576 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9714756500760715 ypos=1.7424447804805054 zpos=2.0754550253006525 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6564209603277917 ypos=2.2656566007091636 zpos=2.087935249921818 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.1123965123228943 ypos=2.2124900020939267 zpos=1.713634692494191 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.0384153897205155 ypos=2.1790980652770373 zpos=2.137799272347858 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.592615882814392 ypos=2.197824437128722 zpos=2.078725692813327 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.8539953619187255 ypos=1.7370575764607095 zpos=1.76279944240112 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.806957630696776 ypos=2.3502282424300147 zpos=2.20896929767705 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.3676335146704584 ypos=1.8311228513812 zpos=2.2140985482179003 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.960626273230546 ypos=2.0811697102817397 zpos=2.32712703361273 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.3006995764448925 ypos=2.2038322527085223 zpos=1.769518932322758 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.8491054438001386 ypos=2.2039415818767396 zpos=2.1668154020280372 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.4563929632261665 ypos=2.2429758324801807 zpos=2.2892635378149584 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.445458132668418 ypos=1.6902418284981469 zpos=2.0320359573884894 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.418292899635475 ypos=1.6383883706051126 zpos=1.904381016434543 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.898667300977612 ypos=2.150286292363303 zpos=2.0250698546481787 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.496228044850708 ypos=2.113968144457923 zpos=1.8889705918367212 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.219342620714523 ypos=2.1251473392760456 zpos=2.3567303568877196 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.821696769462076 ypos=2.07927214240334 zpos=2.0681013195680933 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.126640493789131 ypos=2.1849488032384103 zpos=2.220574748018017 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.353179210826403 ypos=2.159278187856182 zpos=1.9521884749528904 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.01077787088101 ypos=1.9561615712979714 zpos=1.6594206470044852 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.9608313461801936 ypos=2.404042978173958 zpos=2.043247998872221 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.069449673448678 ypos=1.8309270490828196 zpos=1.8523237531904067 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.067564094506036 ypos=1.6553121215229387 zpos=2.011578020065981 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.2503296767341854 ypos=1.731899890755324 zpos=2.0249209452481756 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.2836975673783675 ypos=1.9689481671229463 zpos=1.9884067061390283 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.660220452282417 ypos=2.335133890279086 zpos=2.3111602895043313 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.691010856241055 ypos=2.3974296525499614 zpos=2.1072773266737355 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.6415756420009417 ypos=1.7494598747824788 zpos=1.8255502397148125 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.660765228771405 ypos=1.5621556967840506 zpos=2.1135634043544185 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.9353703353073786 ypos=1.95867621220148 zpos=2.2239205880556754 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.9546346875003 ypos=1.846833059752618 zpos=1.6556671775005547 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.404453363019124 ypos=2.3998654811509743 zpos=2.1180793027642832 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.23333951438358 ypos=1.7182813848433625 zpos=1.972655478325574 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.170184524584811 ypos=1.5166737987466536 zpos=2.075200177889151 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.801107893903021 ypos=1.8681207769994574 zpos=1.7483209538044964 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5511500951549992 ypos=1.8478434915029616 zpos=2.4361532002274497 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4731009856926653 ypos=1.8045807637520221 zpos=2.291392764853974 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.023787435247773 ypos=1.7959446023979935 zpos=2.319056978628826 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.8878914154104525 ypos=1.8286848886763352 zpos=1.937545849513553 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.890777169217452 ypos=1.8950397671139314 zpos=2.10596462037434 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.471803233457407 ypos=1.8101973178016115 zpos=1.6521265447216344 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.5350543672911834 ypos=2.307548590536696 zpos=1.996176297975898 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.348383750290531 ypos=2.1134926836669234 zpos=1.9300485017880877 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8563940618331083 ypos=1.8333393785763277 zpos=2.0061716254221618 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.8546924338557993 ypos=2.1764632379703674 zpos=1.9385236616582584 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.7981579299563077 ypos=1.7901812462401638 zpos=1.8543535735605083 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.129852365028069 ypos=2.2186354571727027 zpos=2.306507565100522 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.736769793947149 ypos=2.0786394451890104 zpos=2.042194074379593 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.162039835397617 ypos=1.5105578284742833 zpos=1.9617493941579236 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.4916606394750382 ypos=2.3732380131633604 zpos=2.2271883939779897 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=5.0056625100498 ypos=1.8789904391767256 zpos=2.190787387524483 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.947714145010302 ypos=2.078231205308241 zpos=1.835360131680563 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=3.644615607813337 ypos=2.329514414930955 zpos=2.2922066982413916 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.7638874454511155 ypos=2.067460820792921 zpos=2.159773431115502 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.440542762125146 ypos=2.2554341914228737 zpos=2.008142625164107 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.118036233300574 ypos=1.8118849697148263 zpos=2.449388622738147 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.773570785925521 ypos=1.7413541997896849 zpos=2.2689789589281126 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.105250702032026 ypos=2.0500760597600953 zpos=2.4506880749743822 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=4.03243448185125 ypos=1.9640350892206726 zpos=1.6858109542934256 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"
FACTORY+="CarlingFish L=0.2 T=1.0 xpos=2.761070401964148 ypos=2.1972494622709853 zpos=1.943923535831162 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -extentx 8.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 100.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.80 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 7 -levelStart 3 -Rtol 0.1 -Ctol 0.01"
OPTIONS+=" -Advection3rdOrder=true"
