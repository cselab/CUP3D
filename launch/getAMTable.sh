for B in 4 8 16 32 64
do
	for R in 0.01 0.1 0.5 0.9 1.0 1.1 2.0 10.0 100.0 1000.0 10000.0
	do
		cat Disk_rhoS${R}_addedmass.dat | awk -v nb=${B} '{ if($3==0 && $1==nb && $2!=1) printf "%f\t",-($4-$5)/$5 }' >> fulltable
	done
	echo '' >> fulltable
done
