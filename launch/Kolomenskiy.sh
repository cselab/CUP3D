cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=double particles=false rk2=true -j
cd ../launch

for L in 1e3 1e6
do
    for C in .01 .05 .1
    do
        for BPD in 8 16 32
        do
            NAME=Kolomenskiy_2010_Lambda${L}_CFL${C}_bpd${BPD}
            FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
            mkdir ${FOLDER}
            cp ../makefiles/simulation ${FOLDER}
            export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .0125 -tend 10. -rhoS 2. -ypos .9 -nu 0.000155188776172763 -sim falling -tdump 0.01 -lambda ${L}
        done
    done
done