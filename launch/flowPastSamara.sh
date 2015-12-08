#module load gcc
export OMP_NUM_THREADS=48

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/
BASENAME_FIXED=Drag_FlowPastFixedSphere_2811_
BASENAME_MOVING=Drag_FlowPastMovingSphere_2811_

cd ../makefiles/;make clean;make config=production bc=periodic precision=single particles=false dlm=true -j;cd ../launch/

for CFL in 0.01
do
for LAMBDA in 1
#1e4 1e5 1e6
do
for B in 4 8 16
do
	OPTIONS=" -bpdx ${B} -bpdy ${B} -bpdz ${B} -CFL ${CFL} -radius .025 -uinf .01 -rhoS 1.00 -lambda ${LAMBDA} -dlm ${LAMBDA}"

	# Re 40
#	NAME_RE40=Re40_DLM${LAMBDA}_CFL${CFL}_$[${B}*32]

#	FOLDER_FIXED=${BASEPATH}${BASENAME_FIXED}${NAME_RE40}
#	mkdir ${FOLDER_FIXED}
#	cp ../makefiles/simulation ${FOLDER_FIXED}
#	cp launchValidationDrag.sh ${FOLDER_FIXED}
#	bsub -n 48 -W 10:00 -o ${BASENAME_FIXED}${NAME_RE40} -J ${BASENAME_FIXED}${NAME_RE40} ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/${BASENAME_FIXED}${NAME_RE40} -sim fixed -Re 40 -tend 8 ${OPTIONS}


#	FOLDER_MOVING=${BASEPATH}${BASENAME_MOVING}${NAME_RE40}
#	mkdir ${FOLDER_MOVING}
#	cp ../makefiles/simulation ${FOLDER_MOVING}
#	cp launchValidationDrag.sh ${FOLDER_MOVING}
#	bsub -n 48 -W 10:00 -o ${BASENAME_MOVING}${NAME_RE40} -J ${BASENAME_MOVING}${NAME_RE40} ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/${BASENAME_MOVING}${NAME_RE40} -sim moving -Re 40 -tend 8 ${OPTIONS}



	# Re 100
	NAME_RE100=Re100_DLM${LAMBDA}_CFL${CFL}_$[${B}*32]

	FOLDER_FIXED=${BASEPATH}${BASENAME_FIXED}${NAME_RE100}
	mkdir ${FOLDER_FIXED}
	cp ../makefiles/simulation ${FOLDER_FIXED}
	cp launchValidationDrag.sh ${FOLDER_FIXED}
	bsub -n 48 -W 40:00 -o ${BASENAME_FIXED}${NAME_RE100} -J ${BASENAME_FIXED}${NAME_RE100} ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/${BASENAME_FIXED}${NAME_RE100} -sim fixed -Re 100 -tend 1000 ${OPTIONS}


	FOLDER_MOVING=${BASEPATH}${BASENAME_MOVING}${NAME_RE100}
	mkdir ${FOLDER_MOVING}
	cp ../makefiles/simulation ${FOLDER_MOVING}
	cp launchValidationDrag.sh ${FOLDER_MOVING}
	bsub -n 48 -W 40:00 -o ${BASENAME_MOVING}${NAME_RE100} -J ${BASENAME_MOVING}${NAME_RE100} ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/${BASENAME_MOVING}${NAME_RE100} -sim moving -Re 100 -tend 1000 ${OPTIONS}



	# Re 1000
#	NAME_RE1000=Re1000_DLM${LAMBDA}_CFL${CFL}_$[${B}*32]

#	FOLDER_FIXED=${BASEPATH}${BASENAME_FIXED}${NAME_RE1000}
#	mkdir ${FOLDER_FIXED}
#	cp ../makefiles/simulation ${FOLDER_FIXED}
#	cp launchValidationDrag.sh ${FOLDER_FIXED}
#	bsub -n 48 -W 10:00 -o ${BASENAME_FIXED}${NAME_RE1000} -J ${BASENAME_FIXED}${NAME_RE1000} ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/${BASENAME_FIXED}${NAME_RE1000} -sim fixed -Re 1000 -tend 5 ${OPTIONS}


#	FOLDER_MOVING=${BASEPATH}${BASENAME_MOVING}${NAME_RE1000}
#	mkdir ${FOLDER_MOVING}
#	cp ../makefiles/simulation ${FOLDER_MOVING}
#	cp launchValidationDrag.sh ${FOLDER_MOVING}
#	bsub -n 48 -W 10:00 -o ${BASENAME_MOVING}${NAME_RE1000} -J ${BASENAME_MOVING}${NAME_RE1000} ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/${BASENAME_MOVING}${NAME_RE1000} -sim moving -Re 1000 -tend 5 ${OPTIONS}
done
done
done