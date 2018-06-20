#!/bin/bash

SETTINGSNAME=$1
BASEPATH=${SCRATCH}/CubismUP3D/

if [ ! -f $SETTINGSNAME ]; then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/Swim3D_open ${FOLDER}

cd $FOLDER

export OMP_NUM_THREADS=12
echo "$OPTIONS" > settings.txt
mpirun -np 1 ./Swim3D_open ${OPTIONS} -factory-content "${FACTORY}"
