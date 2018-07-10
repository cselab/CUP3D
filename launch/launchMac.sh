#!/bin/bash

SETTINGSNAME=$1
BASEPATH=../runs/

if [ ! -f $SETTINGSNAME ]; then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/simulation ${FOLDER}

cd $FOLDER

export OMP_NUM_THREADS=4
echo "$OPTIONS" > settings.txt
mpirun -np 1 ./simulation ${OPTIONS} -factory-content "${FACTORY}"
