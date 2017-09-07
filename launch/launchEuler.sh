#!/bin/bash
SETTINGSNAME=$1
BASEPATH=/cluster/scratch/novatig/CubismUP3D/

INTERACTIVE=0
if [ $# -gt 1 ] ; then
    if [ "${2}" = "node" ]; then
        echo "Running on current node"
        INTERACTIVE=1
    fi
fi

if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

NPROCESSORS=$((${NNODE}*24))
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
cp ${FFACTORY} ${FOLDER}/factory
cp ../makefiles/simulation ${FOLDER}
#cp launchBrutus.sh ${FOLDER}
cp runEuler.sh ${FOLDER}/run.sh

#CURRDIR=`pwd`
cd $FOLDER

if [ $INTERACTIVE -eq 1 ] ; then 
   export OMP_NUM_THREADS=48
   export MV2_ENABLE_AFFINITY=0 
   echo $OPTIONS > settings.txt
   mpirun -np ${NNODE} -ppn 1 ./simulation ${OPTIONS}
else
    bsub -n ${NPROCESSORS} -W 24:00 -J ${BASENAME} < run.sh 
fi

#cd $CURRDIR
#valgrind --tool=memcheck --track-origins=yes --leak-check=yes
