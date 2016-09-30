#!/bin/bash
SETTINGSNAME=$1
BASEPATH=/cluster/scratch/novatig/CubismUP3D/

if [ $# -lt 2 ] ; then
    WCLOCK=24:00
    echo "Setting WCLOCK to ${WCLOCK}"
else
    WCLOCK=$2
fi

INTERACTIVE=0
if [ $# -gt 2 ] ; then
    if [ "${3}" = "node" ]; then
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
cp launchBrutus.sh ${FOLDER}
cp runBrutus.sh ${FOLDER}/run.sh

#CURRDIR=`pwd`
cd $FOLDER

if [ $INTERACTIVE -eq 1 ] ; then 
   export OMP_NUM_THREADS=24
   export MV2_ENABLE_AFFINITY=0 
   echo $OPTIONS > settings.txt
   mpirun -np ${NNODE} -ppn 1 ./simulation ${OPTIONS}
else
    bsub -n ${NPROCESSORS} -W ${WCLOCK} -J ${BASENAME} < run.sh 
fi

#cd $CURRDIR
#valgrind --tool=memcheck --track-origins=yes --leak-check=yes
