#!/bin/bash
#module load gcc
NNODE=$1
NTHREADSPERNODE=48
export OMP_NUM_THREADS=$NTHREADSPERNODE
export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
export PATH=/cluster/home/mavt/chatzidp/usr/mpich3/bin/:$PATH
export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/fftw-3.3.5/lib/:$LD_LIBRARY_PATH
#LD_LIBRARY_PATH=/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.9.2/lib/:$LD_LIBRARY_PATH

CFL=0.1
LAMBDA=1e4
BPDX=40
BPDY=32

OPTIONS=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDY} -nprocsx ${NNODE} -CFL ${CFL} -length 0.1 -uinfx 0.1 -lambda ${LAMBDA} -nu 0.00005"

sort $LSB_DJOB_HOSTFILE | uniq  > lsf_hostfile
mpich_run -np ${NNODE} -ppn 1 -bind-to none -launcher ssh -f lsf_hostfile ./simulation -tend 8 ${OPTIONS}
