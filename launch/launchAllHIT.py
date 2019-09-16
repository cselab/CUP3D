#!/usr/bin/env python3
import os, numpy as np, argparse


def    relFit(nu, eps): return 7.10538 * np.power(eps, 1/6.0) / np.sqrt(nu)

def    etaFit(nu, eps): return np.power(eps, 0.25) * np.power(nu, 0.75)

def lambdaFit(nu, eps): return 5.2623 * np.power(eps,-1/6.0) * np.sqrt(nu);

def runspec(nu, eps, run):
    return "HIT_DNS2_EXT2pi_EPS%.03f_NU%.04f_RUN%d" % (eps, nu, run)

def getSettings(nu, eps):
    return '-bpdx 12 -bpdy 12 -bpdz 12 -extentx 6.2831853072 -CFL 0.02 ' \
           '-dump2D 0 -dump3D 0 -tdump 1 -BC_x periodic -BC_y periodic ' \
           '-BC_z periodic -initCond HITurbulence -spectralForcing 1 ' \
           '-compute-dissipation 1 -nprocsx 1 -nprocsy 1 -nprocsz 1 ' \
           '-spectralIC fromFit -analysis HIT -tAnalysis 0.1 -tend 100 ' \
           '-keepMomentumConstant 1 -nu %f -energyInjectionRate %f ' % (nu, eps)

def launchEuler(nu, eps, run):
    runname = runspec(nu, eps, run)
    print(runname)
    os.system( "export NU=%f  \n \
                export EPS=%f \n \
                echo $NU $EPS \n \
                ./launchEuler.sh settingsHIT_DNS.sh %s " \
                % (nu, eps, runname) )

def launchDaint(nCases):
    f = open('HIT_sbatch','w+')
    f.write('#!/bin/bash -l \n')
    f.write('COMMNAME=DNS_HIT \n')
    f.write('#SBATCH --job-name=${COMMNAME} --time=24:00:00 \n')
    f.write('#SBATCH --output=out.txt --error=err.txt --constraint=gpu \n')
    f.write('#SBATCH --account=929 --partition=normal \n')
    f.write('#SBATCH --array=0-%d --ntasks-per-node=1 \n' % (nCases-1) )
    
    f.write('ind=$SLURM_ARRAY_TASK_ID \n')
    f.write('RUNDIRN=`./launchAllHIT.py --case ${ind} --printName` \n')
    f.write('OPTIONS=`./launchAllHIT.py --case ${ind} --printOptions` \n')

    f.write('mkdir -p ${SCRATCH}/CubismUP3D/${RUNDIRN} \n')
    f.write('cd ${SCRATCH}/CubismUP3D/${RUNDIRN} \n')
    f.write('cp ${HOME}/CubismUP_3D/bin/simulation ./exec \n')

    f.write('export OMP_NUM_THREADS=12 \n')
    f.write('export CRAY_CUDA_MPS=1 \n')
    f.write('srun --ntasks 1 --ntasks-per-node=1 ./exec ${OPTIONS} \n')
    f.close()
    #os.system('sbatch daint_sbatch')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    description = "Compute a target file for RL agent from DNS data.")

    parser.add_argument('--printName', dest='printName', action='store_true',
    help="Only print run name.")
    parser.set_defaults(printName=False)

    parser.add_argument('--printOptions', dest='printOptions', action='store_true',
    help="Only print run options.")
    parser.set_defaults(printOptions=False)

    parser.add_argument('--launchDaint', dest='launchDaint', action='store_true',
    help="Only print run options.")
    parser.set_defaults(launchDaint=False)

    parser.add_argument('--launchEuler', dest='launchEuler', action='store_true',
    help="Only print run options.")
    parser.set_defaults(launchEuler=False)

    parser.add_argument('--case', type = int, default = -1,
    help="Simulation case.")

    args = parser.parse_args()

    NUS, EPS, RUN = [], [], []
    h = 2 * np.pi / 16 / 12
    for nu in np.logspace(np.log10(0.002), np.log10(0.02), 16) :
        for eps in np.logspace(np.log10(0.01), np.log10(2.0), 16) :
            if relFit(nu, eps) > 100 or relFit(nu, eps) < 20: continue
            if lambdaFit(nu, eps) > 0.1 * 2 * np.pi:          continue
            if etaFit(nu, eps) > h or etaFit(nu, eps) < h/8:  continue
            for i in [0, 1, 2] :
                NUS, EPS, RUN = NUS + [nu], EPS + [eps], RUN + [i]

    nCases = len(NUS)
    print('Defined %d cases' % nCases)

    if args.launchDaint: launchDaint(nCases)

    if args.case < 0: cases = range(nCases)
    else: cases = range(args.case, args.case + 1) # im a hack

    
    for i in cases:
        if args.printOptions:
            print( getSettings(NUS[i], EPS[i]) )
        if args.printName:
            print( runspec(NUS[i], EPS[i], RUN[i]) )
        if args.launchEuler:
            launchEuler(NUS[i], EPS[i], RUN[i])

#for nu in [0.002, 0.004, 0.008] :
# for eps in [0.02, 0.04, 0.08, 0.16, 0.32] :
#  tke0 = 2.77578963 * np.power(eps, (2.0/3.0) )
#  for scal in [2, 3] :
#  tke0 = 2.77578963 * np.power(eps, (2.0/3.0) )
#   for scal in [2] :
#  ext = scal * np.pi
#  os.system("\
#  export NU=%f \n\
#  export EPS=%f \n\
#  export TKE0=%f \n\
#  export EXT=%f \n\
#  echo $NU $EPS $TKE0 $EXT \n\
#  ./launchEuler.sh settingsHIT_DNS.sh HIT_DNS_EXT%dpi_EPS%.02f_NU%.03f"
#  % (nu, eps, tke0, ext, scal, eps, nu))
#for nu in [0.001, 0.002, 0.004, 0.008, 0.016] :
# for eps in [0.02, 0.04, 0.08, 0.16, 0.32, 0.64] :
