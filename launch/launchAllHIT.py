import os, numpy as np
nNus = 16
nEps = 16
COEFREL = 7.10538
h = 2 * np.pi / 16 / 12

def    relFit(nu, eps): return 7.10538 * np.power(eps, 1/6.0) / np.sqrt(nu)

def    etaFit(nu, eps): return np.power(eps, 0.25) * np.power(nu, 0.75)

def lambdaFit(nu, eps): return 5.2623 * np.power(eps,-1/6.0) * np.sqrt(nu);

cases = 0
for nu in np.logspace(np.log10(0.002), np.log10(0.02), nNus) :
  for eps in np.logspace(np.log10(0.01), np.log10(2.0), nEps) :
    if relFit(nu,eps) > 100 or relFit(nu,eps) < 20: continue
    if lambdaFit(nu,eps) > 0.1 * 2 * np.pi: continue
    if etaFit(nu,eps) > h or etaFit(nu,eps) < h/8: continue
    print("HIT_RE_EPS%.02f_NU%.04f" % (eps,nu) )
    cases += 1
    '''
    ext = scal * np.pi
    os.system("\
    export NU=%f \n\
    export EPS=%f \n\
    echo $NU $EPS \n\
    ./launchEuler.sh settingsHIT_DNS.sh HIT_DNS2_EXT2_EPS%.03f_NU%.04f"
    % (nu, eps, eps, nu))
    '''
print(cases)

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
