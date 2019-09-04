import os, numpy as np

for nu in [0.002, 0.004, 0.008] :
 for eps in [0.02, 0.04, 0.08, 0.16, 0.32] :
  tke0 = 2.77578963 * np.power(eps, (2.0/3.0) )
  for scal in [2, 3] :
    ext = scal * np.pi
    os.system("\
    export NU=%f \n\
    export EPS=%f \n\
    export TKE0=%f \n\
    export EXT=%f \n\
    echo $NU $EPS $TKE0 $EXT \n\
    ./launchEuler.sh settingsHIT_DNS.sh HIT_DNS_EXT%dpi_EPS%.02f_NU%.03f"
    % (nu, eps, tke0, ext, scal, eps, nu))
