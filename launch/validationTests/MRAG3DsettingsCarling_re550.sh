NNODE=1

# build the settings
SETTINGS=
SETTINGS+=" -dumpfreq 100"
SETTINGS+=" -cfl 0.9"
SETTINGS+=" -jump 2"
SETTINGS+=" -fmm-theta 0.5"
SETTINGS+=" -uniform 0"
SETTINGS+=" -particles 1"
SETTINGS+=" -fmm potential-aggressive-new"
SETTINGS+=" -vtu 1"
SETTINGS+=" -test carlingfish"
SETTINGS+=" -nsteps 100000"
SETTINGS+=" -adaptfreq 5"
SETTINGS+=" -correctionfreq 0"
SETTINGS+=" -lcfl 0.01"
SETTINGS+=" -savefreq 100"
SETTINGS+=" -bpd 8"
SETTINGS+=" -tmax 1000000000.0"
SETTINGS+=" -lmax 7"
#SETTINGS+=" -lmax 6"
SETTINGS+=" -f2z 1"
SETTINGS+=" -rtol 1e-3"
SETTINGS+=" -ctol 1e-5"
SETTINGS+=" -re 550"
SETTINGS+=" -obstacle carling"
SETTINGS+=" -ramp 100"
SETTINGS+=" -D 0.1"
SETTINGS+=" -lambda 10000"
SETTINGS+=" -correction 0"
SETTINGS+=" -xpos 0.65"
SETTINGS+=" -L 0.25"
SETTINGS+=" -T 1.0"
SETTINGS+=" -fadeout 1"
#SETTINGS+=" -mollfactor -1"
SETTINGS+=" -mollfactor 4"
SETTINGS+=" -obstaclepotential potential-aggressive-new"
