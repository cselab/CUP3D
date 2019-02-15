BASENAME=ChannelScaling_04_BX${CUPBX}_BY${CUPBY}_BZ${CUPBZ}

NNODEX=4
NNODEY=1
NNODE=$(($NNODEX * $NNODEY))

FACTORY=''

OPTIONS=
OPTIONS+=" -bpdx 16 -bpdy 16 -bpdz 4"
OPTIONS+=" -dump2D 1 -dump3D 0"
OPTIONS+=" -BC_x periodic -BC_y wall -BC_z periodic"
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1 -tend 10 -uMax_forced 1"
# RE = 8/12 Ly Umax / \nu = 66.666
OPTIONS+=" -nu 0.01"
