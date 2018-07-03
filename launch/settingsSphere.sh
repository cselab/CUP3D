BASENAME=FlowPastSphere
NNODE=4

FACTORY='IF3D_Sphere L=0.075 xpos=0.25 xvel=0.2 computeForces=0 bForcedInSimFrame=1 bFixFrameOfRef=1
'
# for accel and decel start and stop add accel=1 T=time_for_accel
# shift center to shed vortices immediately by ypos=0.250244140625 zpos=0.250244140625

OPTIONS=
#OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -dump2D 1 -dump3D 0"
#OPTIONS+=" -Wtime ${WSECS}"
OPTIONS+=" -nprocsx ${NNODE}"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -lambda 1e5"
OPTIONS+=" -nu 0.00001"
OPTIONS+=" -tend 10"
