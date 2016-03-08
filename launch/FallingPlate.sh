cd ../makefiles/
make clean
make config=production nthreads=12 bc=mixed precision=single particles=false dlm=true hdf=true movingframe=true bs=16 bsx=8 -j
mv simulation simulation_plate
cd ../launch/

sbatch submit_fallingPlate
