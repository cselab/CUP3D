make clean; make bgq=0 qpxemu=1 hdf=1
cd ../ch2hdf
make clean; make
cd -

sbatch submit_daint