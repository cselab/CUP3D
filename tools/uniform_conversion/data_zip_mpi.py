from mpi4py import MPI
import glob, string, h5py

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

myfiles = glob.glob("*-uniform.h5")
myfiles.sort()

tasks = len(myfiles)
my_share = int(tasks / size)
if tasks % size != 0 and rank == size - 1: #last rank gets what's left
    my_share += tasks % size;
my_share = int(my_share)
my_start = int(rank * int(tasks/size))
my_end = int(my_start + my_share)

for i in range(my_start,my_end):
     current_file = myfiles[i]
     Cfile = h5py.File(current_file, "r")
     dims = Cfile["/data"].shape
     print("Compressing: ", current_file, "dimensions: " , dims,flush=True)
     original_data = Cfile['data']
     f = h5py.File(current_file[:len(current_file) - 3]+"-compressed.h5","a" )
     dset = f.create_dataset("data", dims, compression="lzf",data=original_data)
     #dset = f.create_dataset("data", dims, compression="gzip",data=original_data,compression_opts=4)
     Cfile.close()
     f.close()

MPI.Finalize()
