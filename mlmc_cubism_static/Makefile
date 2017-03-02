CFLAGS = -Wall -I../src -I. -DMASTER_ONLY_GENERATES=1 -fopenmp
LDFLAGS = -L../src -lpthread -lm -lz -lhdf5 -lfftw3 -lfftw3_threads -lfftw3_mpi -lgsl -lgslcblas -fopenmp

ifneq "$(findstring daint,$(shell hostname))" ""
	CC = CC
	accfft-inc = /users/novatig/accfft/build_again/include/
	accfft-lib = /users/novatig/accfft/build_again/lib/ 	
endif

LDFLAGS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft

all: ms 

ms: ms.cpp main_util.cpp # ../src/libtorc.a
	$(CC) $(CFLAGS) -o $@ $^ -L. -lsimulation $(LDFLAGS)

clean:
	rm -f ms
