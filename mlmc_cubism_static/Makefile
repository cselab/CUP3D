CFLAGS = -Wall -I../src -I. -DMASTER_ONLY_GENERATES=1 -fopenmp
LDFLAGS = -L../makefiles -lpthread -lm -lz -lhdf5 -lfftw3 -lfftw3_threads -lfftw3_mpi -lgsl -lgslcblas -fopenmp

ifneq "$(findstring daint,$(shell hostname))" ""
	CC = CC
	accfft-inc = /users/novatig/accfft/build_dbg/include/
	accfft-lib = /users/novatig/accfft/build_dbg/lib/ 	
endif

LDFLAGS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft

all: ms 

ms: ms.cpp main_util.cpp cubism # ../src/libtorc.a
	$(CC) $(CFLAGS) -o $@ ms.cpp main_util.cpp -L. -lsimulation $(LDFLAGS)

cubism:
	make -C ../makefiles config=prod lib -j
clean:
	rm -f ms
	make -C ../makefiles clean
