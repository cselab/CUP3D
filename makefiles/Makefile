#CC ?= g++
CC ?= mpic++
NVCC?=nvcc
NVFLAGS?=-O3 -std=c++11 -code=sm_60 -arch=compute_60
LD=$(CC)

# this only affects amount of memory allocated by block processing:
nthreads ?= 24 # needs to be modified if user wants to use more threads per rank
#config ?= prod
config ?= debug
precision ?= double
precision_dump ?= single
bs ?= 16
#bs ?= 32
bsx ?= $(bs)
bsy ?= $(bs)
bsz ?= $(bs)
hdf ?= true
vtk ?= false
rk2 ?= false
accfft ?= false
zlib ?= 0
verbose ?= false
bc ?= open
unboundedfft ?= false
dumpraw ?= false

ifeq "$(dumpraw)" "true"
        CPPFLAGS += -D_DUMP_RAW_
endif

ifeq "$(bc)" "periodicz"
	CPPFLAGS += -DBC_PERIODICZ
endif

ifeq "$(unboundedfft)" "true"
	CPPFLAGS += -D_UNBOUNDED_FFT_
endif

ifeq "$(verbose)" "true"
	CPPFLAGS += -D_VERBOSE_
endif

rltrain ?= false
ifeq "$(rltrain)" "true"
	CPPFLAGS += -DRL_LAYER
	CPPFLAGS += -I${HOME}/smarties/source/ -D__RL_MPI_CLIENT
endif

ifeq "$(zlib)" "1"
	CPPFLAGS += -D_USE_ZLIB_ -D_USE_WAVZ_
	LIBS +=  -lz -lrt -lstdc++
endif

CPPFLAGS+= -D_ALIGNBYTES_=32 -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align
CPPFLAGS+= -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls
CPPFLAGS+= -Wshadow -Woverloaded-virtual -Wuninitialized -pedantic
CPPFLAGS+= -Wno-unused-parameter # does not make much sense with function overloading

ifeq "$(precision)" "single"
	CPPFLAGS += -D_FLOAT_PRECISION_
	NVFLAGS += -D_FLOAT_PRECISION_
endif

ifeq "$(precision_dump)" "double"
	CPPFLAGS += -D_HDF5_DOUBLE_PRECISION_
endif

CPPFLAGS += -D_BLOCKSIZE_=$(bs)
CPPFLAGS += -D_BS_=$(bs)
CPPFLAGS += -D_BLOCKSIZEX_=$(bsx)
CPPFLAGS += -D_BLOCKSIZEY_=$(bsy)
CPPFLAGS += -D_BLOCKSIZEZ_=$(bsz)
CPPFLAGS += -I${HDF5_INCLUDE_OPTS}

ifeq "$(config)" "prod"
	#CPPFLAGS += -DNDEBUG -O2 -fstrict-aliasing -march=core-avx2 -mtune=core-avx2 -ffast-math -falign-functions=32
	CPPFLAGS += -DNDEBUG -O2 -fstrict-aliasing -march=native -mtune=native -ffast-math -falign-functions=32
	#CPPFLAGS += -DNDEBUG -O2 -fstrict-aliasing -ffast-math #-falign-functions=32
else
	CPPFLAGS += -g -O0
	#CPPFLAGS += -Wl,-u,PMI_Init
	#RED='\033[0;31m'
	#NC='\033[0m' # No Color
$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif

ifeq "$(rk2)" "true"
	CPPFLAGS += -D_RK2_
endif

# needed in order to compile on nodes
# include make.daint
ifeq ($(shell uname -s), Darwin)
CPPFLAGS += -dynamic  # Something specific to Darwin linker(?)
CC=mpic++
LD=mpic++
accfft=false
endif

ifneq "$(findstring daint,$(shell hostname))" ""
	#include make.daint_mpich
	include make.daint
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
	include make.brutus
endif

ifneq "$(findstring euler,$(shell hostname))" ""
include make.euler
endif
ifneq "$(findstring eu-,$(shell hostname))" ""
	include make.euler
endif

CPPFLAGS += -std=c++11 -fopenmp
#NVFLAGS += -std=c++11
CPPFLAGS += -DNTHREADS=$(nthreads)
CPPFLAGS += -I../source/ -I../Cubism/source
#CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

LDFLAGS = $(CPPFLAGS)

ifeq "$(vtk)" "true"
CPPFLAGS += -D_VTK_ -I$(vtk-inc)
LIBS +=	-L$(vtk-lib) \
	-lvtkHybrid \
	-lvtkVolumeRendering \
	-lvtkRendering \
	-lvtkIO \
	-lvtkGenericFiltering \
	-lvtkGraphics \
	-lvtkImaging \
	-lvtkFiltering \
	-lvtkCommon \
	-lvtkftgl \
	-lvtkfreetype \
	-lvtkDICOMParser \
	-lvtkexpat \
	-lvtktiff \
	-lvtkpng \
	-lvtkjpeg \
	-lvtkzlib \
	-lvtksys
endif

ifeq "$(hdf)" "true"
CPPFLAGS+= -D_USE_HDF_ #-I$(hdf-inc)
LIBS+= -lhdf5 #-L$(hdf-lib) -lhdf5
endif

#################################################
# FFTW
#################################################

# `FFTW_INC` and `FFTW_DIR` are for custom installation. Names are in fact
# compatible with `module load`, so the `else` statement should work for daint
# too.
ifneq ($(FFTW_INC),)
	CPPFLAGS += -I$(FFTW_INC)
endif
ifeq ($(FFTW_DIR),)
ifeq "$(precision)" "single"
	FFTW_LIBS = -lfftw3f_threads -lfftw3f_mpi -lfftw3f
else
	FFTW_LIBS = -lfftw3_omp -lfftw3_mpi -lfftw3
endif
else
ifeq "$(precision)" "single"
	FFTW_LIBS = $(FFTW_DIR)/libfftw3f_threads.a $(FFTW_DIR)/libfftw3f_mpi.a $(FFTW_DIR)/libfftw3f.a
else
	FFTW_LIBS = $(FFTW_DIR)/libfftw3_threads.a $(FFTW_DIR)/libfftw3_mpi.a $(FFTW_DIR)/libfftw3.a
endif
endif


#LIBS += -L../tools/fpzip/lib -lfpzip

LIBS += -lstdc++ -lgsl -lgslcblas -lm -ldl -lz -lpthread

BUILDDIR = .
SRC_DIR = $(BUILDDIR)/../source/
CUBISM_DIR = $(BUILDDIR)/../Cubism/

OBJECTS =  IF3D_ObstacleFactory.o \
           IF3D_ObstacleOperator.o \
           IF3D_FishLibrary.o \
           IF3D_ObstacleVector.o \
           IF3D_StefanFishOperator.o \
           IF3D_FishOperator.o \
           IF3D_SphereObstacleOperator.o \
           IF3D_PlateObstacleOperator.o \
           IF3D_ExternalObstacleOperator.o \
           IF3D_ElasticFishOperator.o \
           Simulation.o \
           IF3D_NacaOperator.o \
           BufferedLogger.o \
		   IF3D_CarlingFishOperator.o \
           IF3D_DCylinderObstacleOperator.o

ifeq "$(accfft)" "true"
	CPPFLAGS += -D_ACCFFT_ -I$(accfft-inc) -D_CUDA_COMP_ -Wno-float-equal
	LIBS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft
	NVOBJECTS = PoissonSolverScalarACC.o
else
	NVOBJECTS =
	LIBS += $(FFTW_LIBS)
endif

EXEOBJ = $(OBJECTS)  main.o
LIBOBJ = $(OBJECTS) cubism_main.o CubismUP_3D.o
LIBRLOBJ = $(OBJECTS) ext_app_main.o
VPATH := $(VPATH):$(SRC_DIR):$(CUBISM_DIR)
ALLOBJ := $(sort $(EXEOBJ) $(LIBOBJ) $(LIBRLOBJ) $(NVOBJECTS))  # `sort` removes duplicates.
DEPS := $(ALLOBJ:%.o=%.d)

#ifeq "$(rltrain)" "true"
#EXEOBJ += Communicator.o
#VPATH = $(VPATH):${HOME}/smarties/source/
#endif

.DEFAULT_GOAL := ../bin/simulation

../bin/simulation: $(EXEOBJ) $(NVOBJECTS)
	mkdir -p ../bin
	$(LD) $(EXEOBJ) $(NVOBJECTS) $(LIBS) $(LDFLAGS) -o $@

lib: $(LIBOBJ) $(NVOBJECTS)
	mkdir -p ../lib
	ar rs ../lib/libcubismup3d.a $(FFTW_LIBS) $(LIBOBJ) $(NVOBJECTS)

librl: $(LIBRLOBJ) $(NVOBJECTS)
	ar rs libsim_rl.a $(LIBRLOBJ) $(NVOBJECTS)

-include $(DEPS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

PoissonSolverScalarACC.o: PoissonSolverScalarACC.cu
	nvcc $(NVFLAGS) -c $< -o $@

clean:
	rm -f $(EXEOBJ) $(LIBOBJ)
	rm -f $(DEPS)
	rm -f ../bin/simulation Communicator.o ../lib/libcubismup3d.a PoissonSolverScalar*.o
	rmdir ../bin 2> /dev/null || true
	rmdir ../lib 2> /dev/null || true

ifeq "$(config)" "prod"
else
        #RED='\033[0;31m'
        #NC='\033[0m' # No Color
$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif
