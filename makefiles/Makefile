#CC ?= g++
CC ?= mpic++
NVCC?=nvcc
NVFLAGS?=-O3 -code=sm_60 -arch=compute_60
LD=$(CC)

nthreads ?= 24
#config ?= prod
config ?= debug
precision ?= double
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

ifeq "$(bc)" "periodicz"
	CPPFLAGS += -DBC_PERIODICZ
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

CPPFLAGS+= -w -D_ALIGNBYTES_=32 -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align
CPPFLAGS+= -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls
CPPFLAGS+= -Wshadow -Woverloaded-virtual -Wuninitialized -Wno-div-by-zero -pedantic

ifeq "$(precision)" "single"
	CPPFLAGS += -D_SP_COMP_
	NVCCFLAGS += -D_SP_COMP_
endif

CPPFLAGS += -D_BLOCKSIZE_=$(bs)
CPPFLAGS += -D_BS_=$(bs)
CPPFLAGS += -D_BSX_=$(bsx)
CPPFLAGS += -D_BSY_=$(bsy)
CPPFLAGS += -D_BSZ_=$(bsz)
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
CPPFLAGS += -I../source/ -I../Cubism/
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

ifeq "$(precision)" "single"
	LIBS += -lfftw3f -lfftw3f_threads -lfftw3f_mpi
endif
ifeq "$(precision)" "double"
	LIBS += -lfftw3 -lfftw3_threads -lfftw3_mpi
endif

#LIBS += -L../tools/fpzip/lib -lfpzip

LIBS += -lstdc++ -lgsl -lgslcblas -lm -ldl -lz -lpthread

BUILDDIR = .
SRC_DIR = $(BUILDDIR)/../source/
CUBISM_DIR = $(BUILDDIR)/../Cubism/

OBJECTS =  IF3D_VortexOperator.o \
           IF3D_ObstacleFactory.o \
           IF3D_ObstacleOperator.o \
           IF3D_FishLibrary.o \
           IF3D_ObstacleVector.o \
           IF3D_StefanFishOperator.o \
           IF3D_FishOperator.o \
           IF3D_ForcedSphereObstacleOperator.o \
           IF3D_SphereObstacleOperator.o \
           IF3D_PlateObstacleOperator.o \
           Profiler.o \
           Simulation.o \
           ProcessOperatorsOMP.o \
           IF3D_NacaOperator.o \
           IF3D_DCylinderObstacleOperator.o
           #IF3D_DeadFishOperator.o \
           #IF3D_CarlingFishOperator.o \

ifeq "$(accfft)" "true"
	CPPFLAGS += -D_ACCFFT_ -I$(accfft-inc) -D_CUDA_COMP_
	LIBS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft
	NVOBJECTS = PoissonSolverScalarFFTW_ACC.o
else
	NVOBJECTS =
	OBJECTS += PoissonSolverScalarFFTW.o PoissonSolverScalarFFTW_MPI.o
endif

EXEOBJ = $(OBJECTS)  main.o
LIBOBJ = $(OBJECTS) cubism_main.o
LIBRLOBJ = $(OBJECTS) ext_app_main.o
VPATH := $(VPATH):$(SRC_DIR):$(CUBISM_DIR)
ALLOBJ := $(sort $(EXEOBJ) $(LIBOBJ) $(LIBRLOBJ) $(NVOBJECTS))  # `sort` removes duplicates.
DEPS := $(ALLOBJ:%.o=%.d)

#ifeq "$(rltrain)" "true"
#EXEOBJ += Communicator.o
#VPATH = $(VPATH):${HOME}/smarties/source/
#endif

.DEFAULT_GOAL := simulation

simulation: $(EXEOBJ) $(NVOBJECTS)
	$(LD) $(EXEOBJ) $(NVOBJECTS) $(LIBS) $(LDFLAGS) -o $@

lib: $(LIBOBJ) $(NVOBJECTS)
	ar rs libsimulation.a $(LIBOBJ) $(NVOBJECTS)

librl: $(LIBRLOBJ) $(NVOBJECTS)
	ar rs libsim_rl.a $(LIBRLOBJ) $(NVOBJECTS)

-include $(DEPS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) $(LIBS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

PoissonSolverScalarFFTW_ACC.o: PoissonSolverScalarFFTW_ACC.cu
	nvcc $(NVFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPS)
	rm -f simulation Communicator.o libsimulation.a main.o cubism_main.o PoissonSolverScalarFFTW*.o

ifeq "$(config)" "prod"
else
        #RED='\033[0;31m'
        #NC='\033[0m' # No Color
$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif
