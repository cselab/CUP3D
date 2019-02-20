ifeq ($(INTEL_LICENSE_FILE),)
compiler ?= gnu
else
compiler ?= intel
endif

ifeq "$(compiler)" "intel"
CXX ?= icpc
LD = $(CXX)
include make.icc.flags
else
CXX ?= g++
LD = $(CXX)
include make.gcc.flags
endif

NVCC?=nvcc
NVFLAGS?=-std=c++11 -code=sm_60 -arch=compute_60
LD=$(CC)

# this only affects amount of memory allocated by block processing:
nthreads ?= 24 # needs to be modified if user wants to use more threads per rank
#config ?= prod
config ?= debug
precision ?= double
precision_dump ?= single
bs ?= 16
hdf ?= true
rk2 ?= false
accfft ?= false
zlib ?= 0
verbose ?= false
bc ?= open
unboundedfft ?= false
dump_surface_binary ?= false

ifeq "$(dump_surface_binary)" "true"
        CPPFLAGS += -DCUP_DUMP_SURFACE_BINARY
endif

ifeq "$(bc)" "periodicz"
	CPPFLAGS += -DBC_PERIODICZ
endif

ifeq "$(unboundedfft)" "true"
	CPPFLAGS += -DCUP_UNBOUNDED_FFT
endif

ifeq "$(verbose)" "true"
	CPPFLAGS += -DCUP_VERBOSE
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

ifeq "$(rk2)" "true"
	CPPFLAGS += -DCUP_RK2
endif

ifeq "$(precision)" "single"
	CPPFLAGS += -DCUP_SINGLE_PRECISION
	NVFLAGS += -DCUP_SINGLE_PRECISION
endif

ifeq "$(precision_dump)" "double"
	CPPFLAGS += -DCUP_HDF5_DOUBLE_PRECISION
endif

CPPFLAGS += -DCUP_ALIGNMENT=64
CPPFLAGS += -DCUP_BLOCK_SIZE=$(bs)
#CPPFLAGS += -lm

#CPPFLAGS += -I${HDF5_INCLUDE_OPTS}

# needed in order to compile on nodes
# include make.daint
ifeq ($(shell uname -s), Darwin)
CPPFLAGS += -dynamic  # Something specific to Darwin linker(?)
CC=mpic++
LD=mpic++
accfft=false
endif

ifneq "$(findstring nid,$(shell hostname))" ""
	include make.daint
endif
ifneq "$(findstring daint,$(shell hostname))" ""
	include make.daint
endif
ifneq "$(findstring panda,$(shell hostname))" ""
       include make.panda
endif
ifneq "$(findstring euler,$(shell hostname))" ""
include make.euler
endif
ifneq "$(findstring eu-,$(shell hostname))" ""
	include make.euler
endif

CPPFLAGS += -I../source/ -I../Cubism/include/
#CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

LDFLAGS = $(CPPFLAGS)

ifeq "$(hdf)" "true"
CPPFLAGS+= -DCUBISM_USE_HDF #-I$(hdf-inc)
ifneq ($(HDF5_LIB),)
	LIBS+= -L$(HDF5_LIB) -lhdf5 #-L$(hdf-lib) -lhdf5
else
	LIBS+= -lhdf5 #-L$(hdf-lib) -lhdf5
endif
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
	FFTW_LIBS = -lfftw3f_mpi -lfftw3f_omp -lfftw3f
else
	FFTW_LIBS = -lfftw3_mpi -lfftw3_omp -lfftw3
endif
else
ifeq "$(precision)" "single"
	FFTW_LIBS = $(FFTW_DIR)/libfftw3f_threads.a $(FFTW_DIR)/libfftw3f_mpi.a $(FFTW_DIR)/libfftw3f.a
else
	FFTW_LIBS = $(FFTW_DIR)/libfftw3_threads.a $(FFTW_DIR)/libfftw3_mpi.a $(FFTW_DIR)/libfftw3.a
endif
endif

BUILDDIR = .
#SRC_DIR = $(BUILDDIR)/../source/
SRC_DIR = $(sort $(dir $(wildcard $(BUILDDIR)/../source/*/)))
CUBISM_DIR = $(BUILDDIR)/../Cubism/src/


OBJECTS = IF3D_ObstacleFactory.o IF3D_ObstacleOperator.o IF3D_ObstacleVector.o \
	IF3D_FishOperator.o IF3D_StefanFishOperator.o IF3D_CarlingFishOperator.o \
	IF3D_SphereObstacleOperator.o IF3D_PlateObstacleOperator.o BufferedLogger.o \
	IF3D_NacaOperator.o IF3D_FishLibrary.o IF3D_ExternalObstacleOperator.o \
	IF3D_DCylinderObstacleOperator.o IF3D_TestDiffusionObstacleOperator.o \
	SimulationData.o Simulation.o AdvectionDiffusion.o ComputeDissipation.o \
	ExternalForcing.o FadeOut.o InitialConditions.o ObstacleManagement.o \
	PressurePenalization.o PressureRHS.o \
	ArgumentParser.o  # Temporary solution for Cubism .cpp files.
	#IF3D_ElasticFishOperator.o \

ifeq "$(accfft)" "true"
	CPPFLAGS += -D_ACCFFT_ -I$(accfft-inc) -D_CUDA_COMP_ -Wno-float-equal
	NVFLAGS += -I$(accfft-inc) -I$(MPICH_DIR)/include/
	LIBS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft
	NVOBJECTS = PoissonSolverScalarACC.o
else
	NVOBJECTS =
	LIBS += $(FFTW_LIBS)
endif

#LIBS += -L../tools/fpzip/lib -lfpzip
LIBS += -lgsl -lgslcblas -lstdc++

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
	$(LD) $(EXEOBJ) $(NVOBJECTS) $(LDFLAGS) $(LIBS) -o $@

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
