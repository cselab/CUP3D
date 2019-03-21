ifeq ($(shell uname -s), Darwin)
	CPPFLAGS+=-dynamic  # Something specific to Darwin linker(?)
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
hypre ?= false
petsc ?= false
accfft ?= false
zlib ?= 0
verbose ?= false
dump_surface_binary ?= false

ifeq "$(dump_surface_binary)" "true"
	CPPFLAGS += -DCUP_DUMP_SURFACE_BINARY
endif

ifeq "$(verbose)" "true"
	CPPFLAGS += -DCUP_VERBOSE
endif

ifeq "$(zlib)" "1"
	CPPFLAGS += -D_USE_ZLIB_ -D_USE_WAVZ_
	LIBS +=  -lz -lrt -lstdc++
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
CPPFLAGS += -DCUP_NO_MACROS_HEADER
#CPPFLAGS += -lm

#CPPFLAGS += -I${HDF5_INCLUDE_OPTS}


CPPFLAGS += -I../source/ -I../Cubism/include/
NVFLAGS += -I../source/
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
ifneq ($(FFTW_INC),)   # If not empty.
	CPPFLAGS += -I$(FFTW_INC)
endif

ifeq ($(FFTW_DIR),)
ifeq "$(precision)" "single"
	FFTW_LIBS = -lfftw3f_mpi -lfftw3f_omp -lfftw3f
else
	FFTW_LIBS = -lfftw3_mpi -lfftw3_omp -lfftw3
endif
else
	FFTW_LIBS = -L$(FFTW_DIR) -L$(FFTW_DIR)/lib
ifeq "$(precision)" "single"
	FFTW_LIBS += -fftw3f_threads -lfftw3f_mpi -lfftw3f
else
	FFTW_LIBS += -lfftw3_threads -lfftw3_mpi -lfftw3
endif
endif

BUILDDIR = .
#SRC_DIR = $(BUILDDIR)/../source/
SRC_DIR = $(sort $(dir $(wildcard $(BUILDDIR)/../source/*/)))
CUBISM_DIR = $(BUILDDIR)/../Cubism/src/

OBJECTS = ObstacleFactory.o Obstacle.o ObstacleVector.o Ellipsoid.o Cylinder.o \
	Fish.o StefanFish.o CarlingFish.o Sphere.o Plate.o ExternalObstacle.o Naca.o \
	FishLibrary.o BufferedLogger.o SimulationData.o Simulation.o PoissonSolver.o \
	PoissonSolverMixed.o AdvectionDiffusion.o ComputeDissipation.o PressureRHS.o \
	PressureProjection.o Penalization.o InitialConditions.o FluidSolidForces.o \
	ObstaclesCreate.o ObstaclesUpdate.o ExternalForcing.o FadeOut.o \
	FishShapes.o ArgumentParser.o # Temporary solution for Cubism .cpp files.
	#ElasticFishOperator.o

ifeq "$(hypre)" "true"
ifeq "$(compiler)" "gnu"
	# I don't want to see g++ complaining about hypre:
	CPPFLAGS+= -Wno-redundant-decls
endif
LIBS+= -lHYPRE
CPPFLAGS+= -DCUP_HYPRE
OBJECTS+= PoissonSolverHYPREMixed.o
endif

ifeq "$(petsc)" "true"
	CPPFLAGS+= -DCUP_PETSC
	OBJECTS+= PoissonSolverPETSCMixed.o
endif

ifeq "$(accfft)" "true"
	CPPFLAGS += -D_ACCFFT_ -I$(accfft-inc) -Wno-float-equal
	NVFLAGS += -I$(accfft-inc) -I$(MPICH_DIR)/include/
	LIBS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft
	NVOBJECTS = PoissonSolverACC_kernels.o
	OBJECTS += PoissonSolverACCUnbounded.o PoissonSolverACCPeriodic.o
else
	OBJECTS += PoissonSolverUnbounded.o PoissonSolverPeriodic.o
	NVOBJECTS =
	LIBS += $(FFTW_LIBS)
endif

#LIBS += -L../tools/fpzip/lib -lfpzip
ifneq "$(GSL_ROOT_DIR)" ""
	CPPFLAGS += -I$(GSL_ROOT_DIR)/include
	LIBS += -L$(GSL_ROOT_DIR)/lib
endif
LIBS += -lgsl -lgslcblas

ifneq "$(smarties)" "false"
	SMARTIESDIR ?= ../../smarties/source/
	CPPFLAGS += -I$(SMARTIESDIR)/Communicators/ -DSMARTIES_APP
endif

EXEOBJ = $(OBJECTS) main.o
LIBOBJ = $(OBJECTS) cubism_main.o CubismUP_3D.o
VPATH := $(VPATH):$(SRC_DIR):$(CUBISM_DIR)
ALLOBJ := $(sort $(EXEOBJ) $(LIBOBJ) $(NVOBJECTS))  # `sort` removes duplicates.
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

testRL: main_RL_test.o SGS_RL.o $(EXEOBJ) $(NVOBJECTS)
	ar rs libsimulation.a $^

-include $(DEPS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

PoissonSolverACC_kernels.o: PoissonSolverACC_kernels.cu
	$(NVCC) $(NVFLAGS) -c $< -o $@

clean:
	rm -f $(EXEOBJ) $(LIBOBJ)
	rm -f $(DEPS)
	rm -f PoissonSolver*.o PoissonSolver*.d
	rm -f ../bin/simulation Communicator.o ../lib/libcubismup3d.a PoissonSolverScalar*.o
	rmdir ../bin 2> /dev/null || true
	rmdir ../lib 2> /dev/null || true

ifeq "$(config)" "prod"
else
        #RED='\033[0;31m'
        #NC='\033[0m' # No Color
$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif
