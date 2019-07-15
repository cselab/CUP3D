ifeq ($(shell uname -s), Darwin)
	CPPFLAGS+=-dynamic  # Something specific to Darwin linker(?)
	CXX=mpic++
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
ifneq "$(findstring falcon,$(shell hostname))" ""
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
LD=$(CXX)

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
smarties ?= false
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
	FishShapes.o IterativePressurePenalization.o IterativePressureNonUniform.o \
	FixedMassFlux_nonUniform.o SGS.o Analysis.o SpectralManip.o SpectralIcGenerator.o \
	SpectralAnalysis.o SpectralForcing.o ArgumentParser.o \
	#ElasticFishOperator.o # Temporary solution for Cubism .cpp files.

#################################################
# HDF5
#################################################
ifeq "$(hdf)" "true"
	# If you have issues with linking HDF5 so files, add $HDF5_ROOT/lib to your LD_LIBRARY_PATH.
	# For that reason we do not automatically check ../dependencies/build/hdf5-*.
	ifneq ($(HDF5_ROOT),)
		# OK
	else ifneq ($(HDF5ROOT),)
		HDF5_ROOT = HDF5ROOT
	else
		ERROR = "Neither HDF5_ROOT nor HDF5ROOT environment variables set. Did you forget to `module load`?"
	endif
	LIBS     += -L$(HDF5_ROOT)/lib -lhdf5
	CPPFLAGS += -I$(HDF5_ROOT)/include -DCUBISM_USE_HDF
endif

#################################################
# FFTW
#################################################
ifneq ($(FFTW_ROOT),)
	# OK
else ifneq ($(FFTW3_ROOT_DIR),)
	FFTW_ROOT = $(FFTW3_ROOT_DIR)
else ifneq ($(FFTW_DIR),)
	FFTW_ROOT = $(FFTW_DIR)/..
else
	ERROR = "None of the FFTW_ROOT, FFTW_ROOT_DIR or FFTW_DIR environment variables set. Did you forget to `module load`?"
endif
FFTW_LIBS  = -L$(FFTW_ROOT)/lib
CPPFLAGS  += -I$(FFTW_ROOT)/include

ifneq ($(ERROR),)
$(error $(ERROR))  # Easiest way I found to avoid unindenting the code above...
endif

ifeq "$(precision)" "single"
	FFTW_LIBS += -lfftw3f_mpi -lfftw3f_omp -lfftw3f
else
	FFTW_LIBS += -lfftw3_mpi -lfftw3_omp -lfftw3
endif

#################################################
# ACCFFT
#################################################
ifeq "$(accfft)" "true"
	ifneq ($(ACCFFT_ROOT),)
		# OK
	else ifneq (,$(wildcard ../dependencies/build/accfft))
		ACCFFT_ROOT = ../dependencies/build/accfft
	else
		ERROR = "AccFFT not found, aborting."
	endif
	LIBS     += -L$(ACCFFT_ROOT)/lib
	LIBS += -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft
	CPPFLAGS += -I$(ACCFFT_ROOT)/include
	CPPFLAGS += -D_ACCFFT_ -Wno-float-equal
	NVFLAGS  += -I$(ACCFFT_ROOT)/include
	NVFLAGS  += -I$(MPICH_DIR)/include/
	NVOBJECTS = PoissonSolverACC_kernels.o
	OBJECTS += PoissonSolverACCUnbounded.o PoissonSolverACCPeriodic.o
else
	OBJECTS += PoissonSolverUnbounded.o PoissonSolverPeriodic.o
	NVOBJECTS =
	LIBS += $(FFTW_LIBS)
endif

#################################################
# HYPRE & PETSC & MISC
#################################################
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

#################################################

EXEOBJ = $(OBJECTS) main.o
LIBOBJ = $(OBJECTS) cubism_main.o CubismUP_3D.o
VPATH := $(VPATH):$(SRC_DIR):$(CUBISM_DIR)
ALLOBJ := $(sort $(EXEOBJ) $(LIBOBJ) $(NVOBJECTS))  # `sort` removes duplicates.
DEPS := $(ALLOBJ:%.o=%.d)

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
	$(CXX) $(CPPFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

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
