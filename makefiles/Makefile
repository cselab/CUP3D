ifeq ($(shell uname -s), Darwin)
include make.macos
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

# this only affects amount of memory allocated by block processing:
nthreads ?= 24 # needs to be modified if user wants to use more threads per rank
config ?= prod
precision ?= double
precision_dump ?= single
usemap ?= false
bs ?= 8
bsx ?= 8
bsy ?= 8
bsz ?= 8
hdf ?= true
zlib ?= 0
verbose ?= false
dump_surface_binary ?= false

ifeq "$(shell $(CXX) --version | grep -ci 'icpc')" "0"
ifeq "$(shell $(CXX) --version | grep -ci 'llvm')" "0"
compiler=gnu
include make.gcc.flags
else
compiler=llvm
include make.llvm.flags
endif
else
compiler=intel
include make.icc.flags
endif

NVCC?=nvcc
NVFLAGS?=-std=c++14 -code=sm_60 -arch=compute_60
LD=$(CXX)

ifeq "$(usemap)" "true"
	CPPFLAGS += -DCUBISM_USE_MAP
endif

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
#CPPFLAGS += -DCUP_BLOCK_SIZE=$(bs)
CPPFLAGS += -DCUP_BLOCK_SIZEX=$(bsx)
CPPFLAGS += -DCUP_BLOCK_SIZEY=$(bsy)
CPPFLAGS += -DCUP_BLOCK_SIZEZ=$(bsz)
CPPFLAGS += -DCUP_NO_MACROS_HEADER
NVFLAGS  += -DCUP_NO_MACROS_HEADER
#CPPFLAGS += -lm
#CPPFLAGS += -I${HDF5_INCLUDE_OPTS}
CPPFLAGS += -I../Cubism/include/
#CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

LDFLAGS = $(CPPFLAGS)

BUILDDIR = .
#SRC_DIR = $(BUILDDIR)/../source/
SRC_DIR = $(sort $(dir $(wildcard $(BUILDDIR)/../source/*/)))
CUBISM_DIR = $(BUILDDIR)/../Cubism/src/

OBJECTS = ObstacleFactory.o Obstacle.o ObstacleVector.o Ellipsoid.o Cylinder.o \
	Fish.o StefanFish.o CarlingFish.o Sphere.o Plate.o ExternalObstacle.o Naca.o \
	FishLibrary.o BufferedLogger.o SimulationData.o Simulation.o \
	AdvectionDiffusion.o ComputeDissipation.o PressureRHS.o \
	PressureProjection.o Penalization.o InitialConditions.o FluidSolidForces.o \
	ObstaclesCreate.o ObstaclesUpdate.o ExternalForcing.o \
	FishShapes.o \
	SGS.o FixedMassFlux_nonUniform.o Analysis.o \
	ArgumentParser.o MyClock.o \
	Checkpoint.o PoissonSolverAMR.o
	#ElasticFishOperator.o # Temporary solution for Cubism .cpp files.

#################################################
# HDF5
#################################################
ifeq "$(hdf)" "true"
	# If you have issues with linking HDF5 so files, add $HDF5_ROOT/lib to your LD_LIBRARY_PATH.
	# For that reason we do not automatically check ../dependencies/build/hdf5-*.
	ifneq ($(HDF5_ROOT),)
		# OK, do not overwrite HDF5_ROOT
	else ifneq ($(HDF5ROOT),)
		HDF5_ROOT = $(HDF5ROOT)
	endif
	ifneq ($(HDF5_ROOT),)
		LIBS     += -L$(HDF5_ROOT)/lib
		CPPFLAGS += -I$(HDF5_ROOT)/include
	endif
	LIBS     += -lhdf5
	CPPFLAGS += -DCUBISM_USE_HDF
endif

#LIBS += -L../tools/fpzip/lib -lfpzip
ifneq "$(GSL_ROOT_DIR)" ""
	CPPFLAGS += -I$(GSL_ROOT_DIR)/include
	LIBS += -L$(GSL_ROOT_DIR)/lib
endif
LIBS += -lgsl -lgslcblas

#################################################
ifneq ($(ERROR),)
$(error $(ERROR))  # Easiest way I found to avoid unindenting the code above...
endif

EXEOBJ = $(OBJECTS) main.o
LIBOBJ = $(OBJECTS) cubism_main.o CubismUP_3D.o
VPATH := $(VPATH):$(SRC_DIR):$(CUBISM_DIR)
ALLOBJ := $(sort $(EXEOBJ) $(LIBOBJ) $(NVOBJECTS))  # `sort` removes duplicates.
DEPS := $(ALLOBJ:%.o=%.d)

.DEFAULT_GOAL := ../bin/simulation

../bin/simulation: $(EXEOBJ) $(NVOBJECTS)
	mkdir -p ../bin
	$(LD) $(EXEOBJ) $(NVOBJECTS) $(LDFLAGS) $(LIBS) -o $@

-include $(DEPS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@


clean:
	rm -f $(EXEOBJ) $(LIBOBJ)
	rm -f $(DEPS) *.d *.o
	rm -f PoissonSolver*.o PoissonSolver*.d
	rm -f ../bin/simulation ../lib/libcubismup3d.a rlHIT PoissonSolverScalar*.o
	rmdir ../bin 2> /dev/null || true
	rmdir ../lib 2> /dev/null || true

ifeq "$(config)" "prod"
else
	#RED='\033[0;31m'
	#NC='\033[0m' # No Color
	$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif
