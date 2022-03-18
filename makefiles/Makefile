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
config ?= prod
precision ?= double
precision_dump ?= single
usemap ?= true
bs ?= 8
bsx ?= 8
bsy ?= 8
bsz ?= 8
hdf ?= true
verbose ?= false
dump_surface_binary ?= false

precision ?= double
ifeq "$(precision)" "single"
CPPFLAGS += -D_FLOAT_PRECISION_
CPPFLAGS += -Wno-narrowing
endif
ifeq "$(precision)" "double"
CPPFLAGS += -D_DOUBLE_PRECISION_
endif
ifeq "$(precision)" "long_double"
CPPFLAGS += -D_LONG_DOUBLE_PRECISION_
endif


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
	ObstaclesCreate.o ObstaclesUpdate.o \
	FishShapes.o ArgumentParser.o PoissonSolverAMR.o

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

all: ../bin/simulation libcup.a cup.cflags.txt cup.libs.txt
.DEFAULT : all

# COMPILATION INSTRUCTIONS FOR APPLICATION AND LIBRARY
../bin/simulation: $(EXEOBJ) $(NVOBJECTS)
	mkdir -p ../bin
	$(LD) $(EXEOBJ) $(NVOBJECTS) $(LDFLAGS) $(LIBS) -o $@

libcup.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

cup.cflags.txt:
	echo $(CPPFLAGS) > cup.cflags.txt

cup.libs.txt:
	echo $(LIBS) > cup.libs.txt

# COMPILATION INSTRUCTIONS FOR OBJECT FILES
%.o: %.cpp
	$(CXX) $(CPPFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

# COMPILATION INSTRUCTION FOR CLEANING BUILD
clean:
	rm -f $(EXEOBJ) $(LIBOBJ)
	rm -f $(DEPS) *.d *.o
	rm -f PoissonSolver*.o PoissonSolver*.d
	rm -f ../bin/simulation ../lib/libcubismup3d.a rlHIT PoissonSolverScalar*.o
	rm -f cup.cflags.txt cup.libs.txt libcup.a
	rmdir ../bin 2> /dev/null || true
	rmdir ../lib 2> /dev/null || true

# WARN USER FOR NOT USING PRODUCTION FLAGS
ifeq "$(config)" "prod"
else
	#RED='\033[0;31m'
	#NC='\033[0m' # No Color
	$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif
