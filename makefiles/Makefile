ifeq ($(shell uname -s), Darwin)
CXX=mpic++
GSL_ROOT_DIR=$(wildcard /usr/local/Cellar/gsl/*)
HDF5_ROOT=$(wildcard /usr/local/Cellar/hdf5/*)
endif

ifneq "$(findstring nid,$(shell hostname))" ""
CXX=CC
endif
ifneq "$(findstring daint,$(shell hostname))" ""
CXX=CC
endif
ifneq "$(findstring panda,$(shell hostname))" ""
CXX=mpic++
endif
ifneq "$(findstring falcon,$(shell hostname))" ""
CXX=mpic++
endif
ifneq "$(findstring euler,$(shell hostname))" ""
CXX=mpic++
endif
ifneq "$(findstring eu-,$(shell hostname))" ""
CXX=mpic++
endif

# this only affects amount of memory allocated by block processing:
config ?= prod
precision ?= double
usemap ?= true
bs ?= 8
bsx ?= 8
bsy ?= 8
bsz ?= 8
hdf ?= true
verbose ?= false
dump_surface_binary ?= false
gpu ?= false
profile ?= false
onetbb ?= false

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

weno ?= false
ifeq "$(weno)" "true"
	CPPFLAGS += -DWENO
endif
symmetry ?= false
ifeq "$(symmetry)" "true"
	CPPFLAGS += -DPRESERVE_SYMMETRY
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

LD=$(CXX)

ifeq "$(usemap)" "true"
	CPPFLAGS += -DCUBISM_USE_MAP
endif

ifeq "$(precision)" "single"
	CPPFLAGS += -DCUP_SINGLE_PRECISION
	NVFLAGS += -DCUP_SINGLE_PRECISION
endif

CPPFLAGS += -DCUP_ALIGNMENT=64
#CPPFLAGS += -DCUP_BLOCK_SIZE=$(bs)
CPPFLAGS += -DCUP_BLOCK_SIZEX=$(bsx)
CPPFLAGS += -DCUP_BLOCK_SIZEY=$(bsy)
CPPFLAGS += -DCUP_BLOCK_SIZEZ=$(bsz)
CPPFLAGS += -DCUP_NO_MACROS_HEADER
NVFLAGS  += -DCUP_NO_MACROS_HEADER
CPPFLAGS += -I../Cubism/include/

LDFLAGS = $(CPPFLAGS)

BUILDDIR = .
#SRC_DIR = $(BUILDDIR)/../source/
SRC_DIR = $(sort $(dir $(wildcard $(BUILDDIR)/../source/*/)))
CUBISM_DIR = $(BUILDDIR)/../Cubism/src/

OBJECTS = ObstacleFactory.o Obstacle.o Ellipsoid.o Cylinder.o \
	Fish.o StefanFish.o CarlingFish.o Sphere.o Plate.o ExternalObstacle.o Naca.o \
	FishLibrary.o BufferedLogger.o SimulationData.o Simulation.o \
	AdvectionDiffusion.o AdvectionDiffusionImplicit.o ComputeDissipation.o \
	PressureProjection.o Penalization.o InitialConditions.o FluidSolidForces.o \
	ObstaclesCreate.o ObstaclesUpdate.o FishShapes.o ArgumentParser.o \
	PoissonSolverBase.o PoissonSolverAMR.o PoissonSolverAMRKernels.o


#################################################
# oneTBB
#################################################
ifeq "$(onetbb)" "true"
	LIBS     += -L$(ONETBBROOT)/lib64 -ltbb
	CPPFLAGS += -I$(ONETBBROOT)/include
	CPPFLAGS += -DCUBISM_USE_ONETBB
	#LIBS     += -L/users/chatzima/my_installed_onetbb/lib64 -ltbb
	#CPPFLAGS += -I/users/chatzima/my_installed_onetbb/include
endif

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
endif

LIBS += -lgsl -lgslcblas

#################################################
# CUDA
#################################################
NVCC ?= nvcc
NVFLAGS ?= -code=sm_60 -arch=compute_60

ifeq "$(gpu)" "true"
	OBJECTS += PoissonSolverExp.o LocalSpMatDnVec.o BiCGSTAB.o
	CPPFLAGS += -fopenmp -DGPU_POISSON -Wno-shadow -Wno-undef -Wno-float-equal -Wno-redundant-decls
	NVFLAGS = -std=c++17 -O3 --use_fast_math -Xcompiler "$(CPPFLAGS)" -DGPU_POISSON
	LIBS += -lcudart -lcublas -lcusparse
	
	ifeq "$(profile)" "true"
		NVFLAGS += -DBICGSTAB_PROFILER
	endif
endif

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
%.o: %.cu
	$(NVCC) -ccbin=$(CXX) $(NVFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -MMD -MF $(patsubst %.o,%.d,$@) -c $< -o $@

# COMPILATION INSTRUCTION FOR CLEANING BUILD
clean:
	rm -f *.d *.o ../bin/simulation ../lib/libcubismup3d.a cup.cflags.txt cup.libs.txt libcup.a
	rmdir ../bin 2> /dev/null || true
	rmdir ../lib 2> /dev/null || true

# WARN USER FOR NOT USING PRODUCTION FLAGS
ifeq "$(config)" "prod"
else
	#RED='\033[0;31m'
	#NC='\033[0m' # No Color
	$(info *** WARNING! Compiling with debug flags. Use config=prod for optimized code ***)
endif
