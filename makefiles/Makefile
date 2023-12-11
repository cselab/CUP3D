config ?= production
precision ?= double
bs ?= 8
gpu ?= false
openmp ?= false
profile ?= false
onetbb ?= false
symmetry ?= false

# SET FLAGS FOR COMPILER
ifneq ($(MPICXX),)
	CXX=$(MPICXX)
else
	CXX=mpic++
endif

# LOAD FLAGS -- GCC NEEDED
CPPFLAGS+= -std=c++17 -Wall -g
CPPFLAGS+= -Wextra -Wfloat-equal -Wcast-align -Woverloaded-virtual
CPPFLAGS+= -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow
CPPFLAGS+= -Wwrite-strings -Wno-unused-parameter
CPPFLAGS+= -Wno-float-equal
CPPFLAGS+= -Wno-redundant-decls

ifeq "$(openmp)" "true"
	CPPFLAGS+= -fopenmp
endif

# FLAGS FOR EXTERNAL LIBRARIES
LIBS+= -lgsl -lgslcblas -fopenmp -lhdf5

# ADD LIBRARY PATHS IF GIVEN
ifneq ($(HDF5_ROOT),)
	LIBS     += -L$(HDF5_ROOT)/lib
	CPPFLAGS += -I$(HDF5_ROOT)/include
endif

ifneq ($(GSL_ROOT),)
	CPPFLAGS += -I$(GSL_ROOT)/include
	LIBS += -L$(GSL_ROOT)/lib
endif

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

# ENABLE OPTIMIZATION/DEBUG FLAGS IF WISHED
ifeq "$(findstring prod,$(config))" ""
	CPPFLAGS+= -O0
	ifeq "$(config)" "segf"
		CPPFLAGS+= -fsanitize=address
		LIBS+= -fsanitize=address -static-libasan
	endif
	ifeq "$(config)" "nans"
		CPPFLAGS+= -fsanitize=undefined
		LIBS+= -fsanitize=undefined
	endif
else
	CPPFLAGS+= -DNDEBUG -O3 -fstrict-aliasing -march=native -mtune=native -falign-functions -ftree-vectorize -fmerge-all-constants
endif

# SET FLOATING POINT ACCURACY
ifeq "$(precision)" "single"
	CPPFLAGS += -D_FLOAT_PRECISION_
else ifeq "$(precision)" "double"
	CPPFLAGS += -D_DOUBLE_PRECISION_
else ifeq "$(precision)" "long_double"
	CPPFLAGS += -D_LONG_DOUBLE_PRECISION_
endif

# SET VPATH FOR MAKE TO SEARCH FOR FILES
BUILDDIR = .
DIRS = $(sort $(dir $(wildcard ../source/*) $(wildcard ../source/*/)))
VPATH := $(DIRS) $(BUILDDIR)/../Cubism/src/

ifeq "$(symmetry)" "true"
	CPPFLAGS += -DPRESERVE_SYMMETRY
else ifneq "$(findstring prod,$(config))" ""
	CPPFLAGS+= -ffast-math
endif

weno ?= false
ifeq "$(weno)" "true"
	CPPFLAGS += -DWENO
endif

# SET FLAGS FOR CUBISM
CPPFLAGS += -DCUP_ALIGNMENT=64 -DCUP_BLOCK_SIZEX=$(bs) -DCUP_BLOCK_SIZEY=$(bs) -DCUP_BLOCK_SIZEZ=$(bs) -DCUP_NO_MACROS_HEADER
CPPFLAGS += -I$(BUILDDIR)/../Cubism/include/ -DDIMENSION=3

OBJECTS = ObstacleFactory.o Obstacle.o Ellipsoid.o Cylinder.o Pipe.o \
	Fish.o StefanFish.o CarlingFish.o Sphere.o Plate.o ExternalObstacle.o Naca.o \
	FishLibrary.o BufferedLogger.o SimulationData.o Simulation.o \
	AdvectionDiffusion.o AdvectionDiffusionImplicit.o ComputeDissipation.o \
	PressureProjection.o Penalization.o InitialConditions.o FluidSolidForces.o \
	ObstaclesCreate.o ObstaclesUpdate.o FishShapes.o ArgumentParser.o \
	PoissonSolverBase.o PoissonSolverAMR.o PoissonSolverAMRKernels.o DiffusionSolverAMRKernels.o \
	FixMassFlux.o ExternalForcing.o SmartNaca.o CylinderNozzle.o

#################################################
# CUDA
#################################################
NVCC ?= nvcc
NVCCFLAGS ?= -code=sm_60 -arch=compute_60
ifeq ("$(gpu)","true")
	OBJECTS += PoissonSolverExp.o LocalSpMatDnVec.o BiCGSTAB.o
	CPPFLAGS += -fopenmp -DGPU_POISSON -Wno-shadow -Wno-undef -Wno-float-equal -Wno-redundant-decls
	NVCCFLAGS += -std=c++17 -O3 --use_fast_math -Xcompiler "$(CPPFLAGS)" -DGPU_POISSON
	LIBS += -lcudart -lcublas -lcusparse
	ifeq ("$(profile)", "true")
		NVCCFLAGS += -DBICGSTAB_PROFILER
	endif
else
	CPPFLAGS += -Wno-unknown-pragmas
endif

# DEFINE COMPILATION TARGETS
all: ../bin/simulation libcup.a cup.cflags.txt cup.libs.txt
.DEFAULT : all


# COMPILATION INSTRUCTIONS FOR APPLICATION AND LIBRARY
../bin/simulation: main.o $(OBJECTS)
	mkdir -p ../bin
	$(CXX) main.o $(OBJECTS) $(LIBS) -o $@
libcup.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)
cup.cflags.txt:
	echo '$(CPPFLAGS)' > cup.cflags.txt
cup.libs.txt:
	echo '$(LIBS)' > cup.libs.txt

# COMPILATION INSTRUCTIONS FOR OBJECT FILES
%.o: %.cu
	$(NVCC) -ccbin=$(CXX) $(NVCCFLAGS) -c $< -o $@
%.d: %.cu
	$(NVCC) -ccbin=$(CXX) $(NVCCFLAGS) -c -MD $<
%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@
%.d: %.cpp
	$(CXX) $(CPPFLAGS) -c -MD $<

# COMPILATION INSTRUCTION FOR CLEANING BUILD
clean:
	rm -f *.d *.o ../bin/simulation ../lib/libcubismup3d.a cup.cflags.txt cup.libs.txt libcup.a
	rmdir ../bin 2> /dev/null || true
	rmdir ../lib 2> /dev/null || true
