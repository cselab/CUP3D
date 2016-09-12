#CC ?= g++
CC ?= mpi++
LD=$(CC)

nthreads ?= 48
config ?= debug
precision ?= double
poisson ?= split-fftw
bs ?= 32
bsx ?= $(bs)
bsy ?= $(bs)
bsz ?= $(bs)
hdf ?= true
vtk ?= false
rk2 ?= true

CPPFLAGS+= -w -D_ALIGNBYTES_=16 -Wunused

ifeq "$(precision)" "single"
	CPPFLAGS += -D_SP_COMP_
	ISPCFLAGS += -D_SP_COMP_ --target=avx1-i32x8
else
	ISPCFLAGS += --target=avx1-i64x4
endif

CPPFLAGS += -D_BS_=$(bs)
CPPFLAGS += -D_BSX_=$(bsx)
CPPFLAGS += -D_BSY_=$(bsy)
CPPFLAGS += -D_BSZ_=$(bsz)

ifeq "$(config)" "production"
	CPPFLAGS += -DNDEBUG -O3 -march=native
	ISPCFLAGS += -O3
else
	CPPFLAGS += -g
endif

ifeq "$(poisson)" "split-fftw"
	CPPFLAGS += -D_SPLIT_
endif

ifeq "$(hdf)" "true"
	CPPFLAGS += -D_USE_HDF_
	LIBS += -lhdf5
endif

ifeq "$(rk2)" "true"
	CPPFLAGS += -D_RK2_
endif

vtk-inc ?=
vtk-lib ?=
fftw-inc ?=
fftw-lib ?=
hypre-inc ?=
hypre-lib ?=

# needed in order to compile on nodes
include make.daint

ifneq "$(findstring brutus,$(shell hostname))" ""
	include make.brutus
endif

ifneq "$(findstring euler,$(shell hostname))" ""
	include make.euler
endif

CPPFLAGS += -std=c++11 -fopenmp #-fpermissive
CPPFLAGS += -DNTHREADS=$(nthreads)
CPPFLAGS += -I../source/ -I../Cubism/ -I$(fftw-inc)
CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

LDFLAGS = $(CPPFLAGS)

LIBS += -lstdc++ -lm -lz

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

ifeq "$(precision)" "single"
	LIBS += -L$(fftw-lib) -lfftw3 -lfftw3f -lfftw3f_threads -lfftw3_mpi
endif
ifeq "$(precision)" "double"
	LIBS += -L$(fftw-lib) -lfftw3 -lfftw3_threads -lfftw3_mpi
endif

LIBS += -L../tools/fpzip/lib -lfpzip

ifeq "$(poisson)" "hypre"
CC=mpic++
LD=$(CC)
CPPFLAGS += -I$(hypre-inc)
LIBS += -L$(hypre-lib) -lHYPRE
endif

BUILDDIR = .
SRC_DIR = $(BUILDDIR)/../source/
CUBISM_DIR = $(BUILDDIR)/../Cubism/

OBJECTS = GeometryReader.o WavefrontReader.o CTReader.o WaveletCompressor.o PoissonSolverScalarFFTW.o PoissonSolverScalarFFTW_MPI.o ProcessOperatorsOMP.o Simulation.o IF2D_ObstacleFactory.o IF3D_ObstacleOperator.o IF3D_ObstacleVector.o IF3D_CarlingFishOperator.o IF3D_ForcedSphereObstacleOperator.o IF3D_SphereObstacleOperator.o Profiler.o main.o

VPATH := $(SRC_DIR):$(CUBISM_DIR)

.DEFAULT: all;

all: simulation

simulation: $(OBJECTS)
#	echo Linking
#	echo OBJECTS = $(OBJECTS)
#	echo ************************
	$(LD) $(OBJECTS) -o $@ $(LIBS) $(LDFLAGS)
#	echo ************************

clean:
	rm -f simulation
	rm -f $(OBJECTS) $(OBJECTS_TEST)
