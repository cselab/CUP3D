#CC ?= g++
CC ?= mpi++
LD=$(CC)

nthreads ?= 48
config ?= debug
bc ?= mixed
precision ?= double
poisson ?= split-fftw
multiphase ?= false
bs ?= 32
bsx ?= $(bs)
bsy ?= $(bs)
bsz ?= $(bs)
vertexcentered ?= false
dlm ?= false

hdf ?= true
vtk ?= false

constnu ?= true
densitydiff ?= false

movingframe ?= false

rk2 ?= true

CPPFLAGS+= -w -D_ALIGNBYTES_=16

CPPFLAGS+= -Wunused

ifeq "$(bc)" "periodic"
	CPPFLAGS += -D_PERIODIC_
endif
ifeq "$(bc)" "pipe"
	CPPFLAGS += -D_PIPE_
endif
ifeq "$(bc)" "vortex"
	CPPFLAGS += -D_VORTEX_
endif
ifeq "$(bc)" "mixed"
	CPPFLAGS += -D_MIXED_
endif

ifeq "$(precision)" "single"
	CPPFLAGS += -D_SP_COMP_
endif

CPPFLAGS+= -D_BS_=$(bs)
CPPFLAGS+= -D_BSX_=$(bsx)
CPPFLAGS+= -D_BSY_=$(bsy)
CPPFLAGS+= -D_BSZ_=$(bsz)

ifeq "$(config)" "production"
	CPPFLAGS += -DNDEBUG -O3
else
	CPPFLAGS += -g
endif

ifeq "$(poisson)" "split-fftw"
	CPPFLAGS += -D_SPLIT_
endif
ifeq "$(poisson)" "hypre"
	CPPFLAGS += -D_MULTIGRID_
endif

ifeq "$(multiphase)" "true"
	CPPFLAGS += -D_MULTIPHASE_
endif

ifeq "$(vertexcentered)" "true"
	CPPFLAGS += -D_VERTEXCENTERED_
endif

ifeq "$(dlm)" "true"
	CPPFLAGS += -D_DLM_
endif

ifeq "$(hdf)" "true"
	CPPFLAGS += -D_USE_HDF_
	LIBS += -lhdf5
endif

ifeq "$(constnu)" "true"
	CPPFLAGS += -D_CONSTNU_
endif

ifeq "$(densitydiff)" "true"
	CPPFLAGS += -D_DENSITYDIFF_
endif

ifeq "$(rk2)" "true"
	CPPFLAGS += -D_RK2_
endif

ifeq "$(movingframe)" "true"
	CPPFLAGS += -D_MOVING_FRAME_
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

CPPFLAGS += -std=c++11 -fopenmp -fpermissive
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
OBJECTS = GeometryReader.o WavefrontReader.o CTReader.o PoissonSolverScalarFFTW.o PoissonSolverScalarFFTW_MPI.o ProcessOperatorsOMP.o Sim_Jet.o Sim_RayleighTaylor.o Sim_Bubble.o Sim_FSI_Gravity.o Sim_FSI_Rotation.o Sim_FSI_Moving.o Sim_FSI_Fixed.o Profiler.o main.o
OBJECTS_TEST = GeometryReader.o WavefrontReader.o CTReader.o PoissonSolverScalarFFTW.o PoissonSolverScalarFFTW_MPI.o ProcessOperatorsOMP.o TestMPI.o TestBoundaryConditions.o TestAddedMass.o TestPoiseuille.o TestShearLayer.o TestTravelingWave.o TestGeometry.o TestRotation.o TestTranslation.o TestPenalization.o TestGravity.o TestVarCoeffPoisson.o TestPressure.o TestAdvection.o TestDiffusion.o Profiler.o mainTest.o

VPATH := $(SRC_DIR):$(CUBISM_DIR)

.DEFAULT: all;

all: simulation test

simulation: $(OBJECTS)
#	echo Linking
#	echo OBJECTS = $(OBJECTS)
#	echo ************************
	$(LD) $(OBJECTS) -o $@ $(LIBS) $(LDFLAGS)
#	echo ************************

test: $(OBJECTS_TEST)
#	echo Linking
#	echo OBJECTS = $(OBJECTS_TEST)
#	echo ************************
	$(LD) $(OBJECTS_TEST) -o $@ $(LIBS) $(LDFLAGS)
#	echo ************************

%.o: $(SRC_DIR)/%.cpp
#	echo PRODUCING $@
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f simulation test
	rm -f $(OBJECTS) $(OBJECTS_TEST)
