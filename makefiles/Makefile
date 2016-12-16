#CC ?= g++
CC ?= mpic++
LD=$(CC)

nthreads ?= 48
config ?= debug
precision ?= double
bs ?= 16
bsx ?= $(bs)
bsy ?= $(bs)
bsz ?= $(bs)
hdf ?= true
vtk ?= false
rk2 ?= true

CPPFLAGS+= -w -D_ALIGNBYTES_=16 -Wunused

ifeq "$(precision)" "single"
	CPPFLAGS += -D_SP_COMP_
endif

CPPFLAGS += -D_BS_=$(bs)
CPPFLAGS += -D_BSX_=$(bsx)
CPPFLAGS += -D_BSY_=$(bsy)
CPPFLAGS += -D_BSZ_=$(bsz)

ifeq "$(config)" "prod"
	CPPFLAGS += -DNDEBUG -O3 -march=native
else
	CPPFLAGS += -g
endif

ifeq "$(rk2)" "true"
	CPPFLAGS += -D_RK2_
endif

vtk-inc ?=
vtk-lib ?=
fftw-inc ?=
fftw-lib ?=
hdf-inc ?=
hdf-lib ?=
lz4-inc ?=
lz4-lib ?=

# needed in order to compile on nodes
# include make.daint

ifneq "$(findstring daint,$(shell hostname))" ""
	include make.daint
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
	include make.brutus
endif

ifneq "$(findstring euler,$(shell hostname))" ""
	include make.euler
endif

CPPFLAGS += -std=c++11 -fopenmp -fpermissive
CPPFLAGS += -DNTHREADS=$(nthreads)
CPPFLAGS += -I../source/ -I../Cubism/ #-I$(fftw-inc)
CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

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
	CPPFLAGS += -D_USE_HDF_ -I$(hdf-inc)
	LIBS += -lhdf5 #-L$(hdf-lib) -lhdf5
else
	CPPFLAGS += -D_USE_LZ4_ #-I$(lz4-inc) -I$(hdf-inc)
	LIBS += #-L$(lz4-lib) -L$(hdf-lib) -lhdf5
endif

ifeq "$(precision)" "single"
	LIBS += -lfftw3f -lfftw3f_threads -lfftw3f_mpi #-L$(fftw-lib) -lfftw3f -lfftw3f_threads -lfftw3f_mpi
endif
ifeq "$(precision)" "double"
	LIBS += -lfftw3 -lfftw3_threads -lfftw3_mpi # -L$(fftw-lib) -lfftw3 -lfftw3_threads -lfftw3_mpi
endif 

LIBS += -L../tools/fpzip/lib -lfpzip

LIBS += -lstdc++ -lm -ldl -lz

BUILDDIR = .
SRC_DIR = $(BUILDDIR)/../source/
CUBISM_DIR = $(BUILDDIR)/../Cubism/

OBJECTS = PoissonSolverScalarFFTW.o PoissonSolverScalarFFTW_MPI.o ProcessOperatorsOMP.o IF3D_ObstacleFactory.o IF3D_ObstacleOperator.o IF3D_ObstacleVector.o IF3D_CarlingFishOperator.o IF3D_ForcedSphereObstacleOperator.o IF3D_SphereObstacleOperator.o Profiler.o Simulation.o main.o
#GeometryReader.o WavefrontReader.o CTReader.o WaveletCompressor.o
VPATH := $(SRC_DIR):$(CUBISM_DIR)

.DEFAULT_GOAL := simulation

simulation: $(OBJECTS)
	#echo Linking
	#echo CC = $(CC) OBJECTS = $(OBJECTS)
	#echo ************************
	$(LD) $(OBJECTS) -o $@ $(LIBS) $(LDFLAGS)
	#echo ************************
    
%.o: %.cpp
	$(CC) $(CPPFLAGS) $(LIBS) -c $< -o $@

clean:
	rm -f $(OBJECTS)
	rm -f simulation
