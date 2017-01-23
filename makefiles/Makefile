#CC ?= g++
CC ?= mpic++
NVCC?=nvcc
NVFLAGS?=-O3 -code=sm_60 -arch=compute_60
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
rk2 ?= false
accfft ?= true
zlib ?= 0

ifeq "$(zlib)" "1"
	CPPFLAGS += -D_USE_ZLIB_ -D_USE_WAVZ_
	LIBS +=  -lz -lrt -lstdc++
endif

CPPFLAGS+= -w -D_ALIGNBYTES_=16 -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align
CPPFLAGS+= -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls
CPPFLAGS+= -Wshadow -Woverloaded-virtual

ifeq "$(precision)" "single"
	CPPFLAGS += -D_SP_COMP_
	NVCCFLAGS += -D_SP_COMP_
endif

CPPFLAGS += -D_BLOCKSIZE_=$(bs)
CPPFLAGS += -D_BS_=$(bs)
CPPFLAGS += -D_BSX_=$(bsx)
CPPFLAGS += -D_BSY_=$(bsy)
CPPFLAGS += -D_BSZ_=$(bsz)
CPPFLAGS += -I${HDF5_INCLUDE_OPTS}

ifeq "$(config)" "prod"
	CPPFLAGS += -DNDEBUG -O4 -fstrict-aliasing -march=native -mtune=native -ffast-math -falign-functions=16
else
	CPPFLAGS += -g
endif

ifeq "$(rk2)" "true"
	CPPFLAGS += -D_RK2_
endif

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

CPPFLAGS += -std=c++11 -fopenmp #-fpermissive
#NVFLAGS += -std=c++11
CPPFLAGS += -DNTHREADS=$(nthreads)
CPPFLAGS += -I../source/ -I../Cubism/
#CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

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

CPPFLAGS += -D_USE_HDF_ #-I$(hdf-inc)
LIBS += -lhdf5 #-L$(hdf-lib) -lhdf5

ifeq "$(precision)" "single"
	LIBS += -lfftw3f -lfftw3f_threads -lfftw3f_mpi
endif
ifeq "$(precision)" "double"
	LIBS += -lfftw3 -lfftw3_threads -lfftw3_mpi
endif

#LIBS += -L../tools/fpzip/lib -lfpzip

LIBS += -lstdc++ -lgsl -lgslcblas -lm -ldl -lz -lpthread

BUILDDIR = .
SRC_DIR = $(BUILDDIR)/../source/
CUBISM_DIR = $(BUILDDIR)/../Cubism/

OBJECTS = Communicator.o IF3D_ObstacleFactory.o IF3D_ObstacleOperator.o IF3D_ObstacleVector.o IF3D_StefanFishOperator.o IF3D_DeadFishOperator.o IF3D_CarlingFishOperator.o IF3D_FishOperator.o IF3D_ForcedSphereObstacleOperator.o IF3D_SphereObstacleOperator.o Profiler.o Simulation.o Save_splicer.o WaveletCompressor.o main.o
#GeometryReader.o WavefrontReader.o CTReader.o WaveletCompressor.o

ifeq "$(accfft)" "true"
	CPPFLAGS += -D_ACCFFT_ -I$(accfft-inc) -D_CUDA_COMP_
	LIBS += -L$(accfft-lib) -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -lcufft
	NVOBJECTS = PoissonSolverScalarFFTW_ACC.o
else
	NVOBJECTS =
	OBJECTS += PoissonSolverScalarFFTW.o PoissonSolverScalarFFTW_MPI.o
endif

VPATH := $(SRC_DIR):$(CUBISM_DIR)

.DEFAULT_GOAL := simulation

simulation: $(OBJECTS) $(NVOBJECTS)
	#echo Linking
	#echo CC = $(CC) OBJECTS = $(OBJECTS)
	#echo ************************
	$(LD) $(OBJECTS) $(NVOBJECTS) $(LIBS) $(LDFLAGS) -o $@
	#echo ************************

%.o: %.cpp
	$(CC) $(CPPFLAGS) $(LIBS) -c $< -o $@

PoissonSolverScalarFFTW_ACC.o: PoissonSolverScalarFFTW_ACC.cu
	nvcc $(NVFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS)
	rm -f simulation
	rm -f PoissonSolverScalarFFTW*.o
