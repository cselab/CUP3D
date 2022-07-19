CXX=CC
LD=CC
CPPFLAGS+= -std=c++17 -g -fopenmp -Wall
LIBS+= -fopenmp
CPPFLAGS+= -DNDEBUG -O3 -fstrict-aliasing -march=native -mtune=native -ffast-math -falign-functions -ftree-vectorize -fmerge-all-constants


OBJECTS = main.o

#################################################
# HDF5
#################################################
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
#################################################


convert: $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(OBJECTS) $(LIBS) -o $@
	rm main.o


clean:
	rm -f $(OBJECTS) convert
