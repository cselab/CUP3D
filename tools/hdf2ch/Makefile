SHELL := /bin/bash

CC ?= cc
bs ?= 32
ap ?= float
hdf ?= 1
config ?=release

include ../Makefile.config

# CSCS
CPPFLAGS += -I${HDF5_INCLUDE_OPTS} -D_BSX_=$(bs) -D_BSY_=$(bs) -D_BSZ_=$(bs)

#LIBS = $(LDFLAGS)
#ifeq "$(zlib)" "1"
#	MYFLAGS += -D_USE_ZLIB_
#endif
#
#ifeq "$(lz4)" "1"
#	MYFLAGS += -D_USE_LZ4_
#endif


CPPFLAGS += -I../../Cubism/ -I../../source/ -I.  -I../reader

##################
VPATH := ../../Cubism/ ../../source/ .
#.DEFAULT_GOAL := hdf2ch

OBJECTS = main.o WaveletCompressor.o

all: hdf2ch

hdf2ch: $(OBJECTS) 
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

WaveletCompressor.o: ../../source/WaveletCompressor.cpp
	$(CC) $(CPPFLAGS) -c ../../source/WaveletCompressor.cpp -o WaveletCompressor.o

%.o: %.cpp 
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o hdf2ch
