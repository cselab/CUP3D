# module load hdf5-parallel/1.8.11
CC = cc
#CC = CC
bs = 32
#hdf = 1
include ../Makefile.config 

CPPFLAGS += -I../reader 
CPPFLAGS += $(extra)

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
VPATH := ../../Cubism/ ../../source/ .

all: ch2hdf ch2roi

ch2hdf: main.o WaveletCompressor.o
	$(CC) $(CPPFLAGS) WaveletCompressor.o  main.o  -o ch2hdf $(LIBS)

main.o: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o main.o
	
ch2roi: main_roi.o WaveletCompressor.o
	$(CC) $(CPPFLAGS) WaveletCompressor.o  main_roi.o  -o ch2roi $(LIBS)

main_roi.o: main_roi.cpp
	$(CC) $(CPPFLAGS) -c main_roi.cpp -o main_roi.o

WaveletCompressor.o: ../../source/WaveletCompressor.cpp
	$(CC) $(CPPFLAGS) -c ../../source/WaveletCompressor.cpp -o WaveletCompressor.o
	
clean:
	rm -f ch2hdf ch2roi *.o *.lst

cleanall:
	rm -f ch2hdf ch2roi *.o *.lst
