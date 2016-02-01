/*
 *  Test_IO_Compressed.h
 *  CUBISMTests
 *
 *  Created by Babak Hejazialhosseini on 5/16/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Types_single.h"
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#include "Reader_WaveletCompression.h"
#include <ArgumentParser.h>
#include <GridMPI.h>
#define VERBOSE 0

#ifdef _USE_HDF_
#include <hdf5.h>
#endif

#ifdef _FLOAT_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
#else
#define HDF_REAL H5T_NATIVE_DOUBLE
#endif


typedef GridMPI < FluidGrid > G;


class Test_IO_Compressed : public Simulation
{
protected:

	int BPDX, BPDY, BPDZ;
	int XPESIZE, YPESIZE, ZPESIZE;
	int channel;
	int step_id;
	bool VERBOSITY;
	G * grid;
	ArgumentParser parser;
	string inputfile_name, outputfile_name;
	int myrank;
    
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<G, StreamerGridPointIterative> mywaveletdumper;

	void _ic(G& grid)
	{
					
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
//grid.getBlocksInfo();
		Real min_u = 1e8, max_u = -1e8;

#if 1
		int myrank, mypeindex[3], pesize[3];
		bool periodic[3];
		MPI::Cartcomm cartcomm;

		periodic[0] = true;
		periodic[1] = true;
		periodic[2] = true;

		pesize[0] = XPESIZE;
		pesize[1] = YPESIZE;
		pesize[2] = ZPESIZE;

		assert(XPESIZE*YPESIZE*ZPESIZE == MPI::COMM_WORLD.Get_size());

		cartcomm = MPI::COMM_WORLD.Create_cart(3, pesize, periodic, true);
		myrank = cartcomm.Get_rank();

		cartcomm.Get_coords(myrank, 3, mypeindex);
#endif

#if 1
		// open HDF5
		// read data block
		// store in grid block

		int rank;
		char filename[256];
		herr_t status;
		hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int coords[3];
		coords[0] = mypeindex[0]; coords[1] = mypeindex[1]; coords[2] = mypeindex[2];
		//grid.peindex(coords);

		const int NX = BPDX*_BLOCKSIZE_;
		const int NY = BPDY*_BLOCKSIZE_;
		const int NZ = BPDZ*_BLOCKSIZE_;

		printf("coords = %d,%d,%d\n", coords[0], coords[1], coords[2]);
		printf("NX,NY,NZ = %d,%d,%d\n", NX, NY, NZ);

		static const int NCHANNELS = 9;

		Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

		vector<BlockInfo> vInfo_local = grid.getBlocksInfo();

		static const int sX = 0;
		static const int sY = 0;
		static const int sZ = 0;

		const int eX = _BLOCKSIZE_;
		const int eY = _BLOCKSIZE_;
		const int eZ = _BLOCKSIZE_;

		hsize_t count[4] = {
			NX,
			NY,
			NZ, NCHANNELS};

		hsize_t dims[4] = {
			grid.getBlocksPerDimension(0)*_BLOCKSIZE_,
			grid.getBlocksPerDimension(1)*_BLOCKSIZE_,
			grid.getBlocksPerDimension(2)*_BLOCKSIZE_, NCHANNELS};

		hsize_t offset[4] = {
			coords[0]*NX,
			coords[1]*NY,
			coords[2]*NZ, 0};

		sprintf(filename, "%s", inputfile_name.c_str());

		H5open();
		fapl_id = H5Pcreate(H5P_FILE_ACCESS);
		status = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
		file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
		status = H5Pclose(fapl_id);

		dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
		fapl_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

		fspace_id = H5Dget_space(dataset_id);
		H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

		mspace_id = H5Screate_simple(4, count, NULL);
		status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

		printf("before barrier\n"); fflush(0);
		MPI_Barrier(MPI_COMM_WORLD);
		printf("after barrier\n"); fflush(0);
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		//#pragma omp parallel for
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			const int idx[3] = {info.index[0], info.index[1], info.index[2]};
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			//Streamer streamer(b);

			//printf("setting info.index=[%d,%d,%d]\n", info.index[0], info.index[1], info.index[2]);

			for(int ix=sX; ix<eX; ix++)
				for(int iy=sY; iy<eY; iy++)
					for(int iz=sZ; iz<eZ; iz++)
					{
						const int gx = idx[0]*_BLOCKSIZE_ + ix;
						const int gy = idx[1]*_BLOCKSIZE_ + iy;
						const int gz = idx[2]*_BLOCKSIZE_ + iz;

						Real * const ptr_input = array_all + NCHANNELS*(gz + NZ * (gy + NY * gx));
						Real val = ptr_input[channel];

						//streamer.operate(ptr_input, ix, iy, iz);
						b(ix, iy, iz).u = val;
						if (val < min_u) min_u = val;
						if (val > max_u) max_u = val;

					}
		}


		printf("min_u = %lf max_u = %lf\n", min_u, max_u);

		status = H5Pclose(fapl_id);
		status = H5Dclose(dataset_id);
		status = H5Sclose(fspace_id);
		status = H5Sclose(mspace_id);
		status = H5Fclose(file_id);

		H5close();

		delete [] array_all;

	}


	void _setup_mpi_constants(int& xpesize, int& ypesize, int& zpesize)
	{
		xpesize = parser("-xpesize").asInt(1);
		ypesize = parser("-ypesize").asInt(1);
		zpesize = parser("-zpesize").asInt(1);
	}
    
public:
	const bool isroot;
	
	Test_IO_Compressed(const bool isroot, const int argc, const char ** argv) : isroot(isroot) , grid(NULL), parser(argc,argv) { }
    

	void setup()
	{
		parser.unset_strict_mode();

		BPDX = parser("-bpdx").asInt(1);
		BPDY = parser("-bpdy").asInt(1);
		BPDZ = parser("-bpdz").asInt(1);
		
		step_id = parser("-stepid").asInt(0);
		
		channel = parser("-channel").asInt(0);

		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

		inputfile_name = parser("-simdata").asString("none");
		outputfile_name = parser("-outdata").asString("none");

		if ((inputfile_name == "none")||(outputfile_name == "none"))
		{
			printf("usage: %s -simdata <filename> -outdata <filename> -channel <number> [-swap] [-wtype_read <type1>] [-wtype_write <type2>]\n", "tool");
			exit(1);
		}


		_setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);

		if (!isroot)
			VERBOSITY = 0;

		if (VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("////////////      TEST IO Compressed MPI     ///////////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}

		grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);

		assert(grid != NULL);

		_ic(*grid);
		vp(*grid, step_id);
	}


	void vp(G& grid, const int step_id)
	{
		if (isroot) cout << "dumping MPI VP ...\n" ;

		const string path = parser("-fpath").asString(".");
		const int wtype_write = parser("-wtype_write").asInt(1);

		std::stringstream streamer;
		streamer<<path;
		streamer<<"/";
		streamer<<outputfile_name;
		streamer.setf(ios::dec | ios::right);
		streamer.width(5);
		streamer.fill('0');
		streamer<<step_id;

		double threshold = parser("-threshold").asDouble(1e-5);

		mywaveletdumper.verbose();
		printf("setting threshold to %f\n", threshold);
		mywaveletdumper.set_threshold(threshold);
		mywaveletdumper.set_wtype_write(wtype_write);
		mywaveletdumper.Write<0>(grid, streamer.str());

		if (isroot) cout << "done" << endl;
	}

	void dispose()
	{
		if (grid!=NULL) {
			delete grid;
			grid = NULL;
		}
	}
};
