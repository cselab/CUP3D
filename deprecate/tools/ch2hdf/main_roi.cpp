/*
 *  main.cpp
 *  
 *
 *  Created by Panagiotis Chatzidoukas on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <sstream>
#include <mpi.h>
#include <hdf5.h>
#include <omp.h>
#define _TRANSPOSE_DATA_
#define _COLLECTIVE_IO_

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"

int main(int argc, const char **argv)
{
	const double init_t0 = omp_get_wtime();

	/* Initialize MPI */
	MPI::Init_thread(MPI_THREAD_SERIALIZED);

	/* MPI variables */
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;

	const int mpi_rank = mycomm.Get_rank();
	const int mpi_size = mycomm.Get_size();

	const bool isroot = !mpi_rank;

	ArgumentParser argparser(argc, argv);

	if (isroot)
		argparser.loud();
	else
		argparser.mute();

	const string inputfile_name = argparser("-simdata").asString("none");

	if ((inputfile_name == "none"))
	{
		printf("usage: %s -simdata <filename> [-swap] [-wtype <wtype>] [-verb <0|1|2|3>] [-stride <stride>] [-xs <x0>] [-xe <x1>] [-ys <y0>] [-ye <y1>] [-zs <z0>] [-ze <z1>]\n", argv[0]);
		exit(1);
	}

	const int verbose = argparser("-verb").asInt(0);

	const int stride = argparser("-stride").asInt(64);	// 4^3

	const bool swapbytes = argparser.check("-swap");
	const int wtype = argparser("-wtype").asInt(1);

	int Xs = argparser("-xs").asInt(-1);
	int Xe = argparser("-xe").asInt(-1);
	int Ys = argparser("-ys").asInt(-1);
	int Ye = argparser("-ye").asInt(-1);
	int Zs = argparser("-zs").asInt(-1);
	int Ze = argparser("-ze").asInt(-1);

#if 1
	Reader_WaveletCompressionMPI myreader(mycomm, inputfile_name, swapbytes, wtype);
#else
	Reader_WaveletCompression myreader(inputfile_name, swapbytes, wtype);
#endif
	myreader.load_file();
	const double init_t1 = omp_get_wtime();

	const double it0 = omp_get_wtime(); 

	int NBX = myreader.xblocks();
	int NBY = myreader.yblocks();
	int NBZ = myreader.zblocks();

	if ((!mpi_rank) || (verbose == 1))
		fprintf(stdout, "I found in total %dx%dx%d blocks.\n", NBX, NBY, NBZ);


	if (Xs == -1) Xs = 0;
	if (Xe == -1) Xe = NBX-1;
	if (Ys == -1) Ys = 0;
	if (Ye == -1) Ye = NBY-1;
	if (Zs == -1) Zs = 0;
	if (Ze == -1) Ze = NBZ-1;

	int StartX, StartY, StartZ;
	int EndX, EndY, EndZ;
	
	StartX = NBX;
	EndX = -1;

	StartY = NBY;
	EndY = -1;
	
	StartZ = NBZ;
	EndZ = -1;
	
	int set_counter = 0;

	if (verbose == 1)
		fprintf(stdout, "Initial ROI = [%d,%d]x[%d,%d]x[%d,%d]\n", StartX, EndX, StartY, EndY, StartZ, EndZ);	
		
//	int NX = (EndX-StartX+1)*_BLOCKSIZE_;
//	int NY = (EndY-StartY+1)*_BLOCKSIZE_;
//	int NZ = (EndZ-StartZ+1)*_BLOCKSIZE_;

	static Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	
	const int nblocks = NBX*NBY*NBZ;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size; 

#if 0
	int xx[4] = {0,     0,     0,     0};
	int yy[4] = {0,     0, NBY-1, NBY-1};
	int zz[4] = {0, NBZ-1,     0, NBZ-1};

	for (int diag = 0; diag < 4; diag++)
	{
		int x0, y0, z0;
		x0 = xx[diag];
		y0 = yy[diag];
		z0 = zz[diag];

		int dx, dy, dz;	// we need variable stride per direction
		dx = (x0 == 0) ? 1: -1;
		dy = (y0 == 0) ? 1: -1;
		dz = (z0 == 0) ? 1: -1;


		int x, y, z;		

		x = x0; y = y0; z = z0;
		for (int d = 0; d < NBX; d++)	// assuming cube
		{
		
//			if (verbose == 2)
			if (!mpi_rank)
				fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z); 

			double zratio = myreader.load_block2(x, y, z, targetdata);
			
			Real bmin = targetdata[0][0][0];
			Real bmax = targetdata[0][0][0];
			
			for (int zb = 0; zb < _BLOCKSIZE_; zb++)
			{
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int xb = 0; xb < _BLOCKSIZE_; xb++)
					{
						Real d = targetdata[zb][yb][xb];
						if (d > bmax) {
							bmax = d;
						}
						if (d < bmin) {
							bmin = d;
						}
					}
			}

			if (verbose == 3)
				fprintf(stdout, "%.3lf bmin = %.3lf bmax = %.3lf\n", zratio, bmin, bmax);
				
			if (bmin < bmax)	// !=
			{
				if (set_counter == 0)
				{
					StartX = EndX = x;
					StartY = EndY = y;
					StartZ = EndZ = z;					
					set_counter++;
				}
				else
				{
					if (x < StartX) StartX = x;
					if (x > EndX)	EndX = x;
					if (y < StartY) StartY = y;
					if (y > EndY)	EndY = y;
					if (z < StartZ) StartZ = z;
					if (z > EndZ)	EndZ = z;
				}
			}
			x += dx;
			y += dy;
			z += dz;
		}
	}
#endif

#if 1
	for (int b = mpi_rank; b < b_end; b += mpi_size*stride)	
	{
#if defined(_TRANSPOSE_DATA_)
		int x = b / (NBY * NBZ);
		int y = (b / NBZ) % NBY;
		int z = b % NBZ;
#else
		int z = b / (NBY * NBX);
		int y = (b / NBX) % NBY;
		int x = b % NBX;
#endif


/*
		int skip_me;
		if (set_counter > 0)
		{
			if ((StartX <= x) && (x <= EndX) && (StartY <= y) && (y <= EndY) && (StartZ <= z) && (z <= EndZ))
				skip_me = 1;
			else
				skip_me = 0;
		}
		if (skip_me) continue;
*/

                int in_roi = (Xs <= x) && (x <= Xe) && (Ys <= y) && (y <= Ye) && (Zs <= z) && (z <= Ze);
		if (!in_roi)
			continue;

		if ((b < nblocks))
		{
			if (verbose == 2)
				fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z); 

			double zratio = myreader.load_block2(x, y, z, targetdata);
			
			Real bmin = targetdata[0][0][0];
			Real bmax = targetdata[0][0][0];
			
			for (int zb = 0; zb < _BLOCKSIZE_; zb++)
			{
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int xb = 0; xb < _BLOCKSIZE_; xb++)
					{
						Real d = targetdata[zb][yb][xb];
						if (d > bmax) {
							bmax = d;
						}
						if (d < bmin) {
							bmin = d;
						}
					}
			}

			if (verbose == 3)
				fprintf(stdout, "%.3lf bmin = %.3lf bmax = %.3lf\n", zratio, bmin, bmax);
				
			if (bmin < bmax)	// !=
			{
				if (set_counter == 0)
				{
					StartX = EndX = x;
					StartY = EndY = y;
					StartZ = EndZ = z;					
					set_counter++;
				}
				else
				{
					if (x < StartX) StartX = x;
					if (x > EndX)	EndX = x;
					if (y < StartY) StartY = y;
					if (y > EndY)	EndY = y;
					if (z < StartZ) StartZ = z;
					if (z > EndZ)	EndZ = z;
				}
			}

		}
		else {
			//fprintf(stdout, "empty blockk( %d, %d, %d )...\n", x, y, z); 
		}
	}
#endif

	int iStartX, iEndX, iStartY, iEndY, iStartZ, iEndZ;

	MPI_Allreduce(&StartX, &iStartX, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&EndX, &iEndX, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&StartY, &iStartY, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&EndY, &iEndY, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&StartZ, &iStartZ, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&EndZ, &iEndZ, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);


	const double it1 = omp_get_wtime(); 
	if (!mpi_rank)
	{
		fprintf(stdout, "Init ROI = [%d,%d]x[%d,%d]x[%d,%d]\n", iStartX, iEndX, iStartY, iEndY, iStartZ, iEndZ);	
		fprintf(stdout, "Elapsed ROI time = %.3lf seconds\n", it1-it0);
		fflush(0);
	}


	StartX = iStartX;
	EndX = iEndX;
	StartY = iStartY;
	EndY = iEndY;
	StartZ = iStartZ;
	EndZ = iEndZ;
	
	const double t0 = omp_get_wtime(); 

#if 1
	for (int b = mpi_rank; b < b_end; b += mpi_size)	
	{
#if defined(_TRANSPOSE_DATA_)
		int x = b / (NBY * NBZ);
		int y = (b / NBZ) % NBY;
		int z = b % NBZ;
#else
		int z = b / (NBY * NBX);
		int y = (b / NBX) % NBY;
		int x = b % NBX;
#endif

                int in_roi = (Xs <= x) && (x <= Xe) && (Ys <= y) && (y <= Ye) && (Zs <= z) && (z <= Ze);
		if (!in_roi)
			continue;


		int skip_me;
		if (set_counter > 0)
		{
			if ((StartX <= x) && (x <= EndX) && (StartY <= y) && (y <= EndY) && (StartZ <= z) && (z <= EndZ))
				skip_me = 1;
			else
				skip_me = 0;
		}
		if (skip_me) continue;

		if ((b < nblocks))
		{
			if (verbose == 2)
				fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z); 

			double zratio = myreader.load_block2(x, y, z, targetdata);
			
			Real bmin = targetdata[0][0][0];
			Real bmax = targetdata[0][0][0];
			
			for (int zb = 0; zb < _BLOCKSIZE_; zb++)
			{
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int xb = 0; xb < _BLOCKSIZE_; xb++)
					{
						Real d = targetdata[zb][yb][xb];
						if (d > bmax) {
							bmax = d;
						}
						if (d < bmin) {
							bmin = d;
						}
					}
			}

			if (verbose == 3)
				fprintf(stdout, "%.3lf bmin = %.3lf bmax = %.3lf\n", zratio, bmin, bmax);
				
			if (bmin < bmax)	// !=
			{
				if (set_counter == 0)
				{
					StartX = EndX = x;
					StartY = EndY = y;
					StartZ = EndZ = z;					
					set_counter++;
				}
				else
				{
					if (x < StartX) StartX = x;
					if (x > EndX)	EndX = x;
					if (y < StartY) StartY = y;
					if (y > EndY)	EndY = y;
					if (z < StartZ) StartZ = z;
					if (z > EndZ)	EndZ = z;
				}
			}

		}
		else {
			//fprintf(stdout, "empty blockk( %d, %d, %d )...\n", x, y, z); 
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	const double t1 = omp_get_wtime(); 

	if (verbose == 1)
		fprintf(stdout, "Local ROI = [%d,%d]x[%d,%d]x[%d,%d]\n", StartX, EndX, StartY, EndY, StartZ, EndZ);	fflush(0);
	MPI_Barrier(MPI_COMM_WORLD);

	int gStartX, gStartY, gStartZ;
	int gEndX, gEndY, gEndZ;
	

	MPI_Reduce(&StartX, &gStartX, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&EndX, &gEndX, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&StartY, &gStartY, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&EndY, &gEndY, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&StartZ, &gStartZ, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&EndZ, &gEndZ, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (!mpi_rank)
	{
		fprintf(stdout, "Final ROI = [%d,%d]x[%d,%d]x[%d,%d]\n", gStartX, gEndX, gStartY, gEndY, gStartZ, gEndZ);	
		fprintf(stdout, "Init time = %.3lf seconds\n", init_t1-init_t0);
		fprintf(stdout, "Elapsed time = %.3lf seconds\n", t1-t0);
		fflush(0);
	}
#endif
		
	MPI::Finalize();

	return 0;
}
