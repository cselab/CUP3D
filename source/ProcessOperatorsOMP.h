//
//  ProcessOperatorsOMP.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_ProcessOperators_h
#define CubismUP_3D_ProcessOperators_h

#include "Definitions.h"

// -gradp, divergence, advection
template<typename Lab, typename Kernel>
void processOMP(Real dt, vector<BlockInfo>& vInfo, FluidGridMPI & grid)
{
	const Kernel kernel(dt);
 #if 0
	SynchronizerMPI& Synch = grid.sync(kernel);

	vector<BlockInfo> avail0, avail1;

	const int nthreads = omp_get_max_threads();

	LabMPI * labs = new LabMPI[nthreads];

	for(int i = 0; i < nthreads; ++i)
		labs[i].prepare(grid, Synch);

	static int rounds = -1;
	static int one_less = 1;
	if (rounds == -1)
	{
		char *s = getenv("MYROUNDS");
		if (s != NULL)
			rounds = atoi(s);
		else
			rounds = 0;

		char *s2 = getenv("USEMAXTHREADS");
		if (s2 != NULL)
			one_less = !atoi(s2);
	}

	MPI_Barrier(grid.getCartComm());

	avail0 = Synch.avail_inner();
	const int Ninner = avail0.size();
	BlockInfo * ary0 = &avail0.front();

	int nthreads_first;
	if (one_less)
		nthreads_first = nthreads-1;
	else
		nthreads_first = nthreads;

	if (nthreads_first == 0) nthreads_first = 1;

	int Ninner_first = (nthreads_first)*rounds;
	if (Ninner_first > Ninner) Ninner_first = Ninner;
	int Ninner_rest = Ninner - Ninner_first;

	#pragma omp parallel num_threads(nthreads_first)
	{
		int tid = omp_get_thread_num();
		LabMPI& lab = labs[tid];

		#pragma omp for schedule(dynamic,1)
		for(int i=0; i<Ninner_first; i++)
		{
			lab.load(ary0[i], 0);
			kernel(lab, ary0[i], *(FluidBlock*)ary0[i].ptrBlock); // why is this using the local blockInfo? or is it global? is dh correct?
		}
	}

	avail1 = Synch.avail_halo();
	const int Nhalo = avail1.size();
	BlockInfo * ary1 = &avail1.front();

	#pragma omp parallel num_threads(nthreads)
	{
		int tid = omp_get_thread_num();
		LabMPI& lab = labs[tid];

		#pragma omp for schedule(dynamic,1)
		for(int i=-Ninner_rest; i<Nhalo; i++)
		{
			if (i < 0)
			{
				int ii = i + Ninner;
				lab.load(ary0[ii], 0);
				kernel(lab, ary0[ii], *(FluidBlock*)ary0[ii].ptrBlock);
			}
			else
			{
				lab.load(ary1[i], 0);
				kernel(lab, ary1[i], *(FluidBlock*)ary1[i].ptrBlock);
			}
		}
	}

	if(labs!=NULL)
	{
		delete[] labs;
		labs=NULL;
	}

	MPI_Barrier(grid.getCartComm());
 #else
		SynchronizerMPI& Synch = grid.sync(kernel);

		const int nthreads = omp_get_max_threads();
		LabMPI * labs = new LabMPI[nthreads];
		for(int i = 0; i < nthreads; ++i)
			labs[i].prepare(grid, Synch);

		MPI_Barrier(grid.getCartComm());
		vector<BlockInfo> avail0 = Synch.avail_inner();
		const int Ninner = avail0.size();

		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			LabMPI& lab = labs[tid];

			#pragma omp for schedule(dynamic,1)
			for(int i=0; i<Ninner; i++) {
				BlockInfo info = avail0[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				lab.load(info, 0);
				kernel(lab, info, b); // why is this using the local blockInfo? or is it global? is dh correct?
			}
		}

		vector<BlockInfo> avail1 = Synch.avail_halo();
		const int Nhalo = avail1.size();

		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			LabMPI& lab = labs[tid];

			#pragma omp for schedule(dynamic,1)
			for(int i=0; i<Nhalo; i++) {
					BlockInfo info = avail1[i];
					FluidBlock& b = *(FluidBlock*)info.ptrBlock;
					lab.load(info, 0);
					kernel(lab, info, b);
			}
		}

		if(labs!=NULL) {
			delete [] labs;
			labs=NULL;
		}

		MPI_Barrier(grid.getCartComm());
 #endif
}

template<typename Lab, typename Kernel>
void processOMPold(Real dt, vector<BlockInfo>& vInfo, FluidGridMPI & grid)
{
	const int N = vInfo.size();

	#pragma omp parallel
	{
		Lab lab;
		Kernel kernel(dt);
		lab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);

		#pragma omp for schedule(static)
		for (int i=0; i<N; i++) {
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			lab.load(info, 0);
			kernel(lab, info, b);
		}
	}
}

static Real findMaxUOMP(vector<BlockInfo>& myInfo, FluidGridMPI& grid, const Real*const uInf)
{
	Real maxU = 0;
	const int N = myInfo.size();

	#pragma omp parallel for schedule(static) reduction(max:maxU)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			maxU = max(maxU,(Real)abs(b(ix,iy,iz).u +uInf[0]));
			maxU = max(maxU,(Real)abs(b(ix,iy,iz).v +uInf[1]));
			maxU = max(maxU,(Real)abs(b(ix,iy,iz).w +uInf[2]));
		}
	}

	return maxU;
};

static void copyDumpGrid(FluidGridMPI& grid, DumpGridMPI& dump)
{
  vector<BlockInfo> vInfo1 = grid.getBlocksInfo();
  vector<BlockInfo> vInfo2 = dump.getBlocksInfo();
	const int N = vInfo1.size();
	if(vInfo1.size() != vInfo2.size()) {
		 printf("Async dump fail 1.\n");
		 fflush(0);
		 MPI_Abort(grid.getCartComm(), MPI_ERR_OTHER);
 	}
	#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++) {
		BlockInfo info1 = vInfo1[i];
		BlockInfo info2 = vInfo2[i];
		Real p1[3], p2[3];
		info1.pos(p1, 0,0,0);
		info2.pos(p2, 0,0,0);
		if (fabs(p1[0]-p2[0])>info1.h_gridpoint ||
				fabs(p1[1]-p2[1])>info1.h_gridpoint ||
				fabs(p1[2]-p2[2])>info1.h_gridpoint) {
					 printf("Async dump fail 2.\n");
					 fflush(0);
					 MPI_Abort(grid.getCartComm(), MPI_ERR_OTHER);
				}
		FluidBlock& b = *(FluidBlock*)info1.ptrBlock;
		 DumpBlock& d = *( DumpBlock*)info2.ptrBlock;
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			d(ix,iy,iz).u = b(ix,iy,iz).u;
			d(ix,iy,iz).v = b(ix,iy,iz).v;
			d(ix,iy,iz).w = b(ix,iy,iz).w;
			d(ix,iy,iz).chi = b(ix,iy,iz).chi;
		}
	}
};
#endif
