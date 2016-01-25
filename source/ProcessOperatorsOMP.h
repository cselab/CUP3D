//
//  ProcessOperatorsOMP.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_ProcessOperators_h
#define CubismUP_3D_ProcessOperators_h

#include "Shape.h"
#include "Definitions.h"

// -gradp, divergence, advection
template<typename Lab, typename Kernel>
void processOMP(double dt, vector<BlockInfo>& vInfo, FluidGridMPI & grid)
{
	Kernel kernel(dt);
	
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
	
	MPI::COMM_WORLD.Barrier();
	
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
	
	MPI::COMM_WORLD.Barrier();
}

template<typename Lab, typename Kernel>
void processOMPold(double dt, vector<BlockInfo>& vInfo, FluidGridMPI & grid)
{
	BlockInfo * ary = &vInfo.front();
	const int N = vInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(dt);
		
		Lab lab;
		lab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);
		
#pragma omp for schedule(static)
		for (int i=0; i<N; i++)
		{
			lab.load(ary[i], 0);
			kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

double findMaxUOMP(vector<BlockInfo>& myInfo, FluidGridMPI & grid);
#endif
