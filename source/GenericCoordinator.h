//
//  GenericCoordinator.h
//  CubismUP_3D
//
//	This class serves as the interface for a coordinator object
//	A coordinator object schedules the processing of blocks with its operator
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_GenericCoordinator_h
#define CubismUP_3D_GenericCoordinator_h

#include "Definitions.h"

class GenericCoordinator
{
protected:
	FluidGridMPI * grid;
	vector<BlockInfo> vInfo;
	
	inline void check(string infoText)
	{
#ifndef NDEBUG
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						if (std::isnan(b(ix,iy,iz).rho) ||
							std::isnan(b(ix,iy,iz).u) ||
							std::isnan(b(ix,iy,iz).v) ||
							std::isnan(b(ix,iy,iz).w) ||
							std::isnan(b(ix,iy,iz).chi) ||
							std::isnan(b(ix,iy,iz).p) ||
							std::isnan(b(ix,iy,iz).pOld))
							cout << infoText.c_str() << endl;
						
						if (b(ix,iy,iz).rho <= 0)
							cout << infoText.c_str() << endl;
						
						assert(b(ix,iy,iz).rho > 0);
						assert(!std::isnan(b(ix,iy,iz).rho));
						assert(!std::isnan(b(ix,iy,iz).u));
						assert(!std::isnan(b(ix,iy,iz).v));
						assert(!std::isnan(b(ix,iy,iz).w));
						assert(!std::isnan(b(ix,iy,iz).chi));
						assert(!std::isnan(b(ix,iy,iz).p));
						assert(!std::isnan(b(ix,iy,iz).pOld));
						assert(!std::isnan(b(ix,iy,iz).tmpU));
						assert(!std::isnan(b(ix,iy,iz).tmpV));
						assert(!std::isnan(b(ix,iy,iz).tmpW));
						assert(!std::isnan(b(ix,iy,iz).tmp));
						assert(!std::isnan(b(ix,iy,iz).divU));
						assert(b(ix,iy,iz).rho < 1e10);
						assert(b(ix,iy,iz).u < 1e10);
						assert(b(ix,iy,iz).v < 1e10);
						assert(b(ix,iy,iz).w < 1e10);
						assert(b(ix,iy,iz).p < 1e10);
					}
		}
#endif
	}
	
	template <typename Kernel>
	void compute(Kernel kernel)
	{
		SynchronizerMPI& Synch = grid->sync(kernel);
		
		vector<BlockInfo> avail0, avail1;
		
		const int nthreads = omp_get_max_threads();
		
		LabMPI * labs = new LabMPI[nthreads];
		
		for(int i = 0; i < nthreads; ++i)
			labs[i].prepare(*grid, Synch);
		
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
			delete [] labs;
			labs=NULL;
		}
		
		MPI::COMM_WORLD.Barrier();
	}
	
public:
	GenericCoordinator(FluidGridMPI * grid) : grid(grid)
	{
		vInfo = grid->getBlocksInfo();
	}
	
	virtual void operator()(const double dt) = 0;
	
	virtual string getName() = 0;
};

#endif
