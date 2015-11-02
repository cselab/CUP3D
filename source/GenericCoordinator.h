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
	FluidGrid * grid;
	vector<BlockInfo> vInfo;
	
	inline void check(string infoText)
	{
#ifndef NDEBUG
		int rank;
#ifdef _MULTIGRID_
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif // _MULTIGRID_
		
		if (rank==0)
		{
			const int N = vInfo.size();
			
#pragma omp parallel for schedule(static)
			for(int i=0; i<N; i++)
			{
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						if (std::isnan(b(ix,iy).rho) ||
							std::isnan(b(ix,iy).u) ||
							std::isnan(b(ix,iy).v) ||
							std::isnan(b(ix,iy).chi) ||
							std::isnan(b(ix,iy).p) ||
							std::isnan(b(ix,iy).pOld))
							cout << infoText.c_str() << endl;
						
						if (b(ix,iy).rho <= 0)
							cout << infoText.c_str() << endl;
						
						assert(b(ix,iy).rho > 0);
						assert(!std::isnan(b(ix,iy).rho));
						assert(!std::isnan(b(ix,iy).u));
						assert(!std::isnan(b(ix,iy).v));
						assert(!std::isnan(b(ix,iy).chi));
						assert(!std::isnan(b(ix,iy).p));
						assert(!std::isnan(b(ix,iy).pOld));
						assert(!std::isnan(b(ix,iy).tmpU));
						assert(!std::isnan(b(ix,iy).tmpV));
						assert(!std::isnan(b(ix,iy).tmp));
						assert(!std::isnan(b(ix,iy).divU));
						assert(b(ix,iy).rho < 1e10);
						assert(b(ix,iy).u < 1e10);
						assert(b(ix,iy).v < 1e10);
						assert(b(ix,iy).p < 1e10);
					}
			}
		}
#endif
	}
	
public:
	GenericCoordinator(FluidGrid * grid) : grid(grid)
	{
		vInfo = grid->getBlocksInfo();
	}
	
	virtual void operator()(const double dt) = 0;
	
	virtual string getName() = 0;
};

#endif
