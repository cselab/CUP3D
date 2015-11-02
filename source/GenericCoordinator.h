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
