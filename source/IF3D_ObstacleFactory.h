//
//  IF3D_ObstacleFactory.h
//  IF3D_ROCKS
//
//  Created by Wim van Rees on 06/10/14.
//
//

#ifndef __IF3D_ROCKS__IF3D_ObstacleFactory__
#define __IF3D_ROCKS__IF3D_ObstacleFactory__

#include "IF3D_ObstacleOperator.h"

class IF3D_ObstacleFactory
{
    const Real Uinf[3];
    int rank;
    FluidGridMPI * grid;
    int _getlines(std::string filename);
    
public:
    IF3D_ObstacleFactory(FluidGridMPI * grid, const Real* Uinf)
	: grid(grid), Uinf{Uinf[0],Uinf[1],Uinf[2]}
    {
    	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }
    
    ~IF3D_ObstacleFactory()
    {}
    
    std::vector<IF3D_ObstacleOperator * > create(ArgumentParser & parser);
};


#endif /* defined(__IF3D_ROCKS__IF3D_ObstacleFactory__) */
