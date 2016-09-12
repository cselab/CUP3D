//
//  IF2D_ObstacleFactory.h
//  IF2D_ROCKS
//
//  Created by Wim van Rees on 06/10/14.
//
//

#ifndef __IF2D_ROCKS__IF2D_ObstacleFactory__
#define __IF2D_ROCKS__IF2D_ObstacleFactory__

#include "IF2D_ObstacleOperator.h"

class IF2D_ObstacleFactory
{
    const Real Uinf[3];
    FluidGrid * grid;
    int _getlines(std::string filename);
    
public:
    IF2D_ObstacleFactory(FluidGrid * grid, const Real* Uinf)
	: grid(grid), Uinf{Uinf[0],Uinf[1],Uinf[2]}
    {}
    
    ~IF2D_ObstacleFactory()
    {}
    
    std::vector<IF2D_ObstacleOperator * > create(ArgumentParser & parser);
};


#endif /* defined(__IF2D_ROCKS__IF2D_ObstacleFactory__) */
