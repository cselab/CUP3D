//
//  TestBoundaryConditions.h
//  CubismUP_3D
//
//  Created by Christian Conti on 10/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestBoundaryConditions__
#define __CubismUP_3D__TestBoundaryConditions__

#include <stdio.h>
#include "Test.h"

class TestBoundaryConditions : public Test
{
private:
    int offset;
    
    string path2file;
    
    void _ic();
    
public:
    TestBoundaryConditions(const int argc, const char ** argv);
    ~TestBoundaryConditions();
    
    void run();
    void check();
};

#endif /* defined(__CubismUP_3D__TestBoundaryConditions__) */
