//
//  TestBoundaryConditions.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 10/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestBoundaryConditions.h"
#include "Definitions.h"
#include <sstream>
#include <cmath>

void TestBoundaryConditions::_ic()
{
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    const double dh = vInfo[0].h_gridpoint;
    
#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        
        for(int iy=0; iy<FluidBlock::sizeY; iy++)
            for(int ix=0; ix<FluidBlock::sizeX; ix++)
            {
                double p[3];
                info.pos(p, ix, iy);
                p[0] = p[0]*2.-1.;
                p[1] = p[1]*2.-1.;
                
                // this is a test for:
                //	0-dirichlet on x=0
                //	0-neumann on x=1
                
                const int size = 1/dh;
                const int bx = info.index[0]*FluidBlock::sizeX;
                const int by = info.index[1]*FluidBlock::sizeY;
                p[0] = (bx+ix+.5)/(double)size;
                p[1] = (by+iy+.5)/(double)size;
                double x = 4*p[0]*M_PI;
                double y = 3*p[1]*M_PI_2;
                b(ix,iy).divU = cos(y)*sin(x);
            }
    }
    
    
    stringstream ss;
    ss << path2file << "-IC.vti";
    
    dumper.Write(*grid, ss.str());
}

TestBoundaryConditions::TestBoundaryConditions(const int argc, const char ** argv) : Test(argc, argv)
{
    // output settings
    path2file = parser("-file").asString("../data/testBC");
    offset = parser("-offset").asInt(2);
    
    grid = new FluidGrid(2,2,1);
    
    // setup initial condition
    _ic();
}

TestBoundaryConditions::~TestBoundaryConditions()
{
    delete grid;
}

void TestBoundaryConditions::run()
{
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    BlockInfo * ary = &vInfo.front();
    const int N = vInfo.size();
    
    const int stencil_start[3] = {-offset  ,-offset  , 0};
    const int stencil_end[3]   = { offset+1, offset+1, 1};
    
    // 4 channels: N,S,W,E
    
#pragma omp parallel
    {
        Lab lab;
        lab.prepare(*grid, stencil_start, stencil_end, true);
        
#pragma omp for schedule(static)
        for (int i=0; i<N; i++)
        {
            lab.load(ary[i], 0);
            
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            
            const double inv2h = info.h_gridpoint;
            
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    b(ix,iy).u   = lab(ix       ,iy-offset).divU; // N
                    b(ix,iy).v   = lab(ix       ,iy+offset).divU; // S
                    b(ix,iy).p   = lab(ix+offset,iy       ).divU; // W
                    b(ix,iy).chi = lab(ix-offset,iy       ).divU; // E
                }
        }
    }
    
    
    
    stringstream ss;
    ss << path2file << ".vti";
    
    dumper.Write(*grid, ss.str());
}

void TestBoundaryConditions::check()
{
    // this can be checked - the operator is just a translation operator
}