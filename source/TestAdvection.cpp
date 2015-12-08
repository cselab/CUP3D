//
//  TestAdvection.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestAdvection.h"
#include "ProcessOperatorsOMP.h"
#include <sstream>
#include <cmath>

#include "CoordinatorAdvection.h"

/*
 class BS4
 {
 public:
	static inline Real eval(Real x)
	{
 const Real t = fabs(x);
 
 if (t>2) return 0;
 
 if (t>1) return pow(2-t,3)/6;
 
 return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
	}
 };
 */

void TestAdvection::_icLinear()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy, iz);
				
                b(ix, iy, iz).rho = sin(p[0]*8.*M_PI);//*sin(p[1]*2.*M_PI);
				b(ix, iy, iz).u   = 1;
				b(ix, iy, iz).v   = 1;
				b(ix, iy, iz).w   = 1;
				b(ix, iy, iz).chi = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti" ;
	//cout << ss.str() << endl;
	
	dumper.Write(*grid, ss.str());
}

void TestAdvection::_icVortex()
{
	const double center[3] = {.5,.5,.5};
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy, iz);
				
				p[0] = p[0]*2.-1.;
				p[1] = p[1]*2.-1.;
				p[2] = p[2]*2.-1.;
				
                const Real r = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                const Real invR = 1./r;
								
                b(ix, iy, iz).rho = r;
				b(ix, iy, iz).u   =   sin(p[1])*cos(r*M_PI/2)*invR;//-p[1];//
				b(ix, iy, iz).v   =  -sin(p[0])*cos(r*M_PI/2)*invR;// p[0];//
				b(ix, iy, iz).w   =  0;
				b(ix, iy, iz).chi = 0;
				/*
				if (r>.5)
				{
					b(ix,iy).rho = 1;
				}
				/*
				 p[0] = p[0]*2.-1.;
				 p[1] = p[1]*2.-1.;
				 
				 const Real IrI  = sqrt(p[0]*p[0] + p[1]*p[1]);
				 const double strength = 100./(1+IrI*IrI)*BS4::eval(IrI/0.5*2.5);
				 
				 b(ix, iy).rho = p[0];
				 b(ix, iy).u   = -p[1]*strength;
				 b(ix, iy).v   =  p[0]*strength;
				 b(ix, iy).chi = 0;
				 
				/*/
				/*
				const Real dx = p[0] - center[0];
				const Real dy = p[1] - center[1];
				const Real dist = sqrt(dx*dx + dy*dy);
				
				b(ix, iy).rho = abs(p[0]-.5);
				if (dist <= .4)
				{
					const Real amplitude = .5*(cos(dist*5*M_PI+M_PI)+1);
					b(ix, iy).u = -amplitude*dy/dist;
					b(ix, iy).v =  amplitude*dx/dist;
				}
				else
				{
					b(ix, iy).u = 0;
					b(ix, iy).v = 0;
				}
				b(ix, iy).chi = 0;
				
				b(ix, iy).tmpU = 0;
				b(ix, iy).tmpV = 0;
				b(ix, iy).tmp  = 0;
				
				//*/
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	
	dumper.Write(*grid, ss.str());
}

void TestAdvection::_icBurger()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				b(ix, iy, iz).rho = 1;
				b(ix, iy, iz).u   = ix<FluidBlock::sizeX/2 ? p[0] : (1-p[0]);//1-p[0];//1-cos(p[0]*M_PI*2);
				b(ix, iy, iz).v   = 0;
				b(ix, iy, iz).w   = 0;
				b(ix, iy, iz).chi = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti" ;
	//cout << ss.str() << endl;
	
	dumper.Write(*grid, ss.str());
}

TestAdvection::TestAdvection(const int argc, const char ** argv, int testCase, const int bpd, const double dt, const int nsteps) : Test(argc, argv), time(0), testCase(testCase), bpd(bpd), dt(dt), nsteps(nsteps)
{
	grid = new FluidGrid(bpd,bpd,bpd);
	
	// setup initial condition
	if (testCase==0)
	{
		// output settings
		path2file = parser("-file").asString("../data/testAdvectionLinear");
		_icLinear();
	}
	else if (testCase==1)
	{
		// output settings
		path2file = parser("-file").asString("../data/testAdvectionVortex");
		_icVortex();
		//path2file = parser("-file").asString("../data/testAdvectionBurger");
		//_icBurger();
	}
	else
	{
		cout << "unknown test case - aborting!\n";
		abort();
	}
}

TestAdvection::~TestAdvection()
{
	delete grid;
}

void TestAdvection::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	
	int step = 0;
    
    if (nsteps==1)
    {
        CoordinatorTransport<Lab> coordTransport(grid);
        coordTransport(dt);
        time += dt;
        step++;
    }
    else
    {
        CoordinatorTransportTimeTest<Lab> coordTransport(grid);
        while(step<nsteps)
        {
            coordTransport(dt);
            
            time += dt;
            step++;
        }
    }
	
	stringstream ss;
	ss << path2file << "-test" << testCase << "-bpd" << bpd << ".vti";
	dumper.Write(*grid, ss.str());
}

void TestAdvection::check()
{
	const double center[3] = {.5,.5,.5};
	
	//cout << "\tErrors (uLinf, uL1, uL2):\t";
	double uLinf = 0.;
	double uL1 = 0.;
	double uL2 = 0.;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
	
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	const int sizeZ = bpd * FluidBlock::sizeZ;
	
	if (testCase==0)
	{
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy, iz);
					
					double error;
                    error = b(ix, iy).rho - sin((p[0]-time)*8.*M_PI);//*sin((p[1]+dt)*2.*M_PI);
                    b(ix,iy).chi = error;
                    
					
					uLinf = max(uLinf,abs(error));
					uL1 += abs(error);
					uL2 += error*error;
				}
		}
	}
	else
	{
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy, iz);
					
					p[0] = p[0]*2.-1.;
					p[1] = p[1]*2.-1.;
					p[2] = p[2]*2.-1.;
					
					Real r = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
					
                    double error = b(ix, iy, iz).rho - r;
					b(ix,iy).chi = error;
					
					uLinf = max(uLinf,abs(error));
					uL1 += abs(error);
					uL2 += error*error;
				}
		}
		/*
		 #pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		 for (int iy=0; iy<sizeY; iy++)
			for (int ix=0; ix<sizeX; ix++)
			{
		 double error = vorticity(ix,iy)-(*vorticityIC)(ix,iy);
		 vorticityDiff(ix,iy) = error;
		 
		 uLinf = max(uLinf,abs(error));
		 uL1 += abs(error);
		 uL2 += error*error;
			}
		 */
		/*
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					// 1D Burger's
					double p[3];
					info.pos(p, ix, iy);
					double pIC = p[0] - b(ix,iy).u * time; // why this factor 2?
					double error = b(ix, iy).u - (1-cos(pIC*M_PI*2));
					
					b(ix,iy).u = 1-cos(pIC*M_PI*2);
					
					uLinf = max(uLinf,abs(error));
					uL1 += abs(error);
					uL2 += error*error;
				}
		}
		 */
	}
	
	//stringstream sVort;
	//sVort << path2file << "VorticityDiff-" << bpd << ".vti";
	
	stringstream ssol;
	ssol << path2file << "-solution" << testCase << "-bpd" << bpd << ".vti";
	dumper.Write(*grid, ssol.str());
	
	uL1 *= dh*dh*dh;
	uL2 = sqrt(uL2*dh*dh*dh);
	cout << uLinf << "\t" << uL1 << "\t" << uL2 << endl;
	myfile << sizeX << " " << uLinf << " " << uL1 << " " << uL2 << endl;
}
