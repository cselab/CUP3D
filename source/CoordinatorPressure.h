//
//  CoordinatorPressure.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "IF3D_ObstacleVector.h"
#include "PoissonSolverScalarFFTW_MPI.h"

struct PressureObstacleVisitor : public ObstacleVisitor
{
	FluidGridMPI * grid;
    vector<BlockInfo> vInfo;

    PressureObstacleVisitor(FluidGridMPI * grid) : grid(grid)
    {
      vInfo = grid->getBlocksInfo();
    }

     void visit(IF3D_ObstacleOperator* const obstacle)
     {
#pragma omp parallel
         {
             const std::map<int,ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
#pragma omp for schedule(static)
             for(int i=0; i<vInfo.size(); i++) {
            	 BlockInfo info = vInfo[i];
                 const auto pos = obstblocks.find(info.blockID);
                 if(pos == obstblocks.end()) continue;
                 
                 FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;

                 for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
                 for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				 for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					 // what if multiple obstacles share a block??
					 // let's plus equal and wake up during the night to stress about it
					 b(ix,iy,iz).tmpU += pos->second->udef[iz][iy][ix][0];
					 b(ix,iy,iz).tmpV += pos->second->udef[iz][iy][ix][1];
					 b(ix,iy,iz).tmpW += pos->second->udef[iz][iy][ix][2];
				 }
             }
         }
     }
};

struct PressurePenaltyVisitor : ObstacleVisitor
{
	FluidGridMPI * grid;
    const Real dt, lambda;
    Real * const uInf;
    //Real ext_X, ext_Y, ext_Z;
    vector<BlockInfo> vInfo;

    PressurePenaltyVisitor(FluidGridMPI * grid, const Real dt, const Real lambda, Real* const uInf)
    : grid(grid), dt(dt), lambda(lambda), uInf(uInf)
    {
        vInfo = grid->getBlocksInfo();
        /*
        const Real extent = grid->maxextent;
        const unsigned int maxbpd = max(grid->NX*FluidBlock::sizeX,
        							max(grid->NY*FluidBlock::sizeY,
        								grid->NZ*FluidBlock::sizeZ));
        const Real scale[3] = {
        		(Real)(grid->NX*FluidBlock::sizeX)/(Real)maxbpd,
        		(Real)(grid->NY*FluidBlock::sizeY)/(Real)maxbpd,
        		(Real)(grid->NZ*FluidBlock::sizeZ)/(Real)maxbpd
        };
        ext_X = scale[0]*extent;
        ext_Y = scale[1]*extent;
        ext_Z = scale[2]*extent;
        */
    }

     void visit(IF3D_ObstacleOperator* const obstacle)
     {

#pragma omp parallel
         {
            const std::map<int, ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
            Real uBody[3], omegaBody[3], centerOfMass[3];
            obstacle->getCenterOfMass(centerOfMass);
            obstacle->getTranslationVelocity(uBody);
            obstacle->getAngularVelocity(omegaBody);
            /*
			#pragma omp master
            printf("%f %f %f %f %f %f %f %f %f %f %f %f\n",
            		uBody[0],uBody[1],uBody[2],omegaBody[0],omegaBody[1],omegaBody[2],
            		centerOfMass[0],centerOfMass[1],centerOfMass[2],uInf[0],uInf[1],uInf[2]);
			*/
#pragma omp for schedule(static)
            for(int i=0; i<vInfo.size(); i++) {
            	BlockInfo info = vInfo[i];
            	const auto pos = obstblocks.find(info.blockID);
            	if(pos == obstblocks.end()) continue;
            	FluidBlock& b = *(FluidBlock*)info.ptrBlock;

            	for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				if (pos->second->chi[iz][iy][ix] > 0) {
					Real p[3];
					info.pos(p, ix, iy, iz);
					p[0]-=centerOfMass[0];
					p[1]-=centerOfMass[1];
					p[2]-=centerOfMass[2];
					const Real lambdaChi  = lambda * pos->second->chi[iz][iy][ix];
					const Real object_UR[3] = {
							omegaBody[1]*p[2]-omegaBody[2]*p[1],
							omegaBody[2]*p[0]-omegaBody[0]*p[2],
							omegaBody[0]*p[1]-omegaBody[1]*p[0]
					};
					const Real object_UDEF[3] = {
							pos->second->udef[iz][iy][ix][0],
							pos->second->udef[iz][iy][ix][1],
							pos->second->udef[iz][iy][ix][2]
					};
					const Real U_TOT[3] = {
							uBody[0]+object_UR[0]+object_UDEF[0]-uInf[0],
							uBody[1]+object_UR[1]+object_UDEF[1]-uInf[1],
							uBody[2]+object_UR[2]+object_UDEF[2]-uInf[2]
					};
    				b(ix,iy,iz).tmpU += lambdaChi*(U_TOT[0]-b(ix,iy,iz).u);
    				b(ix,iy,iz).tmpV += lambdaChi*(U_TOT[1]-b(ix,iy,iz).v);
    				b(ix,iy,iz).tmpW += lambdaChi*(U_TOT[2]-b(ix,iy,iz).v);
				}
             }
         }
     }
};

class PressRHSOperator : public GenericLabOperator
{
private:
    double dt;

public:
    PressRHSOperator(double dt) : dt(dt)
    {
    	stencil = StencilInfo(-1,-1,-1, 2,2,2,false,6,0,1,2,5,6,7);
        stencil_start[0] = -1;
        stencil_start[1] = -1;
        stencil_start[2] = -1;
        stencil_end[0] = 2;
        stencil_end[1] = 2;
        stencil_end[2] = 2;
    }

    template <typename Lab, typename BlockType>
    void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
    {
    	const Real fac1 = 0.5/(info.h_gridpoint);
    	const Real fac2 = .25/(info.h_gridpoint*info.h_gridpoint);
    	for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    	for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			const FluidElement& phiW = lab(ix-1,iy  ,iz  );
			const FluidElement& phiE = lab(ix+1,iy  ,iz  );
			const FluidElement& phiS = lab(ix  ,iy-1,iz  );
			const FluidElement& phiN = lab(ix  ,iy+1,iz  );
			const FluidElement& phiF = lab(ix  ,iy  ,iz-1);
			const FluidElement& phiB = lab(ix  ,iy  ,iz+1);
			const Real dudx = (phiE.u - phiW.u);
			const Real dudy = (phiN.u - phiS.u);
			const Real dudz = (phiB.u - phiF.u);
			const Real dvdx = (phiE.v - phiW.v);
			const Real dvdy = (phiN.v - phiS.v);
			const Real dvdz = (phiB.v - phiF.v);
			const Real dwdx = (phiE.w - phiW.w);
			const Real dwdy = (phiN.w - phiS.w);
			const Real dwdz = (phiB.w - phiF.w);
			const Real divPen = phiE.tmpU-phiW.tmpU+phiN.tmpV-phiS.tmpV+phiB.tmpW-phiF.tmpW;
			const Real contr = dudx*dudx+dvdy*dvdy+dwdz*dwdz+2*(dudy*dvdx+dudz*dwdx+dvdz*dwdy);
			o(ix,iy,iz).p = fac1*divPen - fac2*contr;
		}
    }
};

class OperatorDivergenceMinusDivTmpU : public GenericLabOperator
{
private:
	Real dt;
	
public:
	OperatorDivergenceMinusDivTmpU(Real dt) : dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 6, 0,1,2,5,6,7);
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	~OperatorDivergenceMinusDivTmpU() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = 0.5/(info.h_gridpoint * dt);
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			// Poisson solver reads field p for the rhs
			const Real uW    = lab(ix-1,iy  ,iz  ).u;
			const Real uE    = lab(ix+1,iy  ,iz  ).u;
			const Real vS    = lab(ix  ,iy-1,iz  ).v;
			const Real vN    = lab(ix  ,iy+1,iz  ).v;
			const Real wF    = lab(ix  ,iy  ,iz-1).w;
			const Real wB    = lab(ix  ,iy  ,iz+1).w;
			const Real uWdef = lab(ix-1,iy  ,iz  ).tmpU;
			const Real uEdef = lab(ix+1,iy  ,iz  ).tmpU;
			const Real vSdef = lab(ix  ,iy-1,iz  ).tmpV;
			const Real vNdef = lab(ix  ,iy+1,iz  ).tmpV;
			const Real wFdef = lab(ix  ,iy  ,iz-1).tmpW;
			const Real wBdef = lab(ix  ,iy  ,iz+1).tmpW;
			o(ix,iy,iz).p = factor * (uE - uW + vN - vS + wB - wF
					-o(ix,iy,iz).chi*(uEdef-uWdef+vNdef-vSdef+wBdef-wFdef));
			o(ix,iy,iz).chi = o(ix,iy,iz).p;
		}
	}
};

class OperatorDivergenceMinusDivTmpU2ndOrder : public GenericLabOperator
{
private:
	Real dt;

public:
	OperatorDivergenceMinusDivTmpU(Real dt) : dt(dt)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 6, 0,1,2,5,6,7);
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	~OperatorDivergenceMinusDivTmpU() {}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = 1./(12.*info.h_gridpoint * dt);

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			// Poisson solver reads field p for the rhs
			const Real uW1    = lab(ix-1,iy  ,iz  ).u;
			const Real uE1    = lab(ix+1,iy  ,iz  ).u;
			const Real vS1    = lab(ix  ,iy-1,iz  ).v;
			const Real vN1    = lab(ix  ,iy+1,iz  ).v;
			const Real wF1    = lab(ix  ,iy  ,iz-1).w;
			const Real wB1    = lab(ix  ,iy  ,iz+1).w;
			const Real uW2    = lab(ix-2,iy  ,iz  ).u;
			const Real uE2    = lab(ix+2,iy  ,iz  ).u;
			const Real vS2    = lab(ix  ,iy-2,iz  ).v;
			const Real vN2    = lab(ix  ,iy+2,iz  ).v;
			const Real wF2    = lab(ix  ,iy  ,iz-2).w;
			const Real wB2    = lab(ix  ,iy  ,iz+2).w;
			const Real uWd1 = lab(ix-1,iy  ,iz  ).tmpU;
			const Real uEd1 = lab(ix+1,iy  ,iz  ).tmpU;
			const Real vSd1 = lab(ix  ,iy-1,iz  ).tmpV;
			const Real vNd1 = lab(ix  ,iy+1,iz  ).tmpV;
			const Real wFd1 = lab(ix  ,iy  ,iz-1).tmpW;
			const Real wBd1 = lab(ix  ,iy  ,iz+1).tmpW;
			const Real uWd2 = lab(ix-2,iy  ,iz  ).tmpU;
			const Real uEd2 = lab(ix+2,iy  ,iz  ).tmpU;
			const Real vSd2 = lab(ix  ,iy-2,iz  ).tmpV;
			const Real vNd2 = lab(ix  ,iy+2,iz  ).tmpV;
			const Real wFd2 = lab(ix  ,iy  ,iz-2).tmpW;
			const Real wBd2 = lab(ix  ,iy  ,iz+2).tmpW;
			const Real tmp1 = 8*(uE1 - uW1 + vN1 - vS1 + wB1 - wF1);
			const Real tmp2 =  -(uE2 - uW2 + vN2 - vS2 + wB2 - wF2);
			const Real tmp3 = 8*(uEd1 - uWd1 + vNd1 - vSd1 + wBd1 - wFd1);
			const Real tmp4 =  -(uEd2 - uWd2 + vNd2 - vSd2 + wBd2 - wFd2);
			o(ix,iy,iz).p = factor*(tmp1-tmp2 -o(ix,iy,iz).chi*(tmp3-tmp4));
			//o(ix,iy,iz).chi = o(ix,iy,iz).p;
		}
	}
};

class OperatorGradP : public GenericLabOperator
{
private:
	Real dt;
	
public:
	OperatorGradP(Real dt) : dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);
		
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorGradP() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real prefactor = -.5 * dt / (info.h_gridpoint);

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			// p contains the pressure correction after the Poisson solver
			const Real divUW = lab(ix-1,iy  ,iz  ).p;
			const Real divUE = lab(ix+1,iy  ,iz  ).p;
			const Real divUS = lab(ix  ,iy-1,iz  ).p;
			const Real divUN = lab(ix  ,iy+1,iz  ).p;
			const Real divUF = lab(ix  ,iy  ,iz-1).p;
			const Real divUB = lab(ix  ,iy  ,iz+1).p;

			o(ix,iy,iz).u += prefactor * (divUE - divUW);//
			o(ix,iy,iz).v += prefactor * (divUN - divUS);//
			o(ix,iy,iz).w += prefactor * (divUB - divUF);//

			assert(!std::isnan(o(ix,iy,iz).u));
			assert(!std::isnan(o(ix,iy,iz).v));
			assert(!std::isnan(o(ix,iy,iz).w));
		}
	}
};

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
protected:
    IF3D_ObstacleVector** const obstacleVector;
//#ifndef _MIXED_
	PoissonSolverScalarFFTW_MPI<FluidGridMPI, StreamerDiv> pressureSolver;
//#else
//	PoissonSolverScalarFFTW_MPI_DCT<FluidGridMPI, StreamerDiv> pressureSolver;
//#endif // _MIXED_

public:
	CoordinatorPressure(FluidGridMPI * grid, IF3D_ObstacleVector** const myobstacles) :
		GenericCoordinator(grid), pressureSolver(NTHREADS,*grid), obstacleVector(myobstacles)
	{
	}
	
	void operator()(const Real dt)
	{
#pragma omp parallel
		{
			const int N = vInfo.size();
#pragma omp for schedule(static)
			for(int i=0; i<vInfo.size(); i++) {
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					b(ix,iy,iz).tmpU = 0;
					b(ix,iy,iz).tmpV = 0;
					b(ix,iy,iz).tmpW = 0; //zero fields, going to contain Udef
				}
			}
		}

	   //store deformation velocities onto tmp fields
		ObstacleVisitor * pressureVisitor = new PressureObstacleVisitor(grid);
		(*obstacleVector)->Accept(pressureVisitor); // accept you son of a french cow
		delete pressureVisitor;
		
		{ 	//place onto p: ( div u^(t+1) - div u^* ) / dt
			//where i want div u^(t+1) to be equal to div udef
			OperatorDivergenceMinusDivTmpU2ndOrder kernelDiv(dt);
			compute(kernelDiv);
		}
		
			pressureSolver.solve(*grid);
		
		{ //pressure correction dudt* = - grad P / rho
			OperatorGradP kernelGradP(dt);
			compute(kernelGradP);
		}
	}
	
	string getName()
	{
		return "Pressure";
	}
};
#endif
