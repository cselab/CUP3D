//
//  IF3D_FishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_VortexOperator.h"
#include "IF3D_VortexLibrary.h"
#include <chrono>
IF3D_VortexOperator::IF3D_VortexOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_ObstacleOperator(grid, parser), created(false)
{
	volume=0;
	for(int i=0;i<3;i++) transVel[i]=0;
	for(int i=0;i<3;i++) angVel[i]=0;
	for(int i=0;i<6;i++) J[i]=0;
}

void IF3D_VortexOperator::create(const int step_id,const Real time, const Real dt, const Real *Uinf)
{
	if(step_id || created) return;

	#pragma omp parallel
	{
		const int N = vInfo.size();
		const Real alpha = 1.25643, fac = 1.39795; //for some reason (ref: ref from Wiki)
		#pragma omp for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
				Real p[3];
				info.pos(p, ix, iy, iz);
				p[0]-=position[0];
				p[1]-=position[1];
				p[2]-=position[2];
				const Real theta = std::atan2(p[1], p[0]);
				const Real r = (p[0]*p[0] + p[1]*p[1])/length;
				const Real vTheta = v_max*fac/r*(1-std::exp(-alpha*std::pow(r,2)))
				b(ix,iy,iz).U += - vTheta*std::sin(theta);
				b(ix,iy,iz).V +=   vTheta*std::cos(theta);
				b(ix,iy,iz).W = 0;
			}
		}
	}
	created = true;
}

void IF3D_VortexOperator::update(const int stepID, const Real t, const Real dt, const Real *Uinf)
{}

void IF3D_VortexOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_ObstacleOperator::_parseArguments(parser);
	parser.set_strict_mode();
	v_max = parser("-vmax").asDouble();
	parser.unset_strict_mode();
}

void IF3D_VortexOperator::computeVelocities(const Real* Uinf) {}
void IF3D_VortexOperator::computeForces(const int stepID, const Real time, const Real dt,
													const Real* Uinf, const Real NU, const bool bDump) {}
