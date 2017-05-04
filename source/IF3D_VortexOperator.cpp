//
//  IF3D_FishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_VortexOperator.h"
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
	if(time<8. || created) return;
	printf("Crating a vortex in %f %f %f core size %f max vel %f\n", position[0],position[1],position[2],length,v_max);
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
				const Real r = std::sqrt(p[0]*p[0] + p[1]*p[1])/length;
				//if(r>3) continue;
				const Real z = std::fabs(p[2])/0.05;
				//if(z>3) continue;
				const Real arg = std::max(-20., -alpha*r*r -z*z);
				const Real vTheta = 0.1*fac/(r+2.2e-16)*(1-std::exp(arg));
				b(ix,iy,iz).u -= vTheta*std::sin(theta);
				b(ix,iy,iz).v += vTheta*std::cos(theta);
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
void IF3D_VortexOperator::Accept(ObstacleVisitor * visitor) {}
void IF3D_VortexOperator::finalize(const int step_id,const Real time,
                                                const Real dt, const Real *Uinf) {}
void IF3D_VortexOperator::execute(Communicator * comm, const int iAgent, const Real time, const int iLabel) {}
void IF3D_VortexOperator::interpolateOnSkin(const Real time, const int stepID, bool dumpWake) {}
void IF3D_VortexOperator::getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL, Real*& pYL, Real*& pXU, Real*& pYU, int& Npts)
{
	pXL = pYL = pXU = pYU = nullptr;
	Npts = 0;
}
void IF3D_VortexOperator::computeDiagnostics(const int stepID, const Real time, const Real* Uinf, const Real lambda) {}
