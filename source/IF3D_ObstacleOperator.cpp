//
//  IF3D_MovingObstacleOperator.h
//  IF3D_ROCKS
//
//  Created by Wim van Rees on 06/10/14.
//
//

#include "IF3D_ObstacleOperator.h"
#include "GenericOperator.h"

struct ForcesOnSkin : public GenericLabOperator
{
    Real t;
    const Real NU, *vel_unit, *Uinf, *CM;
	int stencil_start[3], stencil_end[3];
	array<Real,19>* const measures;
	surfacePoints* const surfData;
	const map<int, pair<int, int>>* const surfaceBlocksFilter;
    std::map<int,ObstacleBlock*>* const obstacleBlocks;

    ForcesOnSkin(const Real NU, const Real* vel_unit, const Real* Uinf, const Real * CM,
    				map<int,ObstacleBlock*>* const obstacleBlocks,    	    //to read udef
					surfacePoints* const surface, 						    //most info I/O
					const map<int,pair<int,int>>* const surfaceBlocksFilter,//skip useless blocks
					array<Real,19>* const measures)     	                //additive quantities
	: t(0), NU(NU), vel_unit(vel_unit), Uinf(Uinf), CM(CM), measures(measures), surfData(surface),
	  surfaceBlocksFilter(surfaceBlocksFilter), obstacleBlocks(obstacleBlocks)
	{
    		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0, 1, 2);
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
	}

    ForcesOnSkin(const ForcesOnSkin& c):
    	t(0), NU(c.NU), vel_unit(c.vel_unit), Uinf(c.Uinf), CM(c.CM), measures(c.measures),
		surfData(c.surface), surfaceBlocksFilter(c.surfaceBlocksFilter), obstacleBlocks(c.obstacleBlocks)
    {
    	abort();
    	stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0, 1, 2);
    	stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
    	stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
    }

    template <typename Lab, typename BlockType>
	void operator()(Lab& lab, const BlockInfo& info, BlockType& b)
	{
    		const auto pos = surfaceBlocksFilter->find(info.blockID);
		if(pos == surfaceBlocksFilter->end()) return;
		const int first  = pos->second.first;
		const int second = pos->second.second;
		const Real _h3 = std::pow(info.h_gridpoint,3);
		const Real _1oH = NU / info.h_gridpoint; // 2 nu / 2 h
		for(int i=first; i<second; i++) { //i now is a fluid >element<
			assert(first<second && surfData->Set.size()>=second);
			Real p[3];
			const int ix = surfData->Set[i]->ix;
			const int iy = surfData->Set[i]->iy;
			const int iz = surfData->Set[i]->iz;
			const auto tempIt = obstacleBlocks->find(info.blockID);
			assert(tempIt != obstacleBlocks->end());
			info.pos(p, ix, iy, iz);
			//shear stresses
			const Real D11 =    _1oH*(lab(ix+1,iy,iz).u - lab(ix-1,iy,iz).u);
			const Real D22 =    _1oH*(lab(ix,iy+1,iz).v - lab(ix,iy-1,iz).v);
			const Real D33 =    _1oH*(lab(ix,iy,iz+1).w - lab(ix,iy,iz-1).w);
			const Real D12 = .5*_1oH*(lab(ix,iy+1,iz).u - lab(ix,iy-1,iz).u
									   +lab(ix+1,iy,iz).v - lab(ix-1,iy,iz).v);
			const Real D13 = .5*_1oH*(lab(ix,iy,iz+1).u - lab(ix,iy,iz-1).u
									   +lab(ix+1,iy,iz).w - lab(ix-1,iy,iz).w);
			const Real D23 = .5*_1oH*(lab(ix,iy+1,iz).w - lab(ix,iy-1,iz).w
									   +lab(ix,iy,iz+1).v - lab(ix,iy,iz-1).v);
			//normals computed with Towers 2009
			const Real normX = surfData->Set[i]->dchidx * _h3;
			const Real normY = surfData->Set[i]->dchidy * _h3;
			const Real normZ = surfData->Set[i]->dchidz * _h3;
			const Real fXV = D11 * normX + D12 * normY + D13 * normZ;
			const Real fYV = D12 * normX + D22 * normY + D23 * normZ;
			const Real fZV = D13 * normX + D23 * normY + D33 * normZ;
			const Real fXP = -b(ix,iy,iz).p * normX; //lab might not contain p!?
			const Real fYP = -b(ix,iy,iz).p * normY;
			const Real fZP = -b(ix,iy,iz).p * normZ;
			const Real fXT = fXV+fXP;
			const Real fYT = fYV+fYP;
			const Real fZT = fZV+fZP;
			surfData->P[i]  = b(ix,iy,iz).p;
			surfData->fX[i] = fXT;  surfData->fY[i] = fYT;  surfData->fZ[i] = fZT;
			surfData->fxP[i] = fXP; surfData->fyP[i] = fYP; surfData->fzP[i] = fZP;
			surfData->fxV[i] = fXV; surfData->fyV[i] = fYV; surfData->fzV[i] = fZV;
			surfData->pX[i] = p[0]; surfData->pY[i] = p[1]; surfData->pZ[i] = p[2];
			//perimeter:
			(*measures)[0] += surfData->Set[i]->delta * _h3;
			//forces (total, visc, pressure):
			(*measures)[1] += fXT; (*measures)[2] += fYT; (*measures)[3] += fZT;
			(*measures)[4] += fXP; (*measures)[5] += fYP; (*measures)[6] += fZP;
			(*measures)[7] += fXV; (*measures)[8] += fYV; (*measures)[9] += fZV;
			//torques:
			(*measures)[16] += (p[1]-CM[1])*fZV - (p[2]-CM[2])*fYT;
			(*measures)[17] += (p[2]-CM[2])*fXT - (p[0]-CM[0])*fZV;
			(*measures)[18] += (p[0]-CM[0])*fYT - (p[1]-CM[1])*fXT;
			//thrust, drag:
			const Real forcePar = fXT*vel_unit[0] + fYT*vel_unit[1] + fZT*vel_unit[2];
			(*measures)[10] += .5*(forcePar + std::abs(forcePar));
			(*measures)[11] -= .5*(forcePar - std::abs(forcePar));
			//save velocities in case of dump:
			surfData->vxDef[i] = tempIt->second->udef[iz][iy][ix][0];
			surfData->vyDef[i] = tempIt->second->udef[iz][iy][ix][1];
			surfData->vzDef[i] = tempIt->second->udef[iz][iy][ix][2];
			surfData->vx[i] = lab(ix,iy,iz).u + Uinf[0];
			surfData->vy[i] = lab(ix,iy,iz).v + Uinf[1];
			surfData->vz[i] = lab(ix,iy,iz).w + Uinf[2];
			//power output (and negative definite variant which ensures no elastic energy absorption)
			const Real powOut = fXT*surfData->vx[i] + fYT*surfData->vy[i] + fZT*surfData->vz[i];
			(*measures)[12] += powOut;
			(*measures)[13] += min((Real)0., powOut);
			//deformation power output (and negative definite variant which ensures no elastic energy absorption)
			const Real powDef = fXT*surfData->vxDef[i] + fYT*surfData->vyDef[i] + fZT*surfData->vzDef[i];
			(*measures)[14] += powDef;
			(*measures)[15] += min((Real)0., powDef);
		}
	}
};

void IF3D_ObstacleOperator::_computeUdefMoments(Real lin_momenta[3], Real ang_momenta[3], const Real CoM[3])
{
	const int N = vInfo.size();
	Real  V(0.0),  J0(0.0),  J1(0.0),  J2(0.0),  J3(0.0),  J4(0.0),  J5(0.0);
	Real  lm0(0.0),  lm1(0.0),  lm2(0.0); //linear momenta
	Real  am0(0.0),  am1(0.0),  am2(0.0); //angular momenta

#pragma omp parallel for schedule(static) reduction(+:V,lm0,lm1,lm2,J0,J1,J2,J3,J4,J5,am0,am1,am2)
	for(int i=0; i<vInfo.size(); i++) {
		BlockInfo info = vInfo[i];
		const auto pos = obstacleBlocks.find(info.blockID);
		if(pos == obstacleBlocks.end()) continue;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			const Real Xs = pos->second->chi[iz][iy][ix];
			if (Xs == 0) continue;
			Real p[3];
			info.pos(p, ix, iy, iz);
			p[0]-=CoM[0];
			p[1]-=CoM[1];
			p[2]-=CoM[2];
			V   += Xs;
			lm0 += Xs *       pos->second->udef[iz][iy][ix][0];
			lm1 += Xs *       pos->second->udef[iz][iy][ix][1];
			lm2 += Xs *       pos->second->udef[iz][iy][ix][2];
			am0 += Xs * (p[1]*pos->second->udef[iz][iy][ix][2] - p[2]*pos->second->udef[iz][iy][ix][1]);
			am1 += Xs * (p[2]*pos->second->udef[iz][iy][ix][0] - p[0]*pos->second->udef[iz][iy][ix][2]);
			am2 += Xs * (p[0]*pos->second->udef[iz][iy][ix][1] - p[1]*pos->second->udef[iz][iy][ix][0]);
			J0  += Xs * (p[1]*p[1]+p[2]*p[2]);
			J1  += Xs * (p[0]*p[0]+p[2]*p[2]);
			J2  += Xs * (p[0]*p[0]+p[1]*p[1]);
			J3  -= Xs *  p[0]*p[1];
			J4  -= Xs *  p[0]*p[2];
			J5  -= Xs *  p[1]*p[2];
		}
	}

    double globals[13];
    double locals[13] = {lm0,lm1,lm2,am0,am1,am2,J0,J1,J2,J3,J4,J5,V};
	MPI::COMM_WORLD.Allreduce(locals, globals, 13, MPI::DOUBLE, MPI::SUM);

    assert(globals[12] > std::numeric_limits<double>::epsilon());
    const Real dv = std::pow(vInfo[0].h_gridpoint,3);
    lin_momenta[0] = globals[0]/globals[12];
    lin_momenta[1] = globals[1]/globals[12];
    lin_momenta[2] = globals[2]/globals[12];
    Real tmp_AV[3] = {
    		globals[3]* dv,
			globals[4]* dv,
			globals[5]* dv
    };
    Real tmp_J[6] = {
    		globals[6] * dv,
			globals[7] * dv,
			globals[8] * dv,
			globals[9] * dv,
			globals[10]* dv,
			globals[11]* dv
    };

    _finalizeAngVel(ang_momenta, tmp_J, tmp_AV[0], tmp_AV[1], tmp_AV[1]);

    if(bFixToPlanar) { //TODO: why this step?
    		ang_momenta[0] = ang_momenta[1] = 0.0;
    		ang_momenta[2] = globals[5]/globals[8]; // av2/j2
    }
}

void IF3D_ObstacleOperator::_makeDefVelocitiesMomentumFree(const Real CoM[3])
{
	_computeUdefMoments(transVel_correction, angVel_correction, CoM);
    printf("Correction of lin: %f %f ang: %f\n", transVel_correction[0], transVel_correction[1], angVel_correction[2]);

#pragma omp parallel for schedule(static)
    for(int i=0; i<vInfo.size(); i++) {
        BlockInfo info = vInfo[i];
        const auto pos = obstacleBlocks.find(info.blockID);
        if(pos == obstacleBlocks.end()) continue;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
            Real p[3];
            info.pos(p, ix, iy, iz);
            p[0]-=CoM[0];
            p[1]-=CoM[1];
            p[2]-=CoM[2];
            const Real correctVel[3] = {
				transVel_correction[0] + (angVel_correction[1]*p[2] - angVel_correction[2]*p[1]),
				transVel_correction[1] + (angVel_correction[2]*p[0] - angVel_correction[0]*p[2]),
				transVel_correction[2] + (angVel_correction[0]*p[1] - angVel_correction[1]*p[0])
            };
            pos->second->udef[iz][iy][ix][0] -= correctVel[0];
            pos->second->udef[iz][iy][ix][1] -= correctVel[1];
            pos->second->udef[iz][iy][ix][2] -= correctVel[2];
        }
    }

#ifndef NDEBUG
    Real dummy_ang[3], dummy_lin[3];
    _computeUdefMoments(dummy_lin, dummy_ang, CoM);
#endif
}

void IF3D_ObstacleOperator::_parseArguments(ArgumentParser & parser)
{
    parser.set_strict_mode();
    position[0] = parser("-xpos").asDouble();
    position[1] = parser("-ypos").asDouble();
    position[2] = parser("-zpos").asDouble();
    length = parser("-L").asDouble();
    parser.unset_strict_mode();
    quaternion[0] = parser("-quat0").asDouble(1.0);
    quaternion[1] = parser("-quat1").asDouble(0.0);
    quaternion[2] = parser("-quat2").asDouble(0.0);
    quaternion[3] = parser("-quat3").asDouble(0.0);
    const Real one = sqrt(quaternion[0]*quaternion[0]
						   +quaternion[1]*quaternion[1]
						   +quaternion[2]*quaternion[2]
						   +quaternion[3]*quaternion[3]);

    if(fabs(one-1.0) > 5*numeric_limits<Real>::epsilon()) {
    	printf("Parsed quaternion length is not equal to one. It really ought to be.\n");
    	abort();
    }
    if(length < 5*numeric_limits<Real>::epsilon()) {
    	printf("Parsed length is equal to zero. It really ought not to be.\n");
    	abort();
    }

    bFixToPlanar = parser("-bFixToPlanar").asBool(false);
    bFixFrameOfRef = parser("-bFixFrameOfRef").asBool(false);
}

void IF3D_ObstacleOperator::_writeComputedVelToFile(const int step_id, const Real t, const Real* Uinf)
{
	if(rank!=0) return;
	stringstream ssR;
	ssR<<"computedVelocity_"<<obstacleID<<".dat";
    std::ofstream savestream(ssR.str(), ios::out | ios::app);
    const std::string tab("\t");
    
    if(step_id==0)
        savestream << "step" << tab << "time" << tab << "CMx" << tab << "CMy" << tab << "CMz" << tab
				   << "quat_0" << tab << "quat_1" << tab << "quat_2" << tab << "quat_3" << tab
				   << "vel_x" << tab << "vel_y" << tab << "vel_z" << tab
				   << "angvel_x" << tab << "angvel_y" << tab << "angvel_z" << tab << "volume" << tab
				   << "J0" << tab << "J1" << tab << "J2" << tab << "J3" << tab << "J4" << tab << "J5" << std::endl;
    
    savestream << step_id << tab;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<float>::digits10 + 1);
    savestream << t << tab << position[0] << tab << position[1] << tab << position[2] << tab
    		   << quaternion[0] << tab << quaternion[1] << tab << quaternion[2] << tab << quaternion[3] << tab
			   << transVel[0]-Uinf[0] << tab << transVel[1]-Uinf[1] << tab << transVel[2]-Uinf[2] << tab
			   << angVel[0] << tab << angVel[1] << tab << angVel[2] << tab << volume << tab
			   << J[0] << tab << J[1] << tab << J[2] << tab << J[3] << tab << J[4] << tab << J[5] << std::endl;
    savestream.close();
}

void IF3D_ObstacleOperator::_writeDiagForcesToFile(const int step_id, const Real t)
{
	if(rank!=0) return;
	stringstream ssR;
	ssR<<"diagnosticsForces_"<<obstacleID<<".dat";
    std::ofstream savestream(ssR.str(), ios::out | ios::app);
    const std::string tab("\t");
    
    if(step_id==0)
        savestream << "step" << tab << "time" << tab << "mass" << tab
				   << "force_x" << tab << "force_y" << tab << "force_z" << tab
				   << "torque_x" << tab << "torque_y" << tab << "torque_z" << std::endl;
    
    savestream << step_id << tab;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<float>::digits10 + 1);
    savestream << t << tab << mass << tab << force[0] << tab << force[1] << tab << force[2] << tab
    		   << torque[0] << tab << torque[1] << tab << torque[2] << std::endl;
    savestream.close();
}

void IF3D_ObstacleOperator::computeDiagnostics(const int stepID, const Real time, const Real* Uinf, const Real lambda)
{
	Real CM[3];
	this->getCenterOfMass(CM);
    const int N = vInfo.size();
    Real _area(0.0), _forcex(0.0), _forcey(0.0), _forcez(0.0), _torquex(0.0), _torquey(0.0), _torquez(0.0);
    Real garea(0.0), gforcex(0.0), gforcey(0.0), gforcez(0.0), gtorquex(0.0), gtorquey(0.0), gtorquez(0.0);

    #pragma omp parallel for schedule(static) reduction(+:_area,_forcex,_forcey,_forcez,_torquex,_torquey,_torquez)
    for(int i=0; i<vInfo.size(); i++) {
        BlockInfo info = vInfo[i];
        std::map<int,ObstacleBlock*>::const_iterator pos = obstacleBlocks.find(info.blockID);
        if(pos == obstacleBlocks.end()) continue;
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
            const Real Xs = pos->second->chi[iz][iy][ix];
            if (Xs == 0) continue;
            Real p[3];
            info.pos(p, ix, iy, iz);
            p[0]-=CM[0];
            p[1]-=CM[1];
            p[2]-=CM[2];
		    const Real object_UR[3] = {
				angVel[1]*p[2]-angVel[2]*p[1],
				angVel[2]*p[0]-angVel[0]*p[2],
				angVel[0]*p[1]-angVel[1]*p[0]
		    };
            const Real object_UDEF[3] = {
                pos->second->udef[iz][iy][ix][0],
                pos->second->udef[iz][iy][ix][1],
                pos->second->udef[iz][iy][ix][2]
            };
            const Real U[3] = {
                b(ix,iy,iz).u+Uinf[0] - (transVel[0]+object_UR[0]+object_UDEF[0]),
                b(ix,iy,iz).v+Uinf[1] - (transVel[1]+object_UR[1]+object_UDEF[1]),
                b(ix,iy,iz).w+Uinf[2] - (transVel[2]+object_UR[2]+object_UDEF[2])
            };
            _area += Xs;
            _forcex += U[0]*Xs;
            _forcey += U[1]*Xs;
            _forcez += U[2]*Xs;
            _torquex += (p[1]*U[2]-p[2]*U[1])*Xs;
            _torquey += (p[2]*U[0]-p[0]*U[2])*Xs;
            _torquez += (p[0]*U[1]-p[1]*U[0])*Xs;
        }
    }

    double globals[7];
    double locals[7] = {_forcex,_forcey,_forcez,_torquex,_torquey,_torquez,_area};
    MPI::COMM_WORLD.Allreduce(locals, globals, 7, MPI::DOUBLE, MPI::SUM);

    const Real dV = std::pow(vInfo[0].h_gridpoint, 3);
    force[0]  = globals[0]*dV*lambda;
    force[1]  = globals[1]*dV*lambda;
    force[2]  = globals[2]*dV*lambda;
    torque[0] = globals[3]*dV*lambda;
    torque[1] = globals[4]*dV*lambda;
    torque[2] = globals[5]*dV*lambda;
    mass      = globals[6]*dV;
    
    _writeDiagForcesToFile(stepID, time);
}

void IF3D_ObstacleOperator::computeVelocities(const Real* Uinf)
{
    Real CM[3];
    this->getCenterOfMass(CM);
    const int N = vInfo.size();

    //local variables
    Real V(0.0), J0(0.0), J1(0.0), J2(0.0), J3(0.0), J4(0.0), J5(0.0);
    Real lm0(0.0), lm1(0.0), lm2(0.0); //linear momenta
    Real am0(0.0), am1(0.0), am2(0.0); //angular momenta

#pragma omp parallel for schedule(static) reduction(+:V,lm0,lm1,lm2,J0,J1,J2,J3,J4,J5,am0,am1,am2)
    for(int i=0; i<vInfo.size(); i++) {
        BlockInfo info = vInfo[i];
        FluidBlock & b = *(FluidBlock*)info.ptrBlock;
        const auto pos = obstacleBlocks.find(info.blockID);
        if(pos == obstacleBlocks.end()) continue;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
            const Real Xs = pos->second->chi[iz][iy][ix];
            if (Xs == 0) continue;
            Real p[3];
            info.pos(p, ix, iy, iz);
            p[0]-=CM[0];
            p[1]-=CM[1];
            p[2]-=CM[2];
            V     += Xs;
            lm0   += Xs * b(ix,iy,iz).u;
            lm1   += Xs * b(ix,iy,iz).v;
            lm2   += Xs * b(ix,iy,iz).w;
            am0 += Xs*(p[1]*b(ix,iy,iz).w - p[2]*b(ix,iy,iz).v);
            am1 += Xs*(p[2]*b(ix,iy,iz).u - p[0]*b(ix,iy,iz).w);
            am2 += Xs*(p[0]*b(ix,iy,iz).v - p[1]*b(ix,iy,iz).u);
            J0+=Xs*(p[1]*p[1]+p[2]*p[2]);
            J1+=Xs*(p[0]*p[0]+p[2]*p[2]);
            J2+=Xs*(p[0]*p[0]+p[1]*p[1]);
            J3-=Xs*p[0]*p[1];
            J4-=Xs*p[0]*p[2];
            J5-=Xs*p[1]*p[2];
        }
    }

    double globals[13];
    double locals[13] = {lm0,lm1,lm2,am0,am1,am2,J0,J1,J2,J3,J4,J5,V};
	MPI::COMM_WORLD.Allreduce(locals, globals, 13, MPI::DOUBLE, MPI::SUM);

    assert(globals[12] > std::numeric_limits<double>::epsilon());
    const Real dv = std::pow(vInfo[0].h_gridpoint,3);
    transVel[0] = globals[0]/globals[12] + Uinf[0];
    transVel[1] = globals[1]/globals[12] + Uinf[1];
    transVel[2] = globals[2]/globals[12] + Uinf[2];
    Real tmp_AV[3] = {
    		globals[3]* dv,
			globals[4]* dv,
			globals[5]* dv
    };
    J[0] = globals[6] * dv;
    J[1] = globals[7] * dv;
    J[2] = globals[8] * dv;
    J[3] = globals[9] * dv;
    J[4] = globals[10]* dv;
    J[5] = globals[11]* dv;
    volume = globals[12] * dv;

    _finalizeAngVel(angVel, J, tmp_AV[0], tmp_AV[1], tmp_AV[1]);

    if(bFixToPlanar) {
        transVel[2] = 0.0;
        angVel[0] = angVel[1] = 0.0;
        angVel[2] = globals[5]/globals[8]; // av2/j2
    }
}

void IF3D_ObstacleOperator::_finalizeAngVel(Real AV[3], const Real _J[6], const Real& gam0, const Real& gam1, const Real& gam2)
{
	// try QR factorization to avoid dealing with determinant
	const Real u1[3] = {_J[0], _J[3], _J[4]};
	const Real magu1sq = u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2];
	// subtract projection of a2 onto u1
	const Real proj1 = u1[0]*_J[3] + u1[1]*_J[1] + u1[2]*_J[5];
	const Real u2[3] = {
			_J[3] - proj1*u1[0]/magu1sq,
			_J[1] - proj1*u1[1]/magu1sq,
			_J[5] - proj1*u1[2]/magu1sq
	};
	const Real magu2sq = u2[0]*u2[0] + u2[1]*u2[1] + u2[2]*u2[2];
	// subtract projection of a3 onto u1
	const Real proj2 = u1[0]*_J[4] + u1[1]*_J[5] + u1[2]*_J[2];
	const Real u3_tmp[3] = {
			_J[4] - proj2*u1[0]/magu1sq,
			_J[5] - proj2*u1[1]/magu1sq,
			_J[2] - proj2*u1[2]/magu1sq
	};
	// subtract projection of u3_tmp onto u2
	const Real proj3 = u2[0]*u3_tmp[0] + u2[1]*u3_tmp[1] + u2[2]*u3_tmp[2];
	const Real u3[3] = {
			u3_tmp[0] - proj3*u2[0]/magu2sq,
			u3_tmp[1] - proj3*u2[1]/magu2sq,
			u3_tmp[2] - proj3*u2[2]/magu2sq
	};
	const Real magu3sq = u3[0]*u3[0] + u3[1]*u3[1] + u3[2]*u3[2];
	const Real magu1 = std::sqrt(magu1sq);
	const Real magu2 = std::sqrt(magu2sq);
	const Real magu3 = std::sqrt(magu3sq);
	const Real Q[3][3] = {
			{u1[0]/magu1, u2[0]/magu2, u3[0]/magu3},
			{u1[1]/magu1, u2[1]/magu2, u3[1]/magu3},
			{u1[2]/magu1, u2[2]/magu2, u3[2]/magu3}
	};
	// find out if Q is orthogonal
	const Real R[3][3] = {
			{Q[0][0]*_J[0] + Q[1][0]*_J[3] + Q[2][0]*_J[4], Q[0][0]*_J[3] + Q[1][0]*_J[1] + Q[2][0]*_J[5], Q[0][0]*_J[4] + Q[1][0]*_J[5] + Q[2][0]*_J[2]},
			{Q[0][1]*_J[0] + Q[1][1]*_J[3] + Q[2][1]*_J[4], Q[0][1]*_J[3] + Q[1][1]*_J[1] + Q[2][1]*_J[5], Q[0][1]*_J[4] + Q[1][1]*_J[5] + Q[2][1]*_J[2]},
			{Q[0][2]*_J[0] + Q[1][2]*_J[3] + Q[2][2]*_J[4], Q[0][2]*_J[3] + Q[1][2]*_J[1] + Q[2][2]*_J[5], Q[0][2]*_J[4] + Q[1][2]*_J[5] + Q[2][2]*_J[2]}
	};
	// d = Q^T b
	const Real d[3] = {
			Q[0][0]*gam0 + Q[1][0]*gam1 + Q[2][0]*gam2,
			Q[0][1]*gam0 + Q[1][1]*gam1 + Q[2][1]*gam2,
			Q[0][2]*gam0 + Q[1][2]*gam1 + Q[2][2]*gam2,
	};
	// bwd subtitution: R x = d
	AV[2] = d[2]/R[2][2];
	AV[1] = (d[1] - R[1][2]*angVel[2])/R[1][1];
	AV[0] = (d[0] - R[0][1]*angVel[1] - R[0][2]*angVel[2])/R[0][0];
}

void IF3D_ObstacleOperator::computeForces(const int stepID, const Real time, const Real* Uinf, const Real NU, const bool bDump)
{   //TODO: improve dumping: gather arrays before writing to file
	//TODO: this kernel does an illegal read on the grid. Guack-a-bug!
    //ugly piece of code that creates a map that allows us to skip the non-surface blocks
    map<int, pair<int, int>> surfaceBlocksFilter;
    {
		vector<int> usefulIDs; //which blocks are occupied by surface
		vector<int> firstInfo; //first element of surfData occupied by each surf. block
							   //(the elems of surfData are grouped together by blockID)
		for(int i=0; i<surfData.Ndata; i++) {
			bool unique = true;
			for(int k=0; k<usefulIDs.size(); k++)
				if (surfData.Set[i]->blockID == usefulIDs[k])
					{ unique = false; break; } //Ive already seen that block
			if (unique) {
				usefulIDs.push_back(surfData.Set[i]->blockID);
				firstInfo.push_back(i);
			}
		}
		firstInfo.push_back(surfData.Ndata);
		for(int i=0; i<usefulIDs.size(); i++) { //now i arrange them in a map because Im lazy
			assert(surfaceBlocksFilter.find(usefulIDs[i]) == surfaceBlocksFilter.end());
			surfaceBlocksFilter[usefulIDs[i]] = make_pair(firstInfo[i],firstInfo[i+1]);

		//const auto pos = surfaceBlocksFilter.find(usefulIDs[i]);
		//printf("%d block id %d (%d) occupies from %d to %d\n",rank,usefulIDs[i],pos->first,pos->second.first,pos->second.second);
		}
    }
    
	Real CM[3];
    this->getCenterOfMass(CM);
    const int nthreads = omp_get_max_threads();
    vector<array<Real,19>> partialSums(nthreads);
    Real vel_unit[3] = {0., 0., 0.};
    const Real velx_tot = transVel[0]-Uinf[0];
    const Real vely_tot = transVel[1]-Uinf[1];
    const Real velz_tot = transVel[2]-Uinf[2];
    const Real vel_norm = std::sqrt(velx_tot*velx_tot + vely_tot*vely_tot + velz_tot*velz_tot);
    if (vel_norm>1e-9) {
        vel_unit[0] = velx_tot/vel_norm;
        vel_unit[1] = vely_tot/vel_norm;
        vel_unit[2] = velz_tot/vel_norm;
    }

    vector<ForcesOnSkin> finalize;
    for(int i = 0; i < nthreads; ++i)
    	finalize.push_back(ForcesOnSkin(NU,vel_unit,Uinf,CM,&obstacleBlocks,&surfData,&surfaceBlocksFilter,&partialSums[i]));

    compute(finalize);
    double localSum[19], globalSum[19];
    for(int i=0; i<nthreads; i++)
    	for(int j=0; j<19; j++)
    		localSum[j] += (double)partialSums[i][j];
    MPI::COMM_WORLD.Allreduce(localSum, globalSum, 19, MPI::DOUBLE, MPI::SUM);

    //additive quantities:
    totChi      = globalSum[0];
    surfForce[0]= globalSum[1];  surfForce[1] = globalSum[2]; surfForce[2]=globalSum[3];
    drag	 	= globalSum[10]; thrust		  = globalSum[11];
    Pout     	= globalSum[12]; PoutBnd      = globalSum[13];
    defPower 	= globalSum[14]; defPowerBnd  = globalSum[15];
    //derived quantities:
    Pthrust    = thrust*vel_norm;
    Pdrag      =   drag*vel_norm;
    EffPDef    = Pthrust/(Pthrust-min(defPower,(Real)0.));
    EffPDefBnd = Pthrust/(Pthrust-    defPowerBnd);
    
    if (bDump)  surfData.print(obstacleID, stepID, rank);

    if(rank==0)
    {
        ofstream fileForce;
        ofstream filePower;
    	stringstream ssF, ssP;
    	ssF<<"forceValues_"<<obstacleID<<".txt";
    	ssP<<"powerValues_"<<obstacleID<<".txt";

        fileForce.open(ssF.str().c_str(), ios::app);
        fileForce<<time<<" "<<surfForce[0] <<" "<<surfForce[1] <<" "<<surfForce[2] <<" "
						    <<globalSum[4] <<" "<<globalSum[5] <<" "<<globalSum[6] <<" "
        		 		    <<globalSum[7] <<" "<<globalSum[8] <<" "<<globalSum[9] <<" "
        		 		    <<globalSum[16]<<" "<<globalSum[17]<<" "<<globalSum[18]<<" "
						    << drag  <<" "<< thrust<<" "<< totChi<<endl;
        fileForce.close();
        filePower.open(ssP.str().c_str(), ios::app);
        filePower<<time<<" "<<Pthrust<<" "<<Pdrag<<" "<<PoutBnd<<" "<<Pout<<" "
        		 <<defPowerBnd<<" "<<defPower<<" "<<EffPDefBnd<<" "<<EffPDef<<endl;
        filePower.close();
    }
}

void IF3D_ObstacleOperator::update(const int step_id, const Real t, const Real dt, const Real* Uinf)
{
    position[0] += dt*transVel[0];
    position[1] += dt*transVel[1];
    position[2] += dt*transVel[2];
    const Real dqdt[4] = {
        0.5*(-angVel[0]*quaternion[1]-angVel[1]*quaternion[2]-angVel[2]*quaternion[3]),
        0.5*( angVel[0]*quaternion[0]+angVel[1]*quaternion[3]-angVel[2]*quaternion[2]),
        0.5*(-angVel[0]*quaternion[3]+angVel[1]*quaternion[0]+angVel[2]*quaternion[1]),
        0.5*( angVel[0]*quaternion[2]-angVel[1]*quaternion[1]+angVel[2]*quaternion[0])
    };

    // normality preserving advection (Simulation of colliding constrained rigid bodies - Kleppmann 2007 Cambridge University, p51)
    // move the correct distance on the quaternion unit ball surface, end up with normalized quaternion
    const Real deltaq[4] = {
        dqdt[0]*dt,
        dqdt[1]*dt,
        dqdt[2]*dt,
        dqdt[3]*dt
    };

    const Real deltaq_length = std::sqrt(deltaq[0]*deltaq[0]+deltaq[1]*deltaq[1]+deltaq[2]*deltaq[2]+deltaq[3]*deltaq[3]);
    if(deltaq_length>std::numeric_limits<Real>::epsilon()) {
        const Real tanfac = std::tan(deltaq_length)/deltaq_length;
        const Real num[4] = {
            quaternion[0]+tanfac*deltaq[0],
            quaternion[1]+tanfac*deltaq[1],
            quaternion[2]+tanfac*deltaq[2],
            quaternion[3]+tanfac*deltaq[3],
        };

        const Real invDenum = 1./(std::sqrt(num[0]*num[0]+num[1]*num[1]+num[2]*num[2]+num[3]*num[3]));
        quaternion[0] = num[0]*invDenum;
        quaternion[1] = num[1]*invDenum;
        quaternion[2] = num[2]*invDenum;
        quaternion[3] = num[3]*invDenum;
    }

#ifndef NDEBUG
    if(rank==0) {
    	std::cout << "POSITION INFO AFTER UPDATE T, DT: " << t << " " << dt << std::endl;
    	std::cout << "POS: " << position[0] << " " << position[1] << " " << position[2] << std::endl;
    	std::cout << "TVL: " << transVel[0] << " " << transVel[1] << " " << transVel[2] << std::endl;
    	std::cout << "QUT: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
    	std::cout << "AVL: " << angVel[0] << " " << angVel[1] << " " << angVel[2] << std::endl;
    }
    const Real q_length=std::sqrt(quaternion[0]*quaternion[0]
								   +quaternion[1]*quaternion[1]
								   +quaternion[2]*quaternion[2]
								   +quaternion[3]*quaternion[3]);
    assert(std::abs(q_length-1.0) < 5*std::numeric_limits<Real>::epsilon());
#endif

    _writeComputedVelToFile(step_id, t, Uinf);
}

void IF3D_ObstacleOperator::characteristic_function()
{
#pragma omp parallel
	{
#pragma omp for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			BlockInfo info = vInfo[i];
			std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(info.blockID);

			if(pos != obstacleBlocks.end()) {
                FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
					b(ix,iy,iz).chi = std::max(pos->second->chi[iz][iy][ix], b(ix,iy,iz).chi);
			}
		}
	}
}

std::vector<int> IF3D_ObstacleOperator::intersectingBlockIDs(const int buffer) const
{
	assert(buffer <= 2); // only works for 2: if different definition of deformation blocks, implement your own
	std::vector<int> retval;
	const int N = vInfo.size();

	for(int i=0; i<N; i++) {
		BlockInfo info = vInfo[i];
		std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(info.blockID);
		if(pos != obstacleBlocks.end()) retval.push_back(info.blockID);
	}
	return retval;
}

void IF3D_ObstacleOperator::getTranslationVelocity(Real UT[3]) const
{
    UT[0]=transVel[0];
    UT[1]=transVel[1];
    UT[2]=transVel[2];
}

void IF3D_ObstacleOperator::setTranslationVelocity(Real UT[3])
{
    transVel[0]=UT[0];
    transVel[1]=UT[1];
    transVel[2]=UT[2];
}

void IF3D_ObstacleOperator::getAngularVelocity(Real W[3]) const
{
    W[0]=angVel[0];
    W[1]=angVel[1];
    W[2]=angVel[2];
}

void IF3D_ObstacleOperator::setAngularVelocity(const Real W[3])
{
	angVel[0]=W[0];
	angVel[1]=W[1];
	angVel[2]=W[2];
}

void IF3D_ObstacleOperator::getCenterOfMass(Real CM[3]) const
{
    CM[0]=position[0];
    CM[1]=position[1];
    CM[2]=position[2];
}

void IF3D_ObstacleOperator::save(const int step_id, const Real t, std::string filename)
{
	if(rank!=0) return;
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename+".txt");
    savestream << t << std::endl;
    savestream << position[0] << "\t" << position[1] << "\t" << position[2] << std::endl;
    savestream << quaternion[0] << "\t" << quaternion[1] << "\t" << quaternion[2] << "\t" << quaternion[3] << std::endl;
    savestream << transVel[0] << "\t" << transVel[1] << "\t" << transVel[2] << std::endl;
    savestream << angVel[0] << "\t" << angVel[1] << "\t" << angVel[2] << std::endl;
}

void IF3D_ObstacleOperator::restart(const Real t, std::string filename)
{
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    
    Real restart_time;
    restartstream >> restart_time;
    assert(std::abs(restart_time-t) < 1e-9);
    
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream.close();
    
    {
        std::cout << "RESTARTED BODY: " << std::endl;
        std::cout << "TIME: \t" << restart_time << std::endl;
        std::cout << "POS : \t" << position[0] << " " << position[1] << " " << position[2] << std::endl;
        std::cout << "ANGLE:\t" << quaternion[0] << " " << quaternion[1] << " " << quaternion[2] << " " << quaternion[3] << std::endl;
        std::cout << "TVEL: \t" << transVel[0] << " " << transVel[1] << " " << transVel[2] << std::endl;
        std::cout << "AVEL: \t" << angVel[0] << " " << angVel[1] << " " << angVel[2] << std::endl;
    }
}

void IF3D_ObstacleOperator::Accept(ObstacleVisitor * visitor)
{
	visitor->visit(this);
}
