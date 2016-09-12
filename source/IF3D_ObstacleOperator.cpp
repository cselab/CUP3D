//
//  IF3D_MovingObstacleOperator.h
//  IF3D_ROCKS
//
//  Created by Wim van Rees on 06/10/14.
//
//

#include "IF3D_ObstacleOperator.h"

void IF3D_ObstacleOperator::_makeDefVelocitiesMomentumFree(const double CoM[3])
{
	BlockInfo * ary = &vInfo.front();
	const int N = vInfo.size();
	double CM[3];
	this->getCenterOfMass(CM);
	const double dv = vInfo[0].h_gridpoint*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint;

	//local variables
	double V(0.0), J0(0.0), J1(0.0), J2(0.0), J3(0.0), J4(0.0), J5(0.0);
	double lm0(0.0), lm1(0.0), lm2(0.0); //linear momenta
	double am0(0.0), am1(0.0), am2(0.0); //angular momenta
	//global allReduce
	double gV(0.0), gJ0(0.0), gJ1(0.0), gJ2(0.0), gJ3(0.0), gJ4(0.0), gJ5(0.0);
	double glm0(0.0), glm1(0.0), glm2(0.0); //linear momenta
	double gam0(0.0), gam1(0.0), gam2(0.0); //angular momenta

#pragma omp parallel for schedule(static) reduction(+:V,lm0,lm1,lm2,J0,J1,J2,J3,J4,J5,am0,am1,am2)
	for(int i=0; i<vInfo.size(); i++) {
		std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(i);
		if(pos == obstacleBlocks.end()) continue;

		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			const double Xs = pos->second->chi[iz][iy][ix];
			if (Xs == 0) continue;
			double p[3];
			info.pos(p, ix, iy, iz);
			p[0]-=CM[0];
			p[1]-=CM[1];
			p[2]-=CM[2];

			const double u[3] = {
					pos->second->udef[iz][iy][ix][0],
					pos->second->udef[iz][iy][ix][1],
					pos->second->udef[iz][iy][ix][2]
			};

			V     += Xs;
			lm0   += Xs * u[0];
			lm1   += Xs * u[1];
			lm2   += Xs * u[2];
			am0 += Xs*(p[1]*u[2] - p[2]*u[1]);
			am1 += Xs*(p[2]*u[0] - p[0]*u[2]);
			am2 += Xs*(p[0]*u[1] - p[1]*u[0]);

			J0+=Xs*(p[1]*p[1]+p[2]*p[2]);
			J1+=Xs*(p[0]*p[0]+p[2]*p[2]);
			J2+=Xs*(p[0]*p[0]+p[1]*p[1]);
			J3-=Xs*p[0]*p[1];
			J4-=Xs*p[0]*p[2];
			J5-=Xs*p[1]*p[2];
		}
	}
	//TODO: faster to place these vars into an array, single reduce and unpack?
	MPI::COMM_WORLD.Allreduce(&lm0, &glm0, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&lm1, &glm1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&lm2, &glm2, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&am0, &gam0, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&am1, &gam1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&am2, &gam2, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J0, &gJ0, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J1, &gJ1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J2, &gJ2, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J3, &gJ3, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J4, &gJ4, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J5, &gJ5, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&V, &gV, 1, MPI::DOUBLE, MPI::SUM);

	transVel_correction[0] = glm0/gV;
	transVel_correction[1] = glm1/gV;
	transVel_correction[2] = glm2/gV;
    const Real _J[6] = {gJ0, gJ1, gJ2, gJ3, gJ4, gJ5};
    _finalizeAngVel(angVel_correction, _J, gam0, gam1, gam2);
    if(bFixToPlanar) { //TODO: why this step?
    	angVel_correction[0] = angVel_correction[1] = 0.0;
    	angVel_correction[2] = gam2/gJ2;
    }

#pragma omp parallel for schedule(static)
    for(int i=0; i<vInfo.size(); i++) {
        std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(i);
        if(pos == obstacleBlocks.end()) continue;
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)ary[i].ptrBlock;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
            double p[3];
            info.pos(p, ix, iy, iz);
            p[0]-=CoM[0];
            p[1]-=CoM[1];
            p[2]-=CoM[2];

            const double correctVel[3] = {
            		transVel_correction[0] + (angVel_correction[1]*p[2] - angVel_correction[2]*p[1]),
					transVel_correction[1] + (angVel_correction[2]*p[0] - angVel_correction[0]*p[2]),
					transVel_correction[2] + (angVel_correction[0]*p[1] - angVel_correction[1]*p[0])
            };

            pos->second->udef[iz][iy][ix][0] -= correctVel[0];
            pos->second->udef[iz][iy][ix][1] -= correctVel[1];
            pos->second->udef[iz][iy][ix][2] -= correctVel[2];
        }
    }
}

void IF3D_ObstacleOperator::_parseArguments(ArgumentParser & parser)
{
    parser.set_strict_mode();
    const Real xpos = parser("-xpos").asDouble();
    parser.unset_strict_mode();
    const Real ypos = parser("-ypos").asDouble(0.5);
    const Real zpos = parser("-zpos").asDouble(0.5);

    quaternion[0] = parser("-quat0").asDouble(0.0);
    quaternion[1] = parser("-quat1").asDouble(0.0);
    quaternion[2] = parser("-quat2").asDouble(0.0);
    quaternion[3] = parser("-quat3").asDouble(0.0);

    bFixToPlanar = parser("-bFixToPlanar").asBool(false);
    
    position[0] = xpos;
    position[1] = ypos;
    position[2] = zpos;
}

void IF3D_ObstacleOperator::_writeComputedVelToFile(const int step_id, const Real t)
{
    const std::string fname = "computedVelocity_"+std::to_string(obstacleID)+".dat";
    std::ofstream savestream(fname, ios::out | ios::app);
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
			   << transVel[0] << tab << transVel[1] << tab << transVel[2] << tab
			   << angVel[0] << tab << angVel[1] << tab << angVel[2] << tab << volume << tab
			   << J[0] << tab << J[1] << tab << J[2] << tab << J[3] << tab << J[4] << tab << J[5] << std::endl;
    savestream.close();
}

void IF3D_ObstacleOperator::_writeDiagForcesToFile(const int step_id, const Real t)
{
    const std::string forcefilename = "diagnosticsForces_"+std::to_string(obstacleID)+".dat";
    std::ofstream savestream(forcefilename, ios::out | ios::app);
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

void IF3D_ObstacleOperator::computeDiagnostics(const int stepID, const double time, const double* Uinf, const double lambda)
{
    BlockInfo * ary = &vInfo.front();
    const int N = vInfo.size();
    const double h = vInfo[0].h_gridpoint;
    double _area(0.0), _forcex(0.0), _forcey(0.0), _forcez(0.0), _torquex(0.0), _torquey(0.0), _torquez(0.0);
    double garea(0.0), gforcex(0.0), gforcey(0.0), gforcez(0.0), gtorquex(0.0), gtorquey(0.0), gtorquez(0.0);

    #pragma omp parallel for schedule(static) reduction(+:_area,_forcex,_forcey,_forcez,_torquex,_torquey,_torquez)
    for(int i=0; i<vInfo.size(); i++)
    {
        std::map<int,ObstacleBlock*>::const_iterator pos = obstacleBlocks.find(i);
        if(pos == obstacleBlocks.end()) continue;
        
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
            const double Xs = pos->second->chi[iy][ix];
            if (Xs == 0) continue;
            
            double p[2];
            info.pos(p, ix, iy, iz);
            p[0]-=position[0];
            p[1]-=position[1];
            p[2]-=position[2];

		    const double object_UR[3] = {
		    		angVel[1]*p[2]-angVel[2]*p[1],
					angVel[2]*p[0]-angVel[0]*p[2],
					angVel[0]*p[1]-angVel[1]*p[0]
		    };
            
            const double object_UDEF[2] = {
                pos->second->udef[iz][iy][ix][0],
                pos->second->udef[iz][iy][ix][1],
                pos->second->udef[iz][iy][ix][2]
            };
            
            const double U[3] = {
                b(ix,iy,iz).u + Uinf[0] - (transVel[0]+object_UR[0]+object_UDEF[0]),
                b(ix,iy,iz).v + Uinf[1] - (transVel[1]+object_UR[1]+object_UDEF[1]),
                b(ix,iy,iz).w + Uinf[2] - (transVel[2]+object_UR[2]+object_UDEF[2])
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
	MPI::COMM_WORLD.Allreduce(&_forcex, &gforcex, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&_forcey, &gforcey, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&_forcez, &gforcez, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&_torquex, &gtorquex, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&_torquey, &gtorquey, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&_torquez, &gtorquez, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&_area, &garea, 1, MPI::DOUBLE, MPI::SUM);

    const Real dA = std::pow(h, 3);
    mass     = garea  *dA;
    force[0] = gforcex*dA*lambda;
    force[1] = gforcey*dA*lambda;
    force[2] = gforcez*dA*lambda;
    torque[0] = gtorquex*dA*lambda;
    torque[1] = gtorquey*dA*lambda;
    torque[2] = gtorquez*dA*lambda;
    
    _writeDiagForcesToFile(stepID, time);
}

void IF3D_ObstacleOperator::computeVelocities(const double* Uinf)
{
    BlockInfo * ary = &vInfo.front();
    const int N = vInfo.size();
    double CM[3];
    this->getCenterOfMass(CM);
    const double dv = vInfo[0].h_gridpoint*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint;

    //local variables
    double V(0.0), J0(0.0), J1(0.0), J2(0.0), J3(0.0), J4(0.0), J5(0.0);
    double lm0(0.0), lm1(0.0), lm2(0.0); //linear momenta
    double am0(0.0), am1(0.0), am2(0.0); //angular momenta

    //global allReduce
    double gV(0.0), gJ0(0.0), gJ1(0.0), gJ2(0.0), gJ3(0.0), gJ4(0.0), gJ5(0.0);
    double glm0(0.0), glm1(0.0), glm2(0.0); //linear momenta
    double gam0(0.0), gam1(0.0), gam2(0.0); //angular momenta

#pragma omp parallel for schedule(static) reduction(+:V,lm0,lm1,lm2,J0,J1,J2,J3,J4,J5,am0,am1,am2)
    for(int i=0; i<vInfo.size(); i++) {
        std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(i);
        if(pos == obstacleBlocks.end()) continue;

        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
            const double Xs = pos->second->chi[iz][iy][ix];
            if (Xs == 0) continue;
            double p[3];
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

	MPI::COMM_WORLD.Allreduce(&lm0, &glm0, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&lm1, &glm1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&lm2, &glm2, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&am0, &gam0, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&am1, &gam1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&am2, &gam2, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J0, &gJ0, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J1, &gJ1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J2, &gJ2, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J3, &gJ3, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J4, &gJ4, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J5, &gJ5, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&V, &gV, 1, MPI::DOUBLE, MPI::SUM);

    transVel[0] = glm0/gV + Uinf[0];
    transVel[1] = glm1/gV + Uinf[1];
    transVel[2] = glm2/gV + Uinf[2];
    const Real _J[6] = {gJ0, gJ1, gJ2, gJ3, gJ4, gJ5};
    _finalizeAngVel(angVel, _J, gam0, gam1, gam2);
    volume = gV * dv;
    J[0] = gJ0 * dv; J[1] = gJ1 * dv; J[2] = gJ2 * dv;
    J[3] = gJ3 * dv; J[4] = gJ4 * dv; J[5] = gJ5 * dv;

    if(bFixToPlanar) {
        transVel[2] = 0.0;
        angVel[0] = angVel[1] = 0.0;
        angVel[2] = gam2/gJ2;
    }
}

void IF3D_ObstacleOperator::_finalizeAngVel(Real& AV[3], const Real& J[6], const Real& gam0, const Real& gam1, const Real& gam2)
{
	// try QR factorization to avoid dealing with determinant
	const double u1[3] = {J[0], J[3], J[4]};
	const double magu1sq = u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2];
	// subtract projection of a2 onto u1
	const double proj1 = u1[0]*J[3] + u1[1]*J[1] + u1[2]*J[5];
	const double u2[3] = {
			J[3] - proj1*u1[0]/magu1sq,
			J[1] - proj1*u1[1]/magu1sq,
			J[5] - proj1*u1[2]/magu1sq
	};
	const double magu2sq = u2[0]*u2[0] + u2[1]*u2[1] + u2[2]*u2[2];
	// subtract projection of a3 onto u1
	const double proj2 = u1[0]*J[4] + u1[1]*J[5] + u1[2]*J[2];
	const double u3_tmp[3] = {
			J[4] - proj2*u1[0]/magu1sq,
			J[5] - proj2*u1[1]/magu1sq,
			J[2] - proj2*u1[2]/magu1sq
	};
	// subtract projection of u3_tmp onto u2
	const double proj3 = u2[0]*u3_tmp[0] + u2[1]*u3_tmp[1] + u2[2]*u3_tmp[2];
	const double u3[3] = {
			u3_tmp[0] - proj3*u2[0]/magu2sq,
			u3_tmp[1] - proj3*u2[1]/magu2sq,
			u3_tmp[2] - proj3*u2[2]/magu2sq
	};
	const double magu3sq = u3[0]*u3[0] + u3[1]*u3[1] + u3[2]*u3[2];
	const double magu1 = std::sqrt(magu1sq);
	const double magu2 = std::sqrt(magu2sq);
	const double magu3 = std::sqrt(magu3sq);
	const double Q[3][3] = {
			{u1[0]/magu1, u2[0]/magu2, u3[0]/magu3},
			{u1[1]/magu1, u2[1]/magu2, u3[1]/magu3},
			{u1[2]/magu1, u2[2]/magu2, u3[2]/magu3}
	};
	// find out if Q is orthogonal
	const double R[3][3] = {
			{Q[0][0]*J[0] + Q[1][0]*J[3] + Q[2][0]*J[4], Q[0][0]*J[3] + Q[1][0]*J[1] + Q[2][0]*J[5], Q[0][0]*J[4] + Q[1][0]*J[5] + Q[2][0]*J[2]},
			{Q[0][1]*J[0] + Q[1][1]*J[3] + Q[2][1]*J[4], Q[0][1]*J[3] + Q[1][1]*J[1] + Q[2][1]*J[5], Q[0][1]*J[4] + Q[1][1]*J[5] + Q[2][1]*J[2]},
			{Q[0][2]*J[0] + Q[1][2]*J[3] + Q[2][2]*J[4], Q[0][2]*J[3] + Q[1][2]*J[1] + Q[2][2]*J[5], Q[0][2]*J[4] + Q[1][2]*J[5] + Q[2][2]*J[2]}
	};
	// d = Q^T b
	const double d[3] = {
			Q[0][0]*gam0 + Q[1][0]*gam1 + Q[2][0]*gam2,
			Q[0][1]*gam0 + Q[1][1]*gam1 + Q[2][1]*gam2,
			Q[0][2]*gam0 + Q[1][2]*gam1 + Q[2][2]*gam2,
	};
	// bwd subtitution: R x = d
	AV[2] = d[2]/R[2][2];
	AV[1] = (d[1] - R[1][2]*angVel[2])/R[1][1];
	AV[0] = (d[0] - R[0][1]*angVel[1] - R[0][2]*angVel[2])/R[0][0];
}

void IF3D_ObstacleOperator::computeForces(const int stepID, const double time, const double* Uinf, const double NU, const bool bDump)
{ //TODO: make sure this works for surfData not ordered by block index bt still "blocked" (it should)
    BlockInfo * ary = &vInfo.front();
    const int N = vInfo.size();
    Real vel_unit[3] = {0.0, 0.0, 0.0};
    const Real velx_tot = transVel[0]-Uinf[0];
    const Real vely_tot = transVel[1]-Uinf[1];
    const Real velz_tot = transVel[2]-Uinf[2];
    const Real vel_norm = std::sqrt(velx_tot*velx_tot + vely_tot*vely_tot + velz_tot*velz_tot);
    if (vel_norm>1e-9) {
        vel_unit[0] = (transVel[0]-Uinf[0])/vel_norm;
        vel_unit[1] = (transVel[1]-Uinf[1])/vel_norm;
        vel_unit[2] = (transVel[2]-Uinf[2])/vel_norm;
    }
    
    const int stencil_start[3] = {-1,-1, 0};
    const int stencil_end[3]   = { 2, 2, 1};
    
    //surfData is processed serially, so points are ordered by block
    vector<int> usefulIDs; //which blocks are occupied by surface
    vector<int> firstInfo; //which entries in surfData correspond to each
    for(int i=0; i<surfData.Ndata; i++) {
        bool unique(true); //if Ive already seen that ID then skip
        for(int k=0; k<usefulIDs.size(); k++)
            if (surfData.Set[i]->blockID == usefulIDs[k])
            { unique = false; break; }
        
        if (unique) {
            usefulIDs.push_back(surfData.Set[i]->blockID);
            firstInfo.push_back(i);
        }
    }
    firstInfo.push_back(surfData.Ndata);
    
    double _totChi(0.0),_totFx(0.0),_totFy(0.0),_totFz(0.0),_totFxP(0.0),_totFyP(0.0),_totFzP(0.0),_totFxV(0.0),_totFyV(0.0),_totFzV(0.0),_drag(0.0),_thrust(0.0),_Pout(0.0),_PoutBnd(0.0),_defPower(0.0),_defPowerBnd(0.0);
    #pragma omp parallel
    {
        Lab lab;
        lab.prepare(*grid, stencil_start, stencil_end, true);
        
        #pragma omp for schedule(static) reduction(+:_totChi,_totFx,_totFy,_totFz,_totFxP,_totFyP,_totFzP,_totFxV,_totFyV,_totFyV,_drag,_thrust,_Pout,_PoutBnd,_defPower,_defPowerBnd)
        for (int j=0; j<usefulIDs.size(); j++) {
            const int k = usefulIDs[j];
            lab.load(ary[k], 0);
            BlockInfo info = vInfo[k];
            
            const double _h3 = info.h_gridpoint*info.h_gridpoint*info.h_gridpoint;
            
            //????????? TO THINK ABOUT:
            const double _1oH = NU / info.h_gridpoint; // 2 nu / 2 h
            
            for(int i=firstInfo[j]; i<firstInfo[j+1]; i++)
            {
                double p[3];
                const int ix = surfData.Set[i]->ix;
                const int iy = surfData.Set[i]->iy;
                const int iz = surfData.Set[i]->iz;
                const auto tempIt = obstacleBlocks.find(k);
                assert(tempIt != obstacleBlocks.end());
                info.pos(p, ix, iy, iz);
                
                const double D11 =    _1oH*(lab(ix+1,iy,iz).u - lab(ix-1,iy,iz).u);
                const double D22 =    _1oH*(lab(ix,iy+1,iz).v - lab(ix,iy-1,iz).v);
                const double D33 =    _1oH*(lab(ix,iy,iz+1).w - lab(ix,iy,iz-1).w);

                const double D12 = .5*_1oH*(lab(ix,iy+1,iz).u - lab(ix,iy-1,iz).u
                                           +lab(ix+1,iy,iz).v - lab(ix-1,iy,iz).v);
                const double D13 = .5*_1oH*(lab(ix,iy,iz+1).u - lab(ix,iy,iz-1).u
                                           +lab(ix+1,iy,iz).w - lab(ix-1,iy,iz).w);
                const double D23 = .5*_1oH*(lab(ix,iy+1,iz).w - lab(ix,iy-1,iz).w
                                           +lab(ix,iy,iz+1).v - lab(ix,iy,iz-1).v);

                const double normX = surfData.Set[i]->dchidx * _h2;
                const double normY = surfData.Set[i]->dchidy * _h2;
                const double normZ = surfData.Set[i]->dchidz * _h2;
                const double fXV = D11 * normX + D12 * normY + D13 * normY;
                const double fYV = D12 * normX + D22 * normY + D23 * normY;
                const double fZV = D13 * normX + D23 * normY + D33 * normY;
                const double fXP = -lab(ix,iy).p * normX;
                const double fYP = -lab(ix,iy).p * normY;
                const double fZP = -lab(ix,iy).p * normZ;
                const double fXT = fXV+fXP;
                const double fYT = fYV+fYP;
                const double fZT = fZV+fZP;
                
                surfData.P[i]  = lab(ix,iy).p;
                surfData.fX[i] = fXT;  surfData.fY[i] = fYT;  surfData.fZ[i] = fZT;
                surfData.fxP[i] = fXP; surfData.fyP[i] = fYP; surfData.fzP[i] = fZP;
                surfData.fxV[i] = fXV; surfData.fyV[i] = fYV; surfData.fzV[i] = fZV;
                surfData.pX[i] = p[0]; surfData.pY[i] = p[1]; surfData.pZ[i] = p[2];
                
                _totChi += surfData.Set[i]->delta * _h3;
                _totFxP += fXP; _totFyP += fYP; _totFzP += fZP;
                _totFxV += fXV; _totFyV += fYV; _totFzV += fZV;
                _totFx  += fXT; _totFy  += fYT; _totFz  += fZT;
                
                
                const double forcePar = fXT*vel_unit[0] + fYT*vel_unit[1] + fZT*vel_unit[2];
                _thrust += .5*(forcePar + std::abs(forcePar));
                _drag   -= .5*(forcePar - std::abs(forcePar));
                
                surfData.vxDef[i] = tempIt->second->udef[iz][iy][ix][0];
                surfData.vyDef[i] = tempIt->second->udef[iz][iy][ix][1];
                surfData.vzDef[i] = tempIt->second->udef[iz][iy][ix][2];
                surfData.vx[i] = lab(ix,iy,iz).u + Uinf[0];
                surfData.vy[i] = lab(ix,iy,iz).v + Uinf[1];
                surfData.vz[i] = lab(ix,iy,iz).w + Uinf[2];
                
                const double powOut = fXT*surfData.vx[i] + fYT*surfData.vy[i] + fZT*surfData.vz[i];
                _Pout   += powOut;
                _PoutBnd+= min(0., powOut);
                
                const double powDef = fXT*surfData.vxDef[i] + fYT*surfData.vyDef[i] + fZT*surfData.vzDef[i];
                _defPower   += powDef;
                _defPowerBnd+= min(0., powDef);
            }
        }
    }

    totChi=_totChi;
    totFx=_totFx; totFy=_totFy; totFz=_totFz;
    drag=_drag; thrust=_thrust;
    Pout=_Pout; PoutBnd=_PoutBnd;
    defPower=_defPower; defPowerBnd=_defPowerBnd;
    Pthrust    = thrust*vel_norm;
    Pdrag      =   drag*vel_norm;
    EffPDef    = Pthrust/(Pthrust-min(defPower,0.));
    EffPDefBnd = Pthrust/(Pthrust-    defPowerBnd);
    
    if (bDump)  surfData.print(obstacleID, stepID);
    if (bDump && pX not_eq nullptr && pY not_eq nullptr && pZ not_eq nullptr && Npts not_eq 0) {
    	surfData.sort(pX,pY,pZ,Npts);
    	surfData.printSorted(obstacleID, stepID);
    }
    
    {
        ofstream filedrag;
        filedrag.open(("forceValues_"+std::to_string(obstacleID)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<_totFxP<<" "<<_totFyP<<" "<<_totFzP<<" "
        		 			   <<_totFxV<<" "<<_totFyV<<" "<<_totFzV<<" "
						   << totFx <<" "<< totFy <<" "<< totFz <<" "
						   << drag  <<" "<< thrust<<" "<< totChi<<endl;
        filedrag.close();
    }
    {
        ofstream filedrag;
        filedrag.open(("powerValues_"+std::to_string(obstacleID)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<Pthrust<<" "<<Pdrag<<" "<<PoutBnd<<" "<<Pout<<" "<<defPowerBnd<<" "<<defPower<<" "<<EffPDefBnd<<" "<<EffPDef<<endl;
        filedrag.close();
    }
}

void IF3D_ObstacleOperator::update(const int step_id, const Real t, const Real dt)
{
    position[0] += dt*transVel[0];
    position[1] += dt*transVel[1];
    position[2] += dt*transVel[2];

    const double dqdt[4] = {
        0.5*(-angVel[0]*quaternion[1]-angVel[1]*quaternion[2]-angVel[2]*quaternion[3]),
        0.5*( angVel[0]*quaternion[0]+angVel[1]*quaternion[3]-angVel[2]*quaternion[2]),
        0.5*(-angVel[0]*quaternion[3]+angVel[1]*quaternion[0]+angVel[2]*quaternion[1]),
        0.5*( angVel[0]*quaternion[2]-angVel[1]*quaternion[1]+angVel[2]*quaternion[0])
    };

    // normality preserving advection (Simulation of colliding constrained rigid bodies - Kleppmann 2007 Cambridge University, p51)
    // move the correct distance on the quaternion unit ball surface, end up with normalized quaternion
    const double deltaq[4] = {
        dqdt[0]*dt,
        dqdt[1]*dt,
        dqdt[2]*dt,
        dqdt[3]*dt
    };

    const double deltaq_length = std::sqrt(deltaq[0]*deltaq[0]+deltaq[1]*deltaq[1]+deltaq[2]*deltaq[2]+deltaq[3]*deltaq[3]);

    if(deltaq_length>std::numeric_limits<double>::epsilon()) {
        const double tanfac = std::tan(deltaq_length)/deltaq_length;
        const double num[4] = {
            quaternion[0]+tanfac*deltaq[0],
            quaternion[1]+tanfac*deltaq[1],
            quaternion[2]+tanfac*deltaq[2],
            quaternion[3]+tanfac*deltaq[3],
        };

        const double invDenum = 1./(std::sqrt(num[0]*num[0]+num[1]*num[1]+num[2]*num[2]+num[3]*num[3]));

        quaternion[0] = num[0]*invDenum;
        quaternion[1] = num[1]*invDenum;
        quaternion[2] = num[2]*invDenum;
        quaternion[3] = num[3]*invDenum;
    }

#ifndef NDEBUG
    std::cout << "POSITION INFO AFTER UPDATE T, DT: " << t << " " << dt << std::endl;
    std::cout << "POS: " << position[0] << " " << position[1] << " " << position[2] << std::endl;
    std::cout << "TVL: " << transVel[0] << " " << transVel[1] << " " << transVel[2] << std::endl;
    std::cout << "QUT: " << quaternion[0] << " " << quaternion[1] << " " << quaternion[2] << " " << quaternion[3] << std::endl;
    std::cout << "AVL: " << angVel[0] << " " << angVel[1] << " " << angVel[2] << std::endl;
    
    const double q_length=std::sqrt(quaternion[0]*quaternion[0]+quaternion[1]*quaternion[1]+quaternion[2]*quaternion[2]+quaternion[3]*quaternion[3]);
    assert(std::abs(q_length-1.0) < 5*std::numeric_limits<double>::epsilon());
#endif
    _writeComputedVelToFile(step_id, t);
}

void IF3D_ObstacleOperator::characteristic_function()
{
#pragma omp parallel
	{
		BlockInfo * ary = &vInfo.front();
#pragma omp for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(i);

			if(pos != obstacleBlocks.end()) {
				FluidBlock& b = *(FluidBlock*)ary[i].ptrBlock;
				for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
					b(ix,iy,iz).chi = std::max(pos->second->chi[iz][iy][ix], b(ix,iy,iz).chi);
			}
		}
	}
}

std::vector<int> IF3D_ObstacleOperator::intersectingBlockIDs(const int buffer)
{
	assert(buffer <= 2); // only works for 2: if different definition of deformation blocks, implement your own
	std::vector<int> retval;
	const int N = vInfo.size();

	for(int i=0; i<N; i++) {
		std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(i);
		if(pos != obstacleBlocks.end())
			retval.push_back(i);
	}
	return retval;
}

void IF3D_ObstacleOperator::getTranslationVelocity(Real UT[3]) const
{
    UT[0]=transVel[0];
    UT[1]=transVel[1];
    UT[3]=transVel[3];
}

void IF3D_ObstacleOperator::setTranslationVelocity(Real UT[3])
{
    transVel[0]=UT[0];
    transVel[1]=UT[1];
    transVel[2]=UT[2];
}

void IF3D_ObstacleOperator::getAngularVelocity(Real & W[3]) const
{
    W[0]=angVel[0];
    W[1]=angVel[1];
    W[3]=angVel[3];
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

double IF3D_ObstacleOperator::getForceX() const
{
    return force[0];
}
double IF3D_ObstacleOperator::getForceY() const
{
    return force[1];
}

void IF3D_ObstacleOperator::save(const int step_id, const Real t, std::string filename)
{
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
    
    if(filename==std::string())
        restartstream.open("restart_IF3D_MovingBody.txt");
    else
        restartstream.open(filename+".txt");
    
    Real restart_time;
    restartstream >> restart_time;
    assert(std::abs(restart_time-t) < std::numeric_limits<Real>::epsilon());
    
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

