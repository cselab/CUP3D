//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_CarlingFishOperator.h"
#include "IF3D_ObstacleLibrary.h"
#include "GenericOperator.h"

Fish::CarlingFishMidlineData::CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext)
: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext)
{
#ifdef BBURST
	ifstream reader("burst_coast_carling_params.txt");
	if (reader.is_open()) {
		reader >> t0;
		reader >> t1;
		reader >> t2;
		reader >> t3;
		reader.close();
	} else {
		cout << "Could not open params.txt" << endl;
		abort();
	}
#endif
}

void Fish::CarlingFishMidlineData::computeMidline(const Real time)
{
	_computeMidlineCoordinates(time);
	_computeMidlineVelocities(time);
	_computeMidlineNormals();
#ifndef NDEBUG
// we dump the profile
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank!=0) return;
	FILE * f = fopen("fish_profile","w");
	for(int i=0;i<Nm;++i)
		fprintf(f,"%d %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",
				i,rS[i],rX[i],rY[i],norX[i],norY[i],vX[i],vY[i],vNorX[i],vNorY[i],width[i],height[i]);
	fclose(f);
	printf("Dumped midline\n");
#endif
}

void Fish::CarlingFishMidlineData::_computeMidlineCoordinates(const Real time)
{
	const Real rampFac = rampFactorSine(time, Tperiod);
	rX[0] = 0.0;
	rY[0] =   rampFac*midline(rS[0], time, length, Tperiod, phaseShift);

	for(int i=1;i<Nm;++i) {
		rY[i]=rampFac*midline(rS[i], time, length, Tperiod, phaseShift);
		const Real dy = rY[i]-rY[i-1];
		const Real ds = rS[i] - rS[i-1];
		const Real dx = std::sqrt(ds*ds-dy*dy);
		rX[i] = rX[i-1] + dx;
	}
}

void Fish::CarlingFishMidlineData::_computeMidlineVelocities(const Real time)
{
	const Real rampFac =    rampFactorSine(time, Tperiod);
	const Real rampFacVel = rampFactorVelSine(time, Tperiod);

	vX[0] = 0.0; //rX[0] is constant
	vY[0] = rampFac*midlineVel(rS[0],time,length,Tperiod, phaseShift) +
			rampFacVel*midline(rS[0],time,length,Tperiod, phaseShift);

	for(int i=1;i<Nm;++i) {
		vY[i]=rampFac*midlineVel(rS[i],time,length,Tperiod, phaseShift) +
				rampFacVel*midline(rS[i],time,length,Tperiod, phaseShift);
		const Real dy = rY[i]-rY[i-1];
		const Real dx = rX[i]-rX[i-1];
		const Real dVy = vY[i]-vY[i-1];
		assert(dx>0); // has to be, otherwise y(s) is multiple valued for a given s
		vX[i] = vX[i-1] - dy/dx * dVy; // use ds^2 = dx^2 + dy^2 --> ddx = -dy/dx*ddy
	}
}

IF3D_CarlingFishOperator::IF3D_CarlingFishOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_ObstacleOperator(grid, parser), theta_internal(0.0), angvel_internal(0.0), sim_time(0.0), sim_dt(0.0), adjTh(adjTh), myFish(nullptr)
{
	_parseArguments(parser);
	const Real target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
	const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
	const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG)+1;
	printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
	fflush(0);
	// multiple of NPPSEG: TODO why?
	myFish = new Fish::CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension);
}

IF3D_CarlingFishOperator::~IF3D_CarlingFishOperator()
{
	if(myFish not_eq nullptr) delete myFish;
}

void IF3D_CarlingFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_FishOperator::_parseArguments(parser);
}
