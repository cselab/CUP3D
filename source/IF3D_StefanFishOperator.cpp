//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_StefanFishOperator.h"
#include "IF3D_FishLibrary.h"

CurvatureDefinedFishData::CurvatureDefinedFishData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext)
: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), l_Tp(Tperiod), timeshift(0), time0(0),
  rK(_alloc(Nm)),vK(_alloc(Nm)),rC(_alloc(Nm)),vC(_alloc(Nm)),rA(_alloc(Nm)),vA(_alloc(Nm)),rB(_alloc(Nm)),vB(_alloc(Nm))
{    	}

void CurvatureDefinedFishData::_correctTrajectory(const Real dtheta, const Real time, const Real dt) override
{
	std::array<Real,6> tmp_curv = std::array<Real,6>();
	for (int i=0; i<tmp_curv.size(); ++i) {tmp_curv[i] = dtheta/M_PI;}
	adjustScheduler.transition(time,time,time+2*dt,tmp_curv, true);
}

void CurvatureDefinedFishData::execute(const Real time, const Real l_tnext, const vector<Real>& input) override
{
	if (input.size()>1) {
		baseScheduler.Turn(input[0], l_tnext);
		//first, shift time to  previous turn node
		timeshift += (l_tnext-time0)/l_Tp;
		time0 = l_tnext;
		l_Tp = Tperiod*(1.+input[1]);
	} else if (input.size()>0) {
		baseScheduler.Turn(input[0], l_tnext);
	}
}

CurvatureDefinedFishData::~CurvatureDefinedFishData()
{
	_dealloc(rK);
	_dealloc(vK);
	_dealloc(rC);
	_dealloc(vC);
	_dealloc(rB);
	_dealloc(vB);
	_dealloc(rA);
	_dealloc(vA);
}

void CurvatureDefinedFishData::computeMidline(const Real time)
{
	const Real _1oL = 1./length;
	const std::array<Real ,6> curvature_values = {
			0.82014*_1oL, 1.46515*_1oL, 2.57136*_1oL,
			3.75425*_1oL, 5.09147*_1oL, 5.70449*_1oL
	};
	const std::array<Real ,6> curvature_points = {
			0., .15*length, .4*length, .65*length, .9*length, length
	};
	const std::array<Real ,7> baseline_points = {
			1.00, 0.75, 0.50, 0.25, 0.00, -0.25, -0.50
	};
	const std::array<Real, 6> curvature_zeros = std::array<Real, 6>();
	curvScheduler.transition(time,0.0,Tperiod,curvature_zeros,curvature_values);

	// query the schedulers for current values
	curvScheduler.gimmeValues(time, curvature_points, Nm, rS, rC, vC);
	baseScheduler.gimmeValues(time, l_Tp, length, baseline_points, Nm, rS, rB, vB);
	adjustScheduler.gimmeValues(time, curvature_points, Nm, rS, rA, vA);

	// construct the curvature
	const Real _1oT = 1./l_Tp;
	for(unsigned int i=0; i<Nm; i++) {
		const Real darg = 2.*M_PI* _1oT;
		const Real arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL) + M_PI*phaseShift;
		rK[i] = rC[i]*(std::sin(arg) + rB[i] + rA[i]);
		vK[i] = vC[i]*(std::sin(arg) + rB[i] + rA[i]) + rC[i]*(std::cos(arg)*darg + vB[i] + vA[i]);
	}

#if 1==0
	{ // we dump the profile points
		FILE * f = fopen("stefan.dat","a");
		std::array<Real, 6> curv,base;
		curvScheduler.ParameterScheduler<6>::gimmeValues(time, curv);
		baseScheduler.ParameterScheduler<6>::gimmeValues(time, base);
		fprintf(f,"%9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e\n",
				time,curv[0],curv[1],curv[2],curv[3],curv[4],curv[5],base[0],base[1],base[2],base[3],base[4],base[5]);
		fclose(f);
	}
	{ // we dump the profile
		FILE * f = fopen("stefan_profile","w");
		for(int i=0;i<Nm;++i)
			fprintf(f,"%d %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e\n",
					i,rS[i],rX[i],rY[i],norX[i],norY[i],vX[i],vY[i],vNorX[i],vNorY[i],width[i]);
		fclose(f);
	}
#endif

	// solve frenet to compute midline parameters
	IF2D_Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
}

IF3D_StefanFishOperator::IF3D_StefanFishOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_FishOperator(grid, parser)
  {
	_parseArguments(parser);
	const Real target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
	const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
	const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG)+1;
	printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
	fflush(0);
	// multiple of NPPSEG: TODO why?
	myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension);
  }

void IF3D_StefanFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_FishOperator::_parseArguments(parser);
	/*
    randomActions = parser("-randomActions").asBool(false);
    if (randomActions) printf("Fish doing random turns\n");
    useLoadedActions = parser("-useLoadedActions").asBool(false);
    if (useLoadedActions) {
        Real dummy_time;
        vector<Real> action(nActions);
        ifstream in("orders_1.txt"); //FUCKING TODO NEED TO USE SOME POINTERS IN THIS SHIT
        std::string line;
        if(in.good()) {
            while (getline(in, line)) {
                istringstream line_in(line);
                line_in >> dummy_time;
                line_in >> action[0];
                if(nActions==2) line_in >> action[1];
                //i want to do pop back later:
                loadedActions.insert(loadedActions.begin(),action);
            }
        } else { printf("Could not load actions from file orders_1.txt\n"); abort(); }
        in.close();
    }
    */
}

void IF3D_StefanFishOperator::execute(Communicator * comm, const int iAgent, const Real time)
{
	/*
    if (time < Tstartlearn) {
        sr->resetAverage();
        sr->t_next_comm = Tstartlearn;

        //TMP: first rnd action is going to be taken after a while
        if (randomActions) sr->t_next_comm = Tstartlearn+6.;

        return;
    }

    if (not bInteractive) {
        if (not randomActions) { sr->t_next_comm=1e6; return; }
        //we might decide to pause turning, then just pick a pause counter
        if (nPauseActions-- > 0) {
            vector<Real> raction(1,0.);
            myFish->execute(time, sr->t_next_comm, raction);
            sr->t_next_comm += .5*myFish->Tperiod;
            printf("pausing from turning at time %f, will still pause %d turns.
            		Next turn at time %f\n",time,nPauseActions,sr->t_next_comm);
            return;
        }
        vector<Real> raction(1,signLastTurn);
        myFish->execute(time, sr->t_next_comm, raction);
        sr->t_next_comm += .5*myFish->Tperiod;
        printf("turning at time %f with modifier %f, turn counter is %d.
        		Next turn at time %f\n",time,signLastTurn,nTurnActions,sr->t_next_comm);

        if (++nTurnActions >= 4) {
            nPauseActions = 10.;
            nTurnActions = 0;
            signLastTurn *= -1.;
        }

    } else if (useLoadedActions) {

        vector<Real> actions(nActions);
        if (loadedActions.size()>1) {
            actions = loadedActions.back();
            loadedActions.pop_back();
        } //else zero actions
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }
        sr->resetAverage();

    } else {

        const Real relT= fmod(time,1.); //1 is Tperiod
#ifdef _NOVISION_
        const int nStates = (nActions==1) ? 20+ 8*NpLatLine : 25+  8*NpLatLine;
#else
        const int nStates = (nActions==1) ? 20+10*NpLatLine : 25+ 10*NpLatLine;
#endif
        vector<Real> state(nStates), actions(nActions);

        int k(0);
        state[k++] = sr->Xrel - GoalDX;
        state[k++] = sr->Yrel;
        state[k++] = sr->RelAng;
        state[k++] = relT;
        state[k++] = new_curv;
        state[k++] = old_curv;

        if(nActions==2) { //this is for backwards compatibility
            state[k++] = new_Tp;
                        //2.*M_PI*((time-time0)/l_Tp +timeshift -rS[i]/length) + M_PI*phaseShift
            Real Fshift = 2.*((-myFish->time0)/myFish->l_Tp +myFish->timeshift)+myFish->phaseShift;
            Fshift = fmod(Fshift,2.0);
            state[k++] = (Fshift<0) ? 2.+Fshift : Fshift;
            state[k++] = sr->VX;
            state[k++] = sr->VY;
            state[k++] = sr->AV;
        }

        state[k++] = sr->Dist;
        state[k++] = sr->Quad;
        state[k++] = sr->VxAvg;
        state[k++] = sr->VyAvg;
        state[k++] = sr->AvAvg;
        state[k++] = sr->Pout;
        state[k++] = sr->defPower;
        state[k++] = sr->EffPDef;
        state[k++] = sr->PoutBnd;
        state[k++] = sr->defPowerBnd;
        state[k++] = sr->EffPDefBnd;
        state[k++] = sr->Pthrust;
        state[k++] = sr->Pdrag;
        state[k++] = sr->ToD;

        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVBelow[j];
        #ifndef _NOVISION_
        for (int j=0; j<2*NpLatLine; j++) state[k++] = sr->raySight[j];
        #endif
        const Real reward = (sr->info==2) ? -10 : sr->EffPDefBnd;
        comm->sendState(iAgent-1, sr->info, state, reward); //TODO
        fflush(0);
        if (sr->info==2) return;

        sr->info = 0;

        comm->recvAction(actions);
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }

        #ifndef TRAINING
        ofstream filedrag;
        filedrag.open(("orders_"+to_string(iAgent)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<new_curv;
        if(nActions==2)
            filedrag<<" "<<new_Tp;
        filedrag<<endl;
        filedrag.close();
        #endif //TRAINING

        sr->resetAverage();
    }
    */
}
