//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_StefanFishOperator.h"
#include "IF3D_FishLibrary.h"

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
  sr->t_next_comm = Tstartlearn - 1/2.; //i want to reset time-averaged quantities before first actual comm
  sr->bForgiving = bForgiving;
  sr->GoalDX = GoalDX;
  sr->thExp = angle;
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

        const int nActions = 2;
        const Real relT= fmod(time,1.); //1 is Tperiod
#ifdef _NOVISION_
        const int nStates = (nActions==1) ? 20+ 8*20 : 25+  8*20;
#else
        const int nStates = (nActions==1) ? 20+10*20 : 25+ 10*20;
#endif
        vector<Real> state(nStates), actions(nActions);

        int k(0);
        state[k++] = sr->Xrel - GoalDX;
        state[k++] = sr->Yrel;
        state[k++] = sr->RelAng;
        state[k++] = relT;
        state[k++] = new_curv;
        state[k++] = old_curv;
/*
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
*/
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
        /*
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
        */
        const Real reward = (sr->info==2) ? -10 : sr->EffPDefBnd;
        comm->sendState(iAgent-1, sr->info, state, reward); //TODO
        fflush(0);
        if (sr->info==2) return;

        sr->info = 0;

        comm->recvAction(actions);
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        //if(nActions==2) {
        //    new_Tp = actions[1];
        //    sr->t_next_comm += .5*myFish->l_Tp;
        //} else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        //}

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
}
