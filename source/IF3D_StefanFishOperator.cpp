//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_StefanFishOperator.h"
#include "IF3D_FishLibrary.h"

void IF3D_StefanFishOperator::save(const int step_id, const Real t, std::string filename)
{
	//assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
	std::ofstream savestream;
	savestream.setf(std::ios::scientific);
	savestream.precision(std::numeric_limits<Real>::digits10 + 1);
	savestream.open(filename + ".txt");

	savestream<<t<<"\t"<<sim_dt<<std::endl;
	savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
	savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
	savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
	savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
	savestream<<theta_internal<<"\t"<<angvel_internal<<"\t"<<adjTh<<std::endl;
	savestream<<_2Dangle<<"\t"<<old_curv<<"\t"<<new_curv<<std::endl;
	savestream.close();

	myFish->curvScheduler.save(filename+"_curv");
	myFish->baseScheduler.save(filename+"_base");
	myFish->adjustScheduler.save(filename+"_adj");
	sr.save(step_id, filename);
}

void IF3D_StefanFishOperator::restart(const Real t, string filename)
{
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal >> adjTh;
    restartstream >> _2Dangle >> old_curv >> new_curv;
    restartstream.close();
		#if 1
		sr.restart(filename);
		myFish->curvScheduler.restart(filename+"_curv");
		myFish->baseScheduler.restart(filename+"_base");
		myFish->adjustScheduler.restart(filename+"_adj");
		#else
		if (bCorrectTrajectory)
		{
			Real velx_tot = 4.67373623945047042549e-02 - transVel[0];
			Real vely_tot = 4.67885955946669776506e-02 - transVel[1];
			Real AngDiff  = std::atan2(vely_tot,velx_tot);
			const Real B = (AngDiff*angVel[2]>0) ? 0.25/M_PI : 0;
			const Real PID = (.5*adjTh +B*AngDiff*fabs(angVel[2]) );
			myFish->adjustScheduler.t0 = sim_time;
			myFish->adjustScheduler.t1 = sim_time + 0.00026;
			for (int i=0; i<6; ++i) {
				myFish->adjustScheduler.parameters_t0[i] = PID/M_PI;
				myFish->adjustScheduler.parameters_t1[i] = PID/M_PI;
			}
			sr.t_next_comm = 1e5;
		}
		else
		{
			std::array<Real,7> tmp_curv = {0, -.592913, 0, .00590721, 0, -.131738, 0};
			myFish->baseScheduler.parameters_t0 = tmp_curv;
			myFish->baseScheduler.t0 = 4.5;
			sr.t_next_comm = 5.;
		}
		#endif
		if(!rank)
		{
		std::cout<<"RESTARTED FISH: "<<std::endl;
		std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
		std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
		std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
		std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
		std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
		std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
		std::cout<<"2D angle: \t"<<_2Dangle<<std::endl;
		}
}

IF3D_StefanFishOperator::IF3D_StefanFishOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_FishOperator(grid, parser)
{
	_parseArguments(parser);
	const int Nextension = NEXTDX*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
	const Real target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
	const Real dx_extension = (1./NEXTDX)*vInfo[0].h_gridpoint;
	const int Nm = (Nextension+1)*(int)std::ceil(target_Nm/(Nextension+1)) + 1;

	printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
	fflush(0);
	myFish = new CurvatureDefinedFishData(Nm, length, Tperiod, phaseShift, dx_extension);

  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, transVel[0], transVel[1], angVel[2]);
}

void IF3D_StefanFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_FishOperator::_parseArguments(parser);
  sr.t_next_comm = Tstartlearn - 1/2.; //i want to reset time-averaged quantities before first actual comm
  sr.GoalDX = GoalDX;
  sr.thExp = _2Dangle;
	/*
    randomActions = parser("-randomActions").asBool(false);
    if (randomActions) printf("Fish doing random turns\n");
  */
	useLoadedActions = parser("-useLoadedActions").asBool(false);
	if (useLoadedActions) {
		printf("Trying to load actionsi %d.\n",nActions);
		fflush(0);
		Real dummy_time;
		vector<Real> action(nActions);
		ifstream in("orders_1.txt");
		std::string line;
		if(in.good()) {
			while (getline(in, line)) {
				istringstream line_in(line);
				if(nActions==2)
					line_in >> dummy_time >> action[0] >> action[1];
				else
					line_in >> dummy_time >> action[0];
				//i want to do pop back later:
				loadedActions.insert(loadedActions.begin(),action);
			}
			in.close();
		} else {
			printf("Could not load actions from file orders_1.txt\n");
			MPI_Abort(grid->getCartComm(), MPI_ERR_OTHER);
		}
	}
}

void IF3D_StefanFishOperator::execute(Communicator * comm, const int iAgent,
																			const Real time, const int iLabel)
{
    if (time < Tstartlearn) {
        sr.resetAverage();
        sr.t_next_comm = Tstartlearn;
        return;
    }

    if (useLoadedActions) {
      vector<Real> actions(nActions);
      if (loadedActions.size()>1) {
          actions = loadedActions.back();
          loadedActions.pop_back();
      } //else zero actions
			else
			 MPI_Abort(grid->getCartComm(), MPI_ERR_OTHER);

      myFish->execute(time, sr.t_next_comm, actions);

      old_curv = new_curv;
      new_curv = actions[0];
      //if(nActions==2) {
      //    new_Tp = actions[1];
      //    sr.t_next_comm += .5*myFish->l_Tp;
      //} else if (nActions==1) {
          sr.t_next_comm += .5*myFish->Tperiod;
      //}
      sr.resetAverage();
    } else if (not bInteractive) {
      sr.t_next_comm=1e6;
      return;
    } else if (comm not_eq nullptr) {
      const Real relT= fmod(time,1.); //1 is Tperiod
      #ifdef _NOVISION_
        const int nStates = (nActions==1) ? 20+ 8*20 : 25+  8*20;
      #else
        const int nStates = (nActions==1) ? 20+10*20 : 25+ 10*20;
      #endif
      vector<Real> state(nStates), actions(nActions);

      int k = 0;
      state[k++] = sr.Xpov - GoalDX;
      state[k++] = sr.Ypov;
      state[k++] = sr.RelAng;
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
          state[k++] = sr.VX;
          state[k++] = sr.VY;
          state[k++] = sr.AV;
      }
      */
      state[k++] = sr.Dist;
      state[k++] = sr.Quad;
      state[k++] = sr.VxAvg;
      state[k++] = sr.VyAvg;
      state[k++] = sr.AvAvg;
      state[k++] = sr.Pout;
      state[k++] = sr.defPower;
      state[k++] = sr.EffPDef;
      state[k++] = sr.PoutBnd;
      state[k++] = sr.defPowerBnd;
      state[k++] = sr.EffPDefBnd;
      state[k++] = sr.Pthrust;
      state[k++] = sr.Pdrag;
      state[k++] = sr.ToD;
      /*
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelNAbove[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelTAbove[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelNBelow[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelTBelow[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FPAbove[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FVAbove[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FPBelow[j];
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FVBelow[j];
      #ifndef _NOVISION_
      for (int j=0; j<2*NpLatLine; j++) state[k++] = sr.raySight[j];
      #endif
      */
      //fflush(0);
      const Real reward = (sr.info==2) ? -10 : sr.EffPDefBnd;
      comm->sendState(iLabel, sr.info, state, reward); //TODO
      if (sr.info==2) return;
      sr.info = 0;

      comm->recvAction(actions);
      myFish->execute(time, sr.t_next_comm, actions);

      old_curv = new_curv;
      new_curv = actions[0];
      //if(nActions==2) {
      //    new_Tp = actions[1];
      //    sr.t_next_comm += .5*myFish->l_Tp;
      //} else if (nActions==1) {
          sr.t_next_comm += .5*myFish->Tperiod;
      //}
      if(!rank) {
        printf("Next action of agent %d at time %g\n", iAgent, sr.t_next_comm);
        ofstream filedrag;
        filedrag.open(("orders_"+to_string(iAgent)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<new_curv;
        //if(nActions==2) filedrag<<" "<<new_Tp;
        filedrag<<endl;
        filedrag.close();
      }

      sr.resetAverage();
      fflush(0);
    }
}
