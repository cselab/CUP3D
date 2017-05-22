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

	const double timeshift = myFish->timeshift;
	const double time0 = myFish->time0;
	const double l_Tp = myFish->l_Tp;

	savestream<<t<<"\t"<<sim_dt<<std::endl;
	savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
	savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
	savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
	savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
	savestream<<theta_internal<<"\t"<<angvel_internal<<"\t"<<adjTh<<std::endl;
	savestream<<timeshift<<"\t"<<time0<<"\t"<<l_Tp<<"\t"<<new_curv<<"\t"<<old_curv<<"\t"<<new_Tp<<std::endl;
	savestream<<_2Dangle<<"\t"<<old_curv<<"\t"<<new_curv<<std::endl;
	savestream.close();

	myFish->curvScheduler.save(filename+"_curv");
	myFish->baseScheduler.save(filename+"_base");
	myFish->adjustScheduler.save(filename+"_adj");
	sr.save(step_id, filename);
}

void IF3D_StefanFishOperator::restart(const Real t, string filename)
{
		double timeshift, time0, l_Tp;
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal >> adjTh;
    restartstream >> timeshift >> time0 >> l_Tp >> new_curv >> old_curv >> new_Tp;
    restartstream >> _2Dangle >> old_curv >> new_curv;
    restartstream.close();

		sr.restart(filename);
		myFish->curvScheduler.restart(filename+"_curv");
		myFish->baseScheduler.restart(filename+"_base");
		myFish->adjustScheduler.restart(filename+"_adj");
		myFish->timeshift = timeshift;
    myFish->time0 = time0;
    myFish->l_Tp = l_Tp;

		if(!rank)
		{
		std::cout<<"RESTARTED FISH: "<<std::endl;
		std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
		std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
		std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
		std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
		std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
		std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
    std::cout<<"TIMESHIFT: "<<timeshift<<" "<<time0<<" "<<l_Tp<<std::endl;
    std::cout<<"ACTIONS: "<<new_curv<<" "<<old_curv<<" "<<new_Tp<<std::endl;
		std::cout<<"2D angle: "<<_2Dangle<<std::endl;
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

	if(!rank) printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
	fflush(0);
	myFish = new CurvatureDefinedFishData(Nm, length, Tperiod, phaseShift, dx_extension);

  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, transVel[0], transVel[1], angVel[2]);
}

void IF3D_StefanFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_FishOperator::_parseArguments(parser);
	/*
    randomActions = parser("-randomActions").asBool(false);
    if (randomActions) printf("Fish doing random turns\n");
  */
	const bool randomStart = parser("-randomStart").asBool(false);
	if (randomStart) {
			double random_starts[4];
			if (!rank) {
				std::random_device rd;
	      std::mt19937 gen(rd());
	      std::uniform_real_distribution<double> dis(-1.,1.);
				random_starts[0]=dis(gen);
				random_starts[1]=dis(gen);
				random_starts[2]=dis(gen);
				random_starts[3]=dis(gen);

				#ifdef __ExploreHalfWake
				//now, let's explore only one half of the domain:
					random_starts[1] = std::fabs(random_starts[1]);
				#endif

	      for (int i = 1; i < size; i++)
          MPI_Send(random_starts, 4, MPI_DOUBLE, i, 35, grid->getCartComm());

	    } else
          MPI_Recv(random_starts, 4, MPI_DOUBLE, 0, 35, grid->getCartComm(), MPI_STATUS_IGNORE);
			printf("Rank %d (out of %d) using seeds %g %g %g %g\n", rank, size,
					random_starts[0],random_starts[1],random_starts[2],random_starts[3]);

      position[0] += .5*length*random_starts[0];
			if (!rank) printf("Assuming leader is in x=0.15\n");
			const Real dX = position[0]-0.15;
			//now adding a shift so that i do not over explore dy = 0;
			const Real shiftDy = random_starts[1]>0 ? (0.25*dX/2.5) : -(0.25*dX/2.5);
      position[1] += .25*length*random_starts[1] + shiftDy;
      absPos[0] = position[0];
      absPos[1] = position[1];
      _2Dangle = .1* M_PI *random_starts[2];

			quaternion[0] = std::cos(0.5*_2Dangle);
			quaternion[1] = 0;
			quaternion[2] = 0;
			quaternion[3] = std::sin(0.5*_2Dangle);
      if(nActions==2) phaseShift = random_starts[3];
    }

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
		const double invlscale  = 1./length;
		const double velscale   = Tperiod*invlscale; //all these are inverse
		const double forcescale = velscale*velscale*invlscale*invlscale; //rho*l^3*l/t^2
		const double powerscale = forcescale*velscale; //rho*l^2*l^2/t^3

		const Real eps = std::numeric_limits<Real>::epsilon();
		assert(std::fabs(std::cos(0.5*_2Dangle)-quaternion[0]) < 1e-6);
		assert(std::fabs(quaternion[1]) < eps);
		assert(std::fabs(quaternion[2]) < eps);
		assert(std::fabs(std::sin(0.5*_2Dangle)-quaternion[3]) < 1e-6);

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
      if(nActions>1) {
          new_Tp = actions[1];
          sr.t_next_comm += .5*myFish->l_Tp;
      } else if (nActions==1) {
          sr.t_next_comm += .5*myFish->Tperiod;
      }
      sr.resetAverage();
    } else if (not bInteractive) {
      sr.t_next_comm=1e6;
      return;
    } else if (comm not_eq nullptr) {
      const Real relT= std::fmod(time,1.); //1 is Tperiod of leader
      const int nStates = (nActions==1) ? 20+10*__NpLatLine : 25+10*__NpLatLine;
      vector<Real> state(nStates), actions(nActions);

      int k = 0;
      //state[k++] = sr.Xpov*invlscale - GoalDX;
			state[k++] = sr.Xpov*invlscale;
      state[k++] = sr.Ypov*invlscale;
      state[k++] = sr.RelAng;
      state[k++] = relT;
      state[k++] = new_curv;
      state[k++] = old_curv;

      if(nActions==2)
      {
          state[k++] = new_Tp;
                      //2.*M_PI*((time-time0)/l_Tp +timeshift -rS[i]/length) + M_PI*phaseShift
          Real Fshift = 2.*((-myFish->time0)/myFish->l_Tp +myFish->timeshift)+myFish->phaseShift;
          Fshift = fmod(Fshift,2.0);
          state[k++] = (Fshift<0) ? 2.+Fshift : Fshift;
          state[k++] = sr.VX*velscale;
          state[k++] = sr.VY*velscale;
          state[k++] = sr.AV*velscale;
      }

      state[k++] = sr.Dist*invlscale;
      state[k++] = sr.Quad;
      state[k++] = sr.VxAvg*velscale;
      state[k++] = sr.VyAvg*velscale;
      state[k++] = sr.AvAvg*velscale;
      state[k++] = sr.Pout*powerscale;
      state[k++] = sr.defPower*powerscale;
      state[k++] = sr.EffPDef;
      state[k++] = sr.PoutBnd*powerscale;
      state[k++] = sr.defPowerBnd*powerscale;
      state[k++] = sr.EffPDefBnd;
      state[k++] = sr.Pthrust*powerscale;
      state[k++] = sr.Pdrag*powerscale;
      state[k++] = sr.ToD;

      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelNAbove[j]*velscale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelTAbove[j]*velscale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelNBelow[j]*velscale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.VelTBelow[j]*velscale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FPAbove[j]*forcescale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FVAbove[j]*forcescale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FPBelow[j]*forcescale;
      for (int j=0; j<NpLatLine; j++) state[k++] = sr.FVBelow[j]*forcescale;
      for (int j=0;j<2*NpLatLine;j++) state[k++] = sr.raySight[j]*invlscale;

      //fflush(0);
      const Real reward = (sr.info==2) ? -10 : sr.EffPDefBnd;
      comm->sendState(iLabel, sr.info, state, reward); //TODO
      if (sr.info==2) return;
      sr.info = 0;

      comm->recvAction(actions);
      myFish->execute(time, sr.t_next_comm, actions);

      old_curv = new_curv;
      new_curv = actions[0];
      if(nActions>1) {
          new_Tp = actions[1];
          sr.t_next_comm += .5*myFish->l_Tp;
      } else if (nActions==1) {
         sr.t_next_comm += .5*myFish->Tperiod;
      }

      #ifndef __RL_TRAINING
      if(!rank) {
        printf("Next action of agent %d at time %g\n", iAgent, sr.t_next_comm);
        ofstream filedrag;
        filedrag.open(("orders_"+to_string(iAgent)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<new_curv;
        //if(nActions==2) filedrag<<" "<<new_Tp;
        filedrag<<endl;
        filedrag.close();
      }
			#endif

      sr.resetAverage();
      fflush(0);
    }
}
