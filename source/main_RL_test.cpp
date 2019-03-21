//
//  main.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <atomic>

#include "Communicator.h"
#include "Simulation.h"
#include "operators/SGS_RL.h"
#include "Cubism/ArgumentParser.h"

#include "mpi.h"
#define FREQ_UPDATE 1

inline bool isTerminal(cubismup3d::SimulationData& sim)
{
  std::atomic<bool> terminal { false };
  const auto& vInfo = sim.vInfo();
  const auto isNotValid = [](const Real val) {
    return std::isnan(val) || std::fabs(val)>1e3;
  };

  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < vInfo.size(); ++i) {
    const cubismup3d::FluidBlock&b= *(cubismup3d::FluidBlock*)vInfo[i].ptrBlock;
    for(int iz=0; iz<cubismup3d::FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<cubismup3d::FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<cubismup3d::FluidBlock::sizeX; ++ix)
      if ( isNotValid(b(ix,iy,iz).u) || isNotValid(b(ix,iy,iz).v) ||
           isNotValid(b(ix,iy,iz).w) || isNotValid(b(ix,iy,iz).p) )
        terminal = true;
  }
  return terminal.load();
}

inline double updateReward(cubismup3d::SimulationData&sim, const double oldRew)
{
  const Real alpha = 0.999; //there are formulas to determne size of integration window
  const Real newReward = 1; // dummy
  return alpha*oldRew + (1-alpha) * newReward;
}


int app_main(
  Communicator*const comm, // communicator with smarties
  MPI_Comm mpicom,         // mpi_comm that mpi-based apps can use
  int argc, char**argv,    // arguments read from app's runtime settings file
  const unsigned numSteps  // number of time steps to run before exit
)
{
  // print received arguments:
  for(int i=0; i<argc; i++) {printf("arg: %s\n",argv[i]); fflush(0);}

  #ifdef CUP_ASYNC_DUMP
    const auto SECURITY = MPI_THREAD_MULTIPLE;
  #else
    const auto SECURITY = MPI_THREAD_FUNNELED;
  #endif
  int provided; MPI_Query_thread(&provided);
  if (provided < SECURITY ) {
    printf("ERROR: MPI implementation does not have required thread support\n");
    fflush(0); MPI_Abort(mpicom, 1);
  }
  int rank; MPI_Comm_rank(mpicom, &rank);

  if (rank==0) {
    std::cout <<
    "=======================================================================\n";
    std::cout <<
    "Cubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
    std::cout <<
    "=======================================================================\n";
    #ifdef NDEBUG
        std::cout << "Running in RELEASE mode!\n";
    #else
        std::cout << "Running in DEBUG mode!\n";
    #endif
  }

  cubism::ArgumentParser parser(argc, argv);

  cubismup3d::Simulation *sim = new cubismup3d::Simulation(mpicom, parser);

  const int nActions = 1, nStates = 27;
  // BIG TROUBLE WITH NAGENTS!
  // If every grid point is an agent: probably will allocate too much memory
  // and crash because smarties allocates a trajectory for each point
  // If only one agent: sequences will be garbled together and cannot
  // send clean Sequences.
  // Also, rememebr that separate agents are thread safe!
  // let's say that each fluid block has one agent
  const int nValidAgents = sim->sim.local_bpdx * sim->sim.local_bpdy * sim->sim.local_bpdz;
  const int nThreadSafetyAgents = omp_get_num_threads();
  comm->update_state_action_dims(nStates, nActions);
  comm->setnAgents(nValidAgents + nThreadSafetyAgents);
  comm->disableDataTrackingForAgents(nValidAgents, nValidAgents + nThreadSafetyAgents);
  const std::vector<double> lower_act_bound{-0.01}, upper_act_bound{0.01};
  comm->set_action_scales(upper_act_bound, lower_act_bound, false);

  if( comm->isTraining() ) { // disable all dumping. comment out for dbg
    sim->sim.b3Ddump = false; sim->sim.muteAll = true;
    sim->sim.b2Ddump = false; sim->sim.saveFreq = 0;
    sim->sim.verbose = false; sim->sim.saveTime = 0;
  }
  char dirname[1024]; dirname[1023] = '\0';
  unsigned sim_id = 0, tot_steps = 0;

  // Terminate loop if reached max number of time steps. Never terminate if 0
  while( numSteps == 0 || tot_steps<numSteps ) // train loop
  {
    if( not comm->isTraining() ) { // avoid too many unneeded folders created
      sprintf(dirname, "run_%08u/", sim_id); // by fast simulations when train
      printf("Starting a new sim in directory %s\n", dirname);
      mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      chdir(dirname);
    }

    double time = 0;
    unsigned step = 0;
    double avgReward = 0;
    bool policyFailed = false;
    sim->reset(); // TODO IMPLEMENT PROPER RESET (EG TIME AND STEP)

    while (true) //simulation loop
    {
      // Assume Re = halfLy * U_mean / nu ~ 2500. Then for channel Re_\tau = 180
      // Also: L_\tau = halfLy/ Re_\tau and T_\tau = L_\tau / U_mean
      // next line assumes U_mean = 1 and halfLy = 1. TODO generalize:
      const double freqUpdateLES = 1.0 / 180.0;
      // A simulation is composed of how many LES updates?
      const unsigned maxNumUpdatesPerSim = 1000; // random number... TODO

      const bool timeOut = step >= maxNumUpdatesPerSim;
      // even if timeOut call updateLES to send all the states of last step:
      cubismup3d::SGS_RL updateLES(sim->sim, comm, timeOut, avgReward);
      updateLES.run();
      if(timeOut) break;

      while ( time < (step+1)*freqUpdateLES )
      {
        const double dt = sim->calcMaxTimestep();
        time += dt;

        if ( sim->timestep( dt ) ) { // if true sim has ended
          printf("Set -tend 0. This file decides the length of train sim.\n");
          assert(false); fflush(0); MPI_Abort(mpicom, 1);
        }
        if ( isTerminal( sim->sim ) ) {
          policyFailed = true;
          break;
        }
      }
      step++;
      tot_steps++;

      avgReward = updateReward(sim->sim, avgReward);

      if ( policyFailed )
      {
        printf("Policy failed\n"); fflush(0);
        const std::vector<double> S_T(nStates, 0); // values in S_T dont matter
        const double R_T = -1000;
        for(int i=0; i<nValidAgents; i++) comm->sendTermState(S_T, R_T, i);
        break;
      }
    } // simulation is done

    if( not comm->isTraining() ) chdir("../"); // matches previous if
    sim_id++;
  }

  delete sim;
  return 0;
}
