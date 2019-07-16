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
#include <unistd.h> // chdir
#include <sys/stat.h> // mkdir options

#include "Communicators/Communicator.h"
#include "Simulation.h"
#include "operators/SGS_RL.h"
#include "operators/SpectralAnalysis.h"
#include "Cubism/ArgumentParser.h"

#include "mpi.h"
#define FREQ_UPDATE 1

struct targetData
{
  int nModes;
  std::vector<double> E_mean;
  std::vector<double> E_std2;

  int nSamplePerMode;
  std::vector<double> E_kde;
  std::vector<double> P_kde;

  double tau_integral;
  double tau_eta;
  double max_rew;

  double getRewardL2(const int modeID, const double msr) const
  {
    double r = (E_mean[modeID] - msr)/E_mean[modeID];
    return -r*r;
  }

  double getRewardKDE(const int modeID, const double msr) const
  {
    if (msr <= E_kde[modeID*nSamplePerMode] or
        msr >= E_kde[modeID*nSamplePerMode + nSamplePerMode - 1]) {
      return 0;
    }
    //size_t lower = modeID*nSamplePerMode;
    //size_t upper = modeID*nSamplePerMode + nSamplePerMode - 1;


    size_t idx = -1;
    for (int i=0; i<nSamplePerMode; i++){
      idx = modeID*nSamplePerMode + i;
      if (msr < E_kde[idx])
        break;
    }
    double  ret = P_kde[idx] + (P_kde[idx-1] - P_kde[idx]) / (E_kde[idx] - E_kde[idx-1]) * (E_kde[idx] - msr);

    return ret;
  }
};


// Here I assume that the target files are given for the same modes
// as the one that we will measure during training.
inline targetData getTarget(cubismup3d::SimulationData& sim, const bool bTrain)
{
  printf("Setting up target data\n");
  targetData ret;

  std::vector<double> E_mean;
  std::vector<double> E_std2;


  std::ifstream inFile;
  std::string fileName_scale = "../scaleTarget.dat";
  std::string fileName_mean  = "../meanTarget.dat";
  std::string fileName_kde   = "../kdeTarget.dat";

  // First get mean and std of the energy spectrum
  inFile.open(fileName_mean);

  if (!inFile) {
    std::cout<<"Cannot open "<<fileName_mean<<std::endl;
    abort();
  }


  const long maxGridSize = std::max({sim.bpdx * cubismup3d::FluidBlock::sizeX,
                                     sim.bpdy * cubismup3d::FluidBlock::sizeY,
                                     sim.bpdz * cubismup3d::FluidBlock::sizeZ});
  ret.nModes =  (int) maxGridSize/2 -1;

  int n=0;
  for(std::string line; std::getline(inFile, line); )
  {
    if (n>ret.nModes) break;
    std::istringstream in(line);
    Real k, mean, std2;
    in >> k >> mean >> std2;
    E_mean.push_back(mean);
    E_std2.push_back(std2);
    n++;
  }

  ret.E_mean = E_mean;
  ret.E_std2 = E_std2;
  inFile.close();


  // Get the Kernel Density Estimations from 'kdeTarget.dat'
  inFile.open(fileName_kde);
  std::string line;
  if (!inFile){
    std::cout<<"Cannot open "<<fileName_kde<<std::endl;
    abort();
  }

  // First line is nSamplePerMode
  std::getline(inFile, line);
  std::istringstream in(line);
  in >> ret.nSamplePerMode;

  size_t size = ret.nSamplePerMode * ret.nModes;
  std::vector<double> E_kde(size, 0.0);
  std::vector<double> P_kde(size, 0.0);

  for (int i=0; i<ret.nModes; i++){
    std::getline(inFile, line);
    std::getline(inFile, line);
    for (int j=0; j<ret.nSamplePerMode; j++){
      std::getline(inFile, line);
      std::istringstream iss(line);
      size_t idx = i * ret.nSamplePerMode + j;
      iss >> E_kde[idx] >> P_kde[idx];
    }
  }

  ret.E_kde = E_kde;
  ret.P_kde = P_kde;
  inFile.close();

  // Get tau_integral and tau_eta from 'scaleTarget.dat'
  inFile.open(fileName_scale);
  if (!inFile){
    std::cout<<"Cannot open "<<fileName_scale<<std::endl;
    abort();
  }
  std::getline(inFile, line);
  in = std::istringstream(line);
  in >> ret.tau_integral;
  std::getline(inFile, line);
  in = std::istringstream(line);
  in >> ret.tau_eta;
  std::getline(inFile, line);
  in = std::istringstream(line);
  in >> ret.max_rew;


  return ret;
}

inline bool isTerminal(cubismup3d::SimulationData& sim)
{
  std::atomic<bool> bSimValid { true };

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
        bSimValid = false;
  }
  int isSimValid = bSimValid.load() == true? 1 : 0; // just for clarity
  MPI_Allreduce(MPI_IN_PLACE, &isSimValid, 1, MPI_INT, MPI_PROD, sim.grid->getCartComm());
  return isSimValid == 0;
}

inline double updateReward(double oldRew, const int nBin, const targetData tgt, const double msr[], const double alpha)
{
  double newRew = 0;
  for (int i=0; i<nBin; i++){
    //newRew += tgt.getRewardL2(i, msr[i]);
    newRew += tgt.getRewardKDE(i, msr[i]);
  }
  return oldRew=0 ? newRew : (1-alpha)*oldRew + alpha*newRew;
}


int app_main(
  smarties::Communicator*const comm, // communicator with smarties
  MPI_Comm mpicom,                  // mpi_comm that mpi-based apps can use
  int argc, char**argv             // args read from app's runtime settings file
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

  std::cout<<"Setting up target data"<<std::endl;
  targetData target = getTarget(sim->sim, comm->isTraining());
  //targetData target;
  #ifdef SGSRL_STATE_INVARIANTS
    const int nActions = 1, nStates = 5;
  #else
    const int nActions = 1, nStates = 9;
  #endif
  // BIG TROUBLE WITH NAGENTS!
  // If every grid point is an agent: probably will allocate too much memory
  // and crash because smarties allocates a trajectory for each point
  // If only one agent: sequences will be garbled together and cannot
  // send clean Sequences.
  // Also, rememebr that separate agents are thread safe!
  // let's say that each fluid block has one agent
  const int nAgentPerBlock = sim->sim.nAgentsPerBlock;
  const int nBlocks = sim->sim.local_bpdx * sim->sim.local_bpdy * sim->sim.local_bpdz;
  const int nValidAgents = nBlocks * nAgentPerBlock;
  const int nThreadSafetyAgents = omp_get_max_threads();
  comm->set_state_action_dims(nStates, nActions);
  comm->set_num_agents(nValidAgents + nThreadSafetyAgents);

  const std::vector<double> lower_act_bound{0.04}, upper_act_bound{0.08};
  comm->set_action_scales(upper_act_bound, lower_act_bound, false);
  comm->disableDataTrackingForAgents(nValidAgents, nValidAgents + nThreadSafetyAgents);

  comm->finalize_problem_description(); // required for thread safety

  if( comm->isTraining() ) { // disable all dumping. comment out for dbg
    sim->sim.b3Ddump = false; sim->sim.muteAll  = true;
    sim->sim.b2Ddump = false; sim->sim.saveFreq = 0;
    sim->sim.verbose = false; sim->sim.saveTime = 0;
  }
  char dirname[1024]; dirname[1023] = '\0';
  unsigned sim_id = 0, tot_steps = 0;

  // Terminate loop if reached max number of time steps. Never terminate if 0
  while(true) // train loop
  {
    if( not comm->isTraining()) { // avoid too many unneeded folders created
      sprintf(dirname, "run_%08u/", sim_id); // by fast simulations when train
      printf("Starting a new sim in directory %s\n", dirname);
      mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      chdir(dirname);
    }

    double time = 0.0;
    double tStart = 5.0;
    int step = 0;
    double avgReward  = 0;

    bool policyFailed = false;
    sim->reset(comm->getPRNG(), tStart, nAgentPerBlock, sim_id, comm->isTraining()); // TODO Does not work for obstacles
    cubismup3d::SpectralAnalysis* sA = new cubismup3d::SpectralAnalysis(sim->sim);

    const Real tau_eta       = target.tau_eta;
    const Real tau_integral  = target.tau_integral;
    const Real timeUpdateLES = 0.5*tau_eta;

    const unsigned int nIntegralTime = 10;
    const int maxNumUpdatesPerSim = (int) (nIntegralTime * tau_integral / timeUpdateLES);

    while (true) //simulation loop
    {
      const bool timeOut = step >= maxNumUpdatesPerSim;
      // even if timeOut call updateLES to send all the states of last step
      bool evalStep = true;
      cubismup3d::SGS_RL updateLES(sim->sim, comm, step, timeOut, evalStep, avgReward, nAgentPerBlock);
      updateLES.run();
      if(timeOut) break;

      while ( time < (step+1)*timeUpdateLES )
      {
        const double dt = sim->calcMaxTimestep();
        time += dt;
        sA->run();
        const double alpha = dt / tau_integral;
        avgReward = updateReward(avgReward, sA->nBin, target, sA->E_msr, alpha);
        if ( sim->timestep( dt ) ) { // if true sim has ended
          printf("Set -tend 0. This file decides the length of train sim.\n");
          assert(false); fflush(0); MPI_Abort(mpicom, 1);
        }
        if ( isTerminal( sim->sim )) {
          policyFailed = true;
          break;
        }
      }
      step++;
      tot_steps++;

      if ( policyFailed )
      {
        const std::vector<double> S_T(nStates, 0); // values in S_T dont matter
        const double R_T = (double) (0.5*target.max_rew*(step - maxNumUpdatesPerSim));
        for(int i=0; i<nValidAgents; i++) comm->sendTermState(S_T, R_T, i);
        break;
      }
    } // simulation is done
    delete sA;
    if( not comm->isTraining()) {
      chdir("../"); // matches previous if
    }
    sim_id++;
  }

  delete sim;
  return 0;
}
