//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and
//  Hugues de Laroussilhe (huguesdelaroussilhe@gmail.com).
//

#include "smarties.h"

#include "Simulation.h"
#include "operators/SGS_RL.h"
#include "spectralOperators/SpectralManip.h"

#include <Cubism/ArgumentParser.h>

#include <sys/stat.h> // mkdir options
#include <unistd.h>  // chdir
#include <sstream>

#define FREQ_UPDATE 1
#define LES_HIT_RL_INIT_T 1
//#define SGSRL_STATE_SCALING

using Real = cubismup3d::Real;

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
    const double scale = (P_kde[idx-1] -P_kde[idx])/(E_kde[idx] -E_kde[idx-1]);
    double  ret = P_kde[idx] + scale * (E_kde[idx] - msr);

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
    fileName_scale = "../../scaleTarget.dat";
    fileName_mean  = "../../meanTarget.dat";
    fileName_kde   = "../../kdeTarget.dat";
    inFile.open(fileName_mean);
    if (!inFile) {
      std::cout<<"Cannot open "<<fileName_mean<<std::endl;
      abort();
    }
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
  if (!inFile){
    std::cout<<"Cannot open "<<fileName_kde<<std::endl;
    abort();
  }

  // First line is nSamplePerMode
  std::string line;
  std::getline(inFile, line);
  std::istringstream in(line);
  in >> ret.nSamplePerMode;

  size_t size = ret.nSamplePerMode * ret.nModes;
  std::vector<double> E_kde(size, 0.0);
  std::vector<double> P_kde(size, 0.0);

  for (int i=0; i<ret.nModes; i++) {
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
    return std::isnan(val) || std::isinf(val);
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

inline void updateReward(const cubismup3d::HITstatistics& stats,
                         const targetData& target,
                         const Real alpha, Real & reward)
{
  double newRew = 0;
  for (int i=0; i<stats.nBin; ++i) {
    //newRew += tgt.getRewardL2(i, msr[i]);
    newRew += target.getRewardKDE(i, stats.E_msr[i]);
  }
  reward = (1-alpha) * reward + alpha * newRew;
}

inline void updateReward(const cubismup3d::HITstatistics& stats,
                         const cubismup3d::SimulationData& sim,
                         const Real alpha, Real & reward)
{
  std::vector<Real> K = std::vector<Real> (stats.nBin, 0.0);
  std::vector<Real> Etgt = std::vector<Real> (stats.nBin, 0.0);
  stats.getTargetSpectrumFit(sim.enInjectionRate, sim.nu, K, Etgt);
  static constexpr Real lnEPS = std::log(std::numeric_limits<Real>::epsilon());

  double newRew = 0;
  for (int i=0; i<stats.nBin; i++)
  {
    const Real logEkTgt = std::log(Etgt[i]), logEk = std::log(stats.E_msr[i]);
    const Real logRk = std::pow((logEkTgt-logEk)/std::max(logEkTgt, lnEPS), 2);
    newRew += logRk;
  }
  newRew = std::exp(- newRew / stats.nBin );
  //printf("Rt : %f\n", newRew);
  reward = (1-alpha) * reward + alpha * newRew;
}

inline void updateGradScaling(const cubismup3d::SimulationData& sim,
                              const cubismup3d::HITstatistics& stats,
                              const Real timeUpdateLES,
                                    Real & scaling_factor)
{
  // Scale gradients with tau_eta corrected with eps_forcing and nu_sgs
  // tau_eta = (nu/eps)^0.5 but nu_tot  = nu + nu_sgs
  // and eps_tot = eps_nu + eps_numerics + eps_sgs = eps_forcing
  // tau_eta_corrected = tau_eta +
  //                   D(tau_eta)/Dnu * delta_nu + D(tau_eta)/Deps * delta_eps
  //                   = tau_eta * 1/2 * tau_eta *  delta_nu  / nu
  //                             - 1/2 * tau_eta * (delta_eps / eps) )
  //                   = tau_eta * (1 + nu_sgs/nu/2 - (eps_forcing-eps)/eps/2)
  // Average gradient scaling over time between LES updates
  const Real beta = sim.dt / timeUpdateLES;
  //const Real turbEnergy = stats.tke;
  const Real viscDissip = stats.eps, totalDissip = sim.dissipationRate;
  const Real avgSGS_nu = sim.nu_sgs, nu = sim.nu;
  const Real tau_eta_sim = std::sqrt(nu / viscDissip);
  const Real correction = 1.0 + 0.5 * (totalDissip - viscDissip) / viscDissip
                              + 0.5 * avgSGS_nu / nu;
                            //+ 0.5 * avgSGS_nu / std::sqrt(nu * viscDissip);
  const Real newScaling = tau_eta_sim * correction;
  if(scaling_factor < 0) scaling_factor = newScaling; // initialization
  else scaling_factor = (1-beta) * scaling_factor + beta * newScaling;
}

inline void sampleNuEps(Real & eps, Real & nu, std::mt19937& gen,
                        const cubismup3d::HITstatistics& stats)
{
  std::uniform_real_distribution<Real> disEps(std::log(0.01), std::log(2.0));
  std::uniform_real_distribution<Real> disNu(std::log(0.002), std::log(0.02));
  const Real refH = 2 * M_PI / 16 / 12;
  while (true) {
    eps = std::exp(disEps(gen));
    nu = std::exp(disNu(gen));
    if( stats.getHITReynoldsFit(eps,nu) < 100 &&
        stats.getHITReynoldsFit(eps,nu) >  20 &&
        stats.getTaylorMicroscaleFit(eps, nu) < 0.1 * stats.L &&
        stats.getKolmogorovL(eps, nu) < refH &&
        stats.getKolmogorovL(eps, nu) > refH/8 )
      break; // eps, nu within problem bounds
  }
}

inline void app_main(
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
  cubismup3d::Simulation sim(mpicom, parser);

  #ifdef CUP3D_USE_FILE_TARGET
    std::cout<<"Setting up target data"<<std::endl;
    targetData target = getTarget(sim.sim, comm->isTraining());
  #endif

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
  const int nAgentPerBlock = 1;
  const int nBlock=sim.sim.local_bpdx * sim.sim.local_bpdy * sim.sim.local_bpdz;
  const int nAgents = nBlock * nAgentPerBlock; // actual learning agents
  const int nThreadSafetyAgents = omp_get_max_threads();
  comm->set_state_action_dims(nStates, nActions);
  comm->set_num_agents(nAgents + nThreadSafetyAgents);

  const std::vector<double> lower_act_bound{0.04}, upper_act_bound{0.06};
  comm->set_action_scales(upper_act_bound, lower_act_bound, false);
  comm->disableDataTrackingForAgents(nAgents, nAgents + nThreadSafetyAgents);

  comm->finalize_problem_description(); // required for thread safety

  if( comm->isTraining() ) { // disable all dumping. comment out for dbg
    sim.sim.b3Ddump = false; sim.sim.muteAll  = true;
    sim.sim.b2Ddump = false; sim.sim.saveFreq = 0;
    sim.sim.verbose = false; sim.sim.saveTime = 0;
  }
  //sim.sim.icFromH5 = "../" + sim.sim.icFromH5; // bcz we will create a run dir

  char dirname[1024]; dirname[1023] = '\0';
  unsigned sim_id = 0, tot_steps = 0;

  // Terminate loop if reached max number of time steps. Never terminate if 0
  while(true) // train loop
  {
    if(sim_id == 0 || ! comm->isTraining()) { // avoid too many unneeded folders
      sprintf(dirname, "run_%08u/", sim_id); // by fast simulations when train
      printf("Starting a new sim in directory %s\n", dirname);
      mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      chdir(dirname);
    }

    assert(sim.sim.spectralManip not_eq nullptr);
    const cubismup3d::HITstatistics& HTstats = sim.sim.spectralManip->stats;

    Real eps, nu;
    while(true) { // initialization loop
      // TODO : FIGURE OUT BOUNDS:
      //sampleNuEps(eps, nu, comm->getPRNG(), HTstats);
      sim.sim.enInjectionRate = 0.1;
      sim.sim.nu = 0.003;
      sim.reset();
      //const double tStart = 2.5 * target.tau_integral; // CUP3D_USE_FILE_TARGET
      const Real tau_integral  = HTstats.getIntegralTimeFit(eps, nu);
      //const Real tau_eta       = HTstats.getKolmogorovT(eps, nu);
      printf("Reset simulation up to time=%g with SGS for eps:%f nu:%f Re:%f\n",
             LES_HIT_RL_INIT_T * tau_integral, eps, nu,
             HTstats.getHITReynoldsFit(eps, nu));
      bool ICsuccess = true;
      while (sim.sim.time < LES_HIT_RL_INIT_T * tau_integral) {
        sim.sim.sgs = "SSM";
        const double dt = sim.calcMaxTimestep();
        sim.timestep(dt);
        if ( isTerminal( sim.sim ) ) {
          ICsuccess = false;
          break;
        }
      }
      if( ICsuccess ) break;
      printf("failed, try new IC\n");
    }

    fflush(0);
    int step = 0;
    double time = 0;
    double avgReward  = 0;
    bool policyFailed = false;
    const Real timeUpdateLES = 0.5 * HTstats.getKolmogorovT(eps, nu);
    const Real tau_integral  = HTstats.getIntegralTimeFit(eps, nu);
    sim.sim.sgs = "RLSM";

    #ifdef SGSRL_STATE_SCALING
      Real scaleGrads = -1;
      updateGradScaling(sim.sim, HTstats, timeUpdateLES, scaleGrads);
    #else
      const Real scaleGrads = std::sqrt(nu / eps);
    #endif

    const unsigned int nIntegralTime = 10;
    const int maxNumUpdatesPerSim= nIntegralTime * tau_integral / timeUpdateLES;
    cubismup3d::SGS_RL updateLES(sim.sim, comm, nAgentPerBlock);

    while (true) //simulation loop
    {
      const bool timeOut = step >= maxNumUpdatesPerSim;
      // even if timeOut call updateLES to send all the states of last step
      updateLES.run(sim.sim.dt, step==0, timeOut, scaleGrads, avgReward);
      if(timeOut) break;

      while ( time < (step+1)*timeUpdateLES )
      {
        const double dt = sim.calcMaxTimestep();
        time += dt;
        // Average reward over integral time
        #ifdef CUP3D_USE_FILE_TARGET
          updateReward(HTstats,  target, dt / tau_integral, avgReward);
        #else
          updateReward(HTstats, sim.sim, dt / tau_integral, avgReward);
        #endif

        #ifdef SGSRL_STATE_SCALING
          updateGradScaling(sim.sim, HTstats, timeUpdateLES, scaleGrads);
        #endif

        if ( sim.timestep( dt ) ) { // if true sim has ended
          printf("Set -tend 0. This file decides the length of train sim.\n");
          assert(false); fflush(0); MPI_Abort(mpicom, 1);
        }
        if ( isTerminal( sim.sim ) ) {
          policyFailed = true;
          break;
        }
      }
      step++;
      tot_steps++;

      if ( policyFailed ) {
        // Agent gets penalized if the simulations blows up. For KDE reward,
        // penal is -0.5 max_reward * (n of missing steps to finish the episode)
        // WARNING: not consistent with L2 norm reward
        const std::vector<double> S_T(nStates, 0); // values in S_T dont matter
        const double R_T = (step - maxNumUpdatesPerSim);
        for(int i=0; i<nAgents; ++i) comm->sendTermState(S_T, R_T, i);
        break;
      }
    } // simulation is done

    if( not comm->isTraining() ) {
      chdir("../"); // matches previous if
    }
    sim_id++;
  }
}

int main(int argc, char**argv)
{
  smarties::Engine e(argc, argv);
  if( e.parse() ) return 1;
  e.run( app_main );
  return 0;
}
