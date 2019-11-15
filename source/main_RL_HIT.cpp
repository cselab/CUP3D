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
#include <sys/unistd.h> // hostname
#include <sstream>

#define FREQ_UPDATE 1
#define LES_HIT_RL_INIT_T 1
//#define SGSRL_STATE_SCALING

using Real = cubismup3d::Real;

struct TargetData
{
  std::vector<std::string> targetFiles_tokens;
  std::string active_token; // specifies eps/nu combo, taken from the list above

  size_t nModes;
  std::vector<double> logE_mean, mode, E_mean;
  std::vector<std::vector<double>> logE_invCov;

  double eps, nu, tKinEn, lambda, Re_lam, tInteg, lInteg, avg_Du, epsVis;

  void sampleParameters(std::mt19937& gen)
  {
    std::uniform_int_distribution<size_t> dist(0, targetFiles_tokens.size()-1);
    active_token = targetFiles_tokens[dist(gen)];
    readScalars(active_token); // for a (eps,nu) read TKE, lambda, epsVisc...
    readMeanSpectrum(active_token); // read also mean energy spectrum
    readInvCovSpectrum(active_token); // read also covariance matrix of logE
  }

  void readScalars(const std::string paramspec)
  {
    std::string line; char arg[32]; double stdev;
    std::ifstream file("../../scalars_" + paramspec);
    if (!file.is_open()) printf("scalars FILE NOT FOUND\n");

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le", arg, &eps);
    assert(strcmp(arg, "eps") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le", arg, &nu);
    assert(strcmp(arg, "nu") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &tKinEn, &stdev);
    assert(strcmp(arg, "tKinEn") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &lambda, &stdev);
    assert(strcmp(arg, "lambda") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &Re_lam, &stdev);
    assert(strcmp(arg, "Re_lam") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &tInteg, &stdev);
    assert(strcmp(arg, "tInteg") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &lInteg, &stdev);
    assert(strcmp(arg, "lInteg") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &avg_Du, &stdev);
    assert(strcmp(arg, "avg_Du") == 0);

    std::getline(file, line);
    sscanf(line.c_str(), "%s %le %le", arg, &epsVis, &stdev);
    assert(strcmp(arg, "epsVis") == 0);

    printf("Params eps:%e nu:%e with mean quantities: tKinEn=%e lambda=%e "
           "Re_lam=%e tInteg=%e lInteg=%e avg_Du=%e epsVis=%e\n", eps, nu,
            tKinEn, lambda, Re_lam, tInteg, lInteg, avg_Du, epsVis);
  }

  void readMeanSpectrum(const std::string paramspec)
  {
    std::string line;
    logE_mean.clear(); mode.clear();
    std::ifstream file("../../spectrumLogE_" + paramspec);
    if (!file.is_open()) printf("spectrumLogE FILE NOT FOUND\n");
    while (std::getline(file, line)) {
        mode.push_back(0); logE_mean.push_back(0);
        sscanf(line.c_str(), "%le, %le", & mode.back(), & logE_mean.back());
    }
    nModes = mode.size();
    //for (size_t i=0; i<nModes; ++i) printf("%f %f\n", mode[i], logE_mean[i]);
    E_mean = std::vector<double>(nModes);
    for(size_t i = 0; i<nModes; ++i) E_mean[i] = std::exp(logE_mean[i]);
  }

  void readInvCovSpectrum(const std::string paramspec)
  {
    std::string line;
    logE_invCov = std::vector<std::vector<double>>(
        nModes, std::vector<double>(nModes,0) );
    std::ifstream file("../../invCovLogE_" + paramspec);
    if (!file.is_open()) printf("invCovLogE FILE NOT FOUND\n");
    size_t j = 0;
    while (std::getline(file, line)) {
        size_t i = 0;
        std::istringstream linestream(line);
        while (std::getline(linestream, line, ','))
          logE_invCov[j][i++] = std::stof(line);
        assert(i==nModes);
        j++;
    }
    assert(j==nModes);
  }

  TargetData(std::string params) : targetFiles_tokens(readTargetSpecs(params))
  {
  }

  static std::vector<std::string> readTargetSpecs(std::string paramsList)
  {
    std::stringstream ss(paramsList);
    std::vector<std::string> tokens;
    std::string item;
    while (getline(ss, item, ',')) tokens.push_back(item);
    return tokens;
  }

  inline void updateReward(const cubismup3d::HITstatistics& stats,
                           const Real alpha, Real & reward)
  {
    std::vector<double> logE(stats.nBin);
    for (int i=0; i<stats.nBin; ++i) logE[i] = std::log(stats.E_msr[i]);

    double dev = 0;
    for (int j=0; j<stats.nBin; ++j)
     for (int i=0; i<stats.nBin; ++i)
      dev += (logE[i]-logE_mean[i])* logE_invCov[j][i] *(logE[j]-logE_mean[j]);

    // normalize with expectation of L2 squared norm of N(0,I) distrib vector:
    // E[X^2] = sum E[x^2] = sum Var[x] = trace I = nBin
    assert(dev >= 0);
    const double arg = dev / stats.nBin;
    const double newRew = arg > 0 ? 1 / arg : std::exp(1-arg);
    //printf("Rt : %f\n", newRew);
    reward = (1-alpha) * reward + alpha * newRew;
  }
};

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
  //int wrank; MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

  if (rank==0) {
    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    std::cout <<
    "=======================================================================\n";
    std::cout <<
    "Cubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
    std::cout <<
    "=======================================================================\n";
    #ifdef NDEBUG
        std::cout<<"Running on "<<hostname<<"in RELEASE mode!\n";
    #else
        std::cout<<"Running on "<<hostname<<"in DEBUG mode!\n";
    #endif
  }

  cubism::ArgumentParser parser(argc, argv);
  cubismup3d::Simulation sim(mpicom, parser);
  TargetData target(parser("-initCondFileTokens").asString());

  const int nActions = 1, nStates = 9;
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
  comm->setStateActionDims(nStates, nActions);
  comm->setNumAgents(nAgents + nThreadSafetyAgents);

  const std::vector<double> lower_act_bound{0.04}, upper_act_bound{0.06};
  comm->setActionScales(upper_act_bound, lower_act_bound, false);
  comm->disableDataTrackingForAgents(nAgents, nAgents + nThreadSafetyAgents);

  comm->finalizeProblemDescription(); // required for thread safety

  if( comm->isTraining() ) { // disable all dumping. //  && wrank != 1
    sim.sim.b3Ddump = false; sim.sim.muteAll  = true;
    sim.sim.b2Ddump = false; sim.sim.saveFreq = 0;
    sim.sim.verbose = false; sim.sim.saveTime = 0;
  }

  char dirname[1024]; dirname[1023] = '\0';
  unsigned sim_id = 0, tot_steps = 0;

  // Terminate loop if reached max number of time steps. Never terminate if 0
  while(true) // train loop
  {
    // avoid too many unneeded folders:
    if(sim_id == 0 || not comm->isTraining()) { //  || wrank == 1
      sprintf(dirname, "run_%08u/", sim_id);
      printf("Starting a new sim in directory %s\n", dirname);
      mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      chdir(dirname);
    }

    assert(sim.sim.spectralManip not_eq nullptr);
    const cubismup3d::HITstatistics& HTstats = sim.sim.spectralManip->stats;

    target.sampleParameters(comm->getPRNG());
    sim.sim.enInjectionRate = target.eps;
    sim.sim.nu = target.nu;
    sim.sim.spectralIC = "fromFile";
    sim.sim.initCondModes = target.mode;
    sim.sim.initCondSpectrum = target.E_mean;
    const Real tau_integral = target.tInteg;
    const Real tInit = LES_HIT_RL_INIT_T * target.tInteg;
    printf("Reset simulation up to time=%g with SGS for eps:%f nu:%f Re:%f\n",
           tInit, target.eps, target.nu, target.Re_lam);
    //const Real tau_eta       = HTstats.getKolmogorovT(eps, nu);

    while(true) { // initialization loop
      sim.reset();
      bool ICsuccess = true;
      while (sim.sim.time < tInit) {
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
    const Real timeUpdateLES = HTstats.getKolmogorovT(target.epsVis,target.nu)/2;
    sim.sim.sgs = "RLSM";

      //Real scaleGrads = -1;
      //updateGradScaling(sim.sim, HTstats, timeUpdateLES, scaleGrads);

    const auto nIntegralTime = 10;
    const int maxNumUpdatesPerSim= nIntegralTime * tau_integral / timeUpdateLES;
    cubismup3d::SGS_RL updateLES(sim.sim, comm, nAgentPerBlock);

    while (true) //simulation loop
    {
      const bool timeOut = step >= maxNumUpdatesPerSim;
      // even if timeOut call updateLES to send all the states of last step
      updateLES.run(sim.sim.dt, step==0, timeOut, target.tKinEn, target.epsVis, avgReward);
      if(timeOut) break;

      while ( time < (step+1)*timeUpdateLES )
      {
        const double dt = sim.calcMaxTimestep();
        time += dt;
        // Average reward over integral time:
        target.updateReward(HTstats, dt / tau_integral, avgReward);

        //updateGradScaling(sim.sim, HTstats, timeUpdateLES, scaleGrads);

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

    if(not comm->isTraining()) { //  || wrank == 1
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
