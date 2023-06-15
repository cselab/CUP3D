//
//  CubismUP_3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Definitions.h"
#include <array>
#include <vector>
#include <random>
#include <memory>

namespace cubism {
  class Profiler;
  class ArgumentParser;
}

CubismUP_3D_NAMESPACE_BEGIN

class Operator;
class Obstacle;
class ObstacleVector;
class PoissonSolverBase;

struct SimulationData
{
  // MPI
  MPI_Comm comm;
  int rank;

  // Profiler
  cubism::Profiler * profiler = nullptr;

  // declare grids
  ScalarGrid * chi  = nullptr;
  ScalarGrid * pres = nullptr;
  VectorGrid * vel  = nullptr;
  VectorGrid * tmpV = nullptr;
  ScalarGrid * lhs  = nullptr;

  // mesh refinement
  ScalarAMR *  chi_amr;
  ScalarAMR * pres_amr;
  VectorAMR *  vel_amr;
  VectorAMR * tmpV_amr;
  ScalarAMR *  lhs_amr;

  // Get blocks on current rank
  inline std::vector<cubism::BlockInfo>&  chiInfo() const {return  chi->getBlocksInfo();}
  inline std::vector<cubism::BlockInfo>& presInfo() const {return pres->getBlocksInfo();}
  inline std::vector<cubism::BlockInfo>&  velInfo() const {return  vel->getBlocksInfo();}
  inline std::vector<cubism::BlockInfo>& tmpVInfo() const {return tmpV->getBlocksInfo();}
  inline std::vector<cubism::BlockInfo>&  lhsInfo() const {return  lhs->getBlocksInfo();}

  // Container holding the obstacles
  ObstacleVector * obstacle_vector = nullptr;

  // Operator Pipeline
  std::vector<std::shared_ptr<Operator>> pipeline;

  // Pressure solver for PressureProjection
  std::shared_ptr<PoissonSolverBase> pressureSolver;

  // Timestepping
  Real dt      = 0;//current timestep
  Real dt_old  = 0;//previous timestep
  Real CFL     = 0;//Courant number
  Real time    = 0;//current time
  int step     = 0;//currect step number
  Real endTime = 0;//stop simulation at t=endTime (=0 means inactive)
  int nsteps   = 0;//stop simulation after nsteps (=0 means inactive)
  int rampup;        //exponential CFL rampup for the first 'rampup' steps
  int step_2nd_start;//explicit Euler for the first 'step_2nd_start' steps (to initialize u_{n-1} for n=1)
  Real coefU[3] = {1.5,-2.0,0.5};//used for 2nd order time integration of obstacle positions

  //AMR & simulation domain
  int bpdx, bpdy, bpdz;              //blocks per dimension at refinement level 0
  int levelStart;                    //initial refinement level
  int levelMax;                      //max refinement level
  Real Rtol;                         //mesh refinement tolerance
  Real Ctol;                         //mesh compression tolerance
  std::array<Real, 3> extents;       //simulation cubic domain extents
  Real maxextent ;                   //max(extents[0],extents[1],extents[2])
  Real hmin, hmax;                   //max and min grid spacing
  std::array<Real, 3> uinf = {0,0,0};//velocity of Frame of Reference
  int levelMaxVorticity;             //mesh refinement due to vorticity magnitude is allowed only up to levelMaxVorticity levels (default value is levelMax)

  //Other stuff
  Real uMax_measured = 0;         //max velocity at current timestep
  Real uMax_allowed;              //if uMax_measured > uMax_allowed simulation will abort
  Real nu;                        //fluid kinematic viscosity
  Real lambda;                    //penalisation coefficient
  bool bImplicitPenalization = true;//explicit/implicit Penalisation
  Real DLM=0;                     // if DLM>0 then lambda=DLM/dt
  Real PoissonErrorTol;           //Poisson solver absolute error tolerance
  Real PoissonErrorTolRel;        //Poisson solver relative error tolerance
  std::string poissonSolver;      // Type of backed solver for poisson equation
  bool bCollision = false;          //indicator for collision between obstacles
  std::vector<int> bCollisionID;    //vector with indices of colliding obstacles
  BCflag BCx_flag = freespace;      //boundary conditions in X
  BCflag BCy_flag = freespace;      //boundary conditions in Y
  BCflag BCz_flag = freespace;      //boundary conditions in Z
  int bMeanConstraint = 1;          // if 0, no zero mean constraint for Poisson
                                    // if 1, replace one equation with zero mean contraint
				    // if 2, add mean to LHS. Ie. solve nabla^2 P + mean(P) = RHS
				    // which is mathematically equivalent to the zero-mean 
				    // solution of nabla^2 P = RHS (for zero Neumann BCs of course)
				    // if >2, we set one grid point to p=0
  bool StaticObstacles = false; //if true, obstacles do not change shape or location.
                                //Useful when obstacle creation is expensive (for ExternalObstacle for example)
                                //as it skips obstacle creation whenever the mesh is not refined
  bool MeshChanged = true; //becomes true when the mesh is refined/compressed or when load balancing occurs

  // Initial conditions
  std::string initCond = "zero";

  // uMax Channel flow
  Real uMax_forced = 0;
  bool bFixMassFlux = false;

  // Dump Settings
  int freqDiagnostics = 0;
  int freqProfiler = 0;
  int saveFreq = 0;
  bool bDump = false;
  bool verbose = false;
  bool muteAll = false;
  Real dumpTime=0;
  Real nextSaveTime=0;
  std::string path4serialization = "./";
  bool dumpP;
  bool dumpChi;
  bool dumpOmega,dumpOmegaX,dumpOmegaY,dumpOmegaZ;
  bool dumpVelocity,dumpVelocityX,dumpVelocityY,dumpVelocityZ;

  bool implicitDiffusion;
  Real DiffusionErrorTol;
  Real DiffusionErrorTolRel;

  void startProfiler(std::string name) const;
  void stopProfiler() const;
  void printResetProfiler();
  void _preprocessArguments();
  void writeRestartFiles();
  void readRestartFiles();
  ~SimulationData();
  SimulationData() = delete;
  SimulationData(const SimulationData &) = delete;
  SimulationData(SimulationData &&) = delete;
  SimulationData &operator=(const SimulationData &) = delete;
  SimulationData &operator=(SimulationData &&) = delete;
  SimulationData(MPI_Comm mpicomm, cubism::ArgumentParser &parser);
};

CubismUP_3D_NAMESPACE_END
