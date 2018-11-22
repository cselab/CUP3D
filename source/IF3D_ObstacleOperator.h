//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef IF3D_ROCKS_IF3D_ObstacleOperator_h
#define IF3D_ROCKS_IF3D_ObstacleOperator_h

#include "Definitions.h"
//#include "IF3D_ObstacleLibrary.h"
#include "OperatorComputeForces.h"

#include <array>
#include <fstream>

// forward declaration of derived class for visitor

namespace cubism { class ArgumentParser; }

class IF3D_ObstacleOperator;
class IF3D_ObstacleVector;


/*
 * Return the size of the domain in the physical space.
 * All returned values are in the range (0, 1], and at least one is equal to 1.
 */
std::array<double, 3> getGridExtent(const FluidGridMPI &grid);


/*
 * Structure containing all externally configurable parameters of a base obstacle.
 */
struct ObstacleArguments {
  double length = 0.0;
  std::array<double, 3> position = {{0.0, 0.0, 0.0}};
  std::array<double, 4> quaternion = {{0.0, 0.0, 0.0, 0.0}};
  std::array<double, 3> enforcedVelocity = {{0.0, 0.0, 0.0}};  // Only if bForcedInSimFrame.
  std::array<bool, 3> bForcedInSimFrame = {{false, false, false}};
  std::array<bool, 3> bFixFrameOfRef = {{false, false, false}};
  bool bFixToPlanar = false;
  bool bComputeForces = true;

  ObstacleArguments() = default;

  /* Convert human-readable format into internal representation of parameters. */
  ObstacleArguments(const FluidGridMPI &grid, ArgumentParser &parser);
};


struct ObstacleVisitor
{
  virtual ~ObstacleVisitor() {}

  virtual void visit(IF3D_ObstacleOperator* const obstacle) = 0;
  //virtual void visit(IF3D_ObstacleVector  * const obstacle) {}
};

class IF3D_ObstacleOperator
{
protected:
  #ifdef RL_LAYER
    StateReward sr;
  #endif
  FluidGridMPI * grid;
  const Real*const _uInf;
  std::vector<BlockInfo> vInfo;
  std::map<int,ObstacleBlock*> obstacleBlocks;
  std::vector<double> sum;
  bool printedHeaderVels = false;
  bool isSelfPropelled = false;
public:
  int obstacleID=0, rank=0, size=0;
  bool bFixToPlanar=1, bInteractive=0, bHasSkin=0, bForces=0;
  double quaternion[4] = {1,0,0,0}, _2Dangle = 0, phaseShift=0; //orientation
  double position[3] = {0,0,0}, absPos[3] = {0,0,0}, transVel[3] = {0,0,0};
  double angVel[3] = {0,0,0}, volume = 0, J[6] = {0,0,0,0,0,0}; //mom of inertia
  //from diagnostics:
  double mass=0, force[3] = {0,0,0}, torque[3] = {0,0,0};
  //from compute forces: perimeter, circulation and forces
  double totChi=0, gamma[3]={0,0,0}, surfForce[3]={0,0,0};
  //pressure and viscous contribution from compute forces:
  double presForce[3]={0,0,0}, viscForce[3]={0,0,0}, surfTorque[3]={0,0,0};
  double drag=0, thrust=0, Pout=0, PoutBnd=0, pLocom=0;
  double defPower=0, defPowerBnd=0, Pthrust=0, Pdrag=0, EffPDef=0, EffPDefBnd=0;
  double transVel_correction[3]={0,0,0}, angVel_correction[3]={0,0,0}, length;
  //forced obstacles:
  double transVel_computed[3]= {0,0,0}, angVel_computed[3]= {0,0,0};
  double transVel_imposed[3]= {0,0,0};
  double ext_X, ext_Y, ext_Z;

  // stuff dealing with frame of reference:
  std::array<bool, 3> bFixFrameOfRef = {{false, false, false}};
  std::array<bool, 3> bForcedInSimFrame = {{false, false, false}};
  std::array<bool, 3> bBlockRotation = {{false, false, false}};
  bool isMPIBarrierOnChiCompute = false;
protected:
  virtual void _writeComputedVelToFile(const int step_id, const double t, const Real * uInf);
  virtual void _writeDiagForcesToFile(const int step_id, const double t);
  virtual void _writeSurfForcesToFile(const int step_id, const double t);
  void _makeDefVelocitiesMomentumFree(const double CoM[3]);
  void _computeUdefMoments(double lin_momenta[3], double ang_momenta[3], const double CoM[3]);
  //void _finalizeAngVel(Real AV[3], const Real J[6], const Real& gam0, const Real& gam1, const Real& gam2);

public:
  IF3D_ObstacleOperator(FluidGridMPI *g, const ObstacleArguments &args, const Real *uInf);
  IF3D_ObstacleOperator(FluidGridMPI *g, ArgumentParser &parser, const Real *uInf)
      : IF3D_ObstacleOperator(g, ObstacleArguments(*g, parser), uInf) { }

  IF3D_ObstacleOperator(FluidGridMPI * g, const Real *uInf = nullptr) : grid(g), _uInf(uInf)
  {
    MPI_Comm_rank(grid->getCartComm(),&rank);
    MPI_Comm_size(grid->getCartComm(),&size);
    vInfo = grid->getBlocksInfo();
    updateSRextents();
  }

  void updateSRextents()
  {
    std::array<double, 3> extent = getGridExtent(*grid);
    ext_X = extent[0];
    ext_Y = extent[1];
    ext_Z = extent[2];
    if(!rank) printf("Got sim extents %g %g %g\n", ext_X, ext_Y, ext_Z);
    #ifdef RL_LAYER
      sr.ext_X = ext_X;
      sr.ext_Y = ext_Y;
      sr.ext_Z = ext_Z;
    #endif
  }

  virtual void Accept(ObstacleVisitor * visitor);
  virtual Real getD() const {return length;}

  virtual void computeDiagnostics(const int stepID, const double time, const Real* Uinf, const double lambda) ;
  virtual void computeVelocities(const Real* Uinf);
  virtual void computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump);
  virtual void update(const int step_id, const double t, const double dt, const Real* Uinf);
  virtual void save(const int step_id, const double t, std::string filename = std::string());
  virtual void restart(const double t, std::string filename = std::string());

  // some non-pure methods
  virtual void create(const int step_id,const double time, const double dt, const Real *Uinf);
  virtual void computeChi(const int step_id, const double time, const double dt, const Real *Uinf, int& mpi_status);
  virtual void finalize(const int step_id,const double time, const double dt, const Real *Uinf);

  //methods that work for all obstacles
  std::map<int,ObstacleBlock*> getObstacleBlocks() const
  {
      return obstacleBlocks;
  }
  std::map<int,ObstacleBlock*>* getObstacleBlocksPtr()
  {
      return &obstacleBlocks;
  }

  void getObstacleBlocks(std::map<int,ObstacleBlock*>*& obstblock_ptr)
  {
      obstblock_ptr = &obstacleBlocks;
  }
  virtual void characteristic_function();

  virtual std::vector<int> intersectingBlockIDs(const int buffer) const;

  virtual ~IF3D_ObstacleOperator()
  {
    for(auto & entry : obstacleBlocks) {
      if(entry.second != nullptr) {
          delete entry.second;
          entry.second = nullptr;
      }
    }
    obstacleBlocks.clear();
  }

  virtual void getTranslationVelocity(double UT[3]) const;
  virtual void getAngularVelocity(double W[3]) const;
  virtual void getCenterOfMass(double CM[3]) const;
  virtual void setTranslationVelocity(double UT[3]);
  virtual void setAngularVelocity(const double W[3]);

  // driver to execute finite difference kernels either on all points relevant
  // to the mass of the obstacle (where we have char func) or only on surface
  enum INTEGRAL { VOLUME, SURFACE };
  template <typename Kernel, INTEGRAL integral>
  void compute(const std::vector<Kernel*>& kernels)
  {
    SynchronizerMPI<Real>& Synch = grid->sync(*(kernels[0]));

    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i)  labs[i].prepare(*grid, Synch);

    MPI_Barrier(grid->getCartComm());
    std::vector<BlockInfo> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel num_threads(nthreads)
    {
      const int tid = omp_get_thread_num();
      Kernel& kernel=*(kernels[tid]);
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic, 1)
      for(int i=0; i<Ninner; i++)
      {
        BlockInfo info = avail0[i];
        const auto pos = obstacleBlocks.find(info.blockID);
        if(pos == obstacleBlocks.end()) continue;
        if(integral==SURFACE && pos->second->nPoints == 0) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        lab.load(info, 0);
        kernel(lab, info, b, pos->second);
      }
    }

    std::vector<BlockInfo> avail1 = Synch.avail_halo();
    const int Nhalo = avail1.size();

    #pragma omp parallel num_threads(nthreads)
    {
      const int tid = omp_get_thread_num();
      Kernel& kernel=*(kernels[tid]);
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic, 1)
      for(int i=0; i<Nhalo; i++)
      {
        BlockInfo info = avail1[i];
        const auto pos = obstacleBlocks.find(info.blockID);
        if(pos == obstacleBlocks.end()) continue;
        if(integral==SURFACE && pos->second->nPoints == 0) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        lab.load(info, 0);
        kernel(lab, info, b, pos->second);
      }
    }

    if(labs!=NULL) {
      delete [] labs;
      labs=NULL;
    }

    MPI_Barrier(grid->getCartComm());
  }


  #ifdef RL_LAYER
    StateReward* _getData() { return &sr; }
    virtual void execute(const int iAgent, const double time, const std::vector<double> action);

    virtual void getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL,
                                Real*& pYL, Real*& pXU, Real*& pYU, int& Npts);

    virtual void interpolateOnSkin(const double time, const int stepID, bool dumpWake=false);
  #endif
};

#endif
