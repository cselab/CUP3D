//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef IF3D_ROCKS_IF3D_ObstacleOperator_h
#define IF3D_ROCKS_IF3D_ObstacleOperator_h

#include "Definitions.h"
//#include "IF3D_ObstacleLibrary.h"
#include "IF2D_FactoryFileLineParser.h"
#include "OperatorComputeForces.h"
#include <fstream>

// forward declaration of derived class for visitor

class IF3D_ObstacleOperator;
class IF3D_ObstacleVector;

struct ObstacleVisitor
{
    virtual void visit(IF3D_ObstacleOperator* const obstacle) = 0;
  //virtual void visit(IF3D_ObstacleVector  * const obstacle) {}
};

class IF3D_ObstacleOperator
{
protected:
  StateReward sr;
  FluidGridMPI * grid;
  const Real*const _uInf;
  vector<BlockInfo> vInfo;
  std::map<int,ObstacleBlock*> obstacleBlocks;
  vector<double> sum;

public:
  int obstacleID=0, rank, size;
  bool bFixToPlanar=1, bInteractive=0, bHasSkin=0, bForces=0;
  double quaternion[4] = {1,0,0,0}, _2Dangle = 0, phaseShift=0; //orientation
  double position[3] = {0,0,0}, absPos[3] = {0,0,0}, transVel[3] = {0,0,0};
  double angVel[3] = {0,0,0}, volume = 0, J[6] = {0,0,0,0,0,0}; //mom of inertia
  //from diagnostics:
  double mass=0, force[3] = {0,0,0}, torque[3] = {0,0,0};
  //from compute forces:
  double totChi=0, surfForce[3]={0,0,0}, drag=0, thrust=0, Pout=0, PoutBnd=0;
  double defPower=0, defPowerBnd=0, Pthrust=0, Pdrag=0, EffPDef=0, EffPDefBnd=0;
  double transVel_correction[3]={0,0,0}, angVel_correction[3]={0,0,0}, length;
  //forced obstacles:
  double transVel_computed[3]= {0,0,0}, angVel_computed[3]= {0,0,0};
  double ext_X, ext_Y, ext_Z;
  double torqueZsection=0.0; 
 
  // stuff dealing with frame of reference:
  bool bFixFrameOfRef[3] = {false, false, false};
  bool bForcedInSimFrame[3] = {false, false, false};

protected:
  virtual void _parseArguments(ArgumentParser & parser);
  virtual void _writeComputedVelToFile(const int step_id, const double t, const Real * uInf);
  virtual void _writeDiagForcesToFile(const int step_id, const double t);
  void _makeDefVelocitiesMomentumFree(const double CoM[3]);
  void _computeUdefMoments(double lin_momenta[3], double ang_momenta[3], const double CoM[3]);
  //void _finalizeAngVel(Real AV[3], const Real J[6], const Real& gam0, const Real& gam1, const Real& gam2);

public:
  IF3D_ObstacleOperator(FluidGridMPI*g, ArgumentParser&parser, const Real*const uInf) : grid(g), _uInf(uInf)
  {
    MPI_Comm_rank(grid->getCartComm(),&rank);
    MPI_Comm_size(grid->getCartComm(),&size);
    vInfo = grid->getBlocksInfo();
    const double extent = 1;//grid->maxextent;
    const double NFE[3] = {
        (double)grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
        (double)grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
        (double)grid->getBlocksPerDimension(2)*FluidBlock::sizeZ,
    };
    const double maxbpd = max(NFE[0], max(NFE[1], NFE[2]));
    const double scale[3] = { NFE[0]/maxbpd, NFE[1]/maxbpd, NFE[2]/maxbpd };
    sr.ext_X = ext_X = scale[0]*extent;
    sr.ext_Y = ext_Y = scale[1]*extent;
    sr.ext_Z = ext_Z = scale[2]*extent;
    if(!rank)
    printf("Got sim extents %g %g %g\n", sr.ext_X, sr.ext_Y, sr.ext_Z);
    _parseArguments(parser);
  }

  IF3D_ObstacleOperator(FluidGridMPI * g) : grid(g), _uInf(nullptr)
  {
    MPI_Comm_rank(grid->getCartComm(),&rank);
    MPI_Comm_size(grid->getCartComm(),&size);
    vInfo = grid->getBlocksInfo();
    const double extent = 1;//grid->maxextent;
    const double NFE[3] = {
        (double)grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
        (double)grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
        (double)grid->getBlocksPerDimension(2)*FluidBlock::sizeZ,
    };
    const double maxbpd = max(NFE[0], max(NFE[1], NFE[2]));
    const double scale[3] = { NFE[0]/maxbpd, NFE[1]/maxbpd, NFE[2]/maxbpd };
    sr.ext_X = ext_X = scale[0]*extent;
    sr.ext_Y = ext_Y = scale[1]*extent;
    sr.ext_Z = ext_Z = scale[2]*extent;
    if(!rank)
    printf("Got sim extents %g %g %g\n", sr.ext_X, sr.ext_Y, sr.ext_Z);
  }

  virtual void Accept(ObstacleVisitor * visitor);
  virtual void dumpWake(const int stepID, const double t, const Real* Uinf);
  virtual Real getD() const {return length;}

  virtual void computeDiagnostics(const int stepID, const double time, const Real* Uinf, const double lambda) ;
  virtual void computeVelocities(const Real* Uinf);
  virtual void computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump);
  virtual void update(const int step_id, const double t, const double dt, const Real* Uinf);
  virtual void save(const int step_id, const double t, std::string filename = std::string());
  virtual void restart(const double t, std::string filename = std::string());
  virtual void interpolateOnSkin(const double time, const int stepID, bool dumpWake=false);
  virtual void execute(const int iAgent, const double time, const vector<double> action);
  StateReward* _getData() { return &sr; }
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

  virtual void getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL,
                              Real*& pYL, Real*& pXU, Real*& pYU, int& Npts);
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

  enum INTEGRAL { VOLUME, SURFACE };
  template <typename Kernel, INTEGRAL integral>
  void compute(const vector<Kernel*>& kernels)
  {
    SynchronizerMPI& Synch = grid->sync(*(kernels[0]));

    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    for(int i = 0; i < nthreads; ++i)
      labs[i].prepare(*grid, Synch);

    MPI_Barrier(grid->getCartComm());
    vector<BlockInfo> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel num_threads(nthreads)
    {
      const int tid = omp_get_thread_num();
      Kernel& kernel=*(kernels[tid]);
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic,1)
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

    vector<BlockInfo> avail1 = Synch.avail_halo();
    const int Nhalo = avail1.size();

    #pragma omp parallel num_threads(nthreads)
    {
      const int tid = omp_get_thread_num();
      Kernel& kernel=*(kernels[tid]);
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic,1)
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
};

#endif
