//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef CubismUP_3D_Fish_h
#define CubismUP_3D_Fish_h

#include "Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

class FishMidlineData;
struct VolumeSegment_OBB;

class Fish: public Obstacle
{
 protected:
  FishMidlineData * myFish = nullptr;
  // Arguments read from parser
  double Tperiod, phaseShift;
  bool bCorrectTrajectory, bCorrectPosition;
  // Rest
  double volume_internal=0, J_internal=0;
  double CoM_internal[2]={0,0}, vCoM_internal[2]={0,0};
  double theta_internal=0, angvel_internal=0, angvel_internal_prev=0;
  double angvel_integral[3] = {0,0,0};


  void integrateMidline();
  //void apply_pid_corrections();

  // first how to create blocks of segments:
  typedef std::vector<VolumeSegment_OBB> vecsegm_t;
  vecsegm_t prepare_vSegments();
  // second how to intersect those blocks of segments with grid blocks:
  // (override to create special obstacle blocks for local force balances)
  typedef std::vector<std::vector<VolumeSegment_OBB*>> intersect_t;
  virtual intersect_t prepare_segPerBlock(vecsegm_t& vSeg);
  // third how to interpolate on the grid given the intersections:
  virtual void writeSDFOnBlocks(std::vector<VolumeSegment_OBB> & vSegments);

 public:
  Fish(SimulationData&s, cubism::ArgumentParser&p);
  ~Fish() override;
  void save(std::string filename = std::string()) override;
  void restart(std::string filename = std::string()) override;

  virtual void update() override;

  virtual void create() override;
  virtual void finalize() override;

  //MPI stuff, for ObstaclesCreate
  struct BlockID
  {
    //BlockID(double a_h, double a_ox, double a_oy, double a_oz)
    //{
    //  h = a_h;
    //  origin_x = a_ox;
    //  origin_y = a_oy;
    //  origin_z = a_oz;
    //  blockID = 0;
    //}
    double h;
    double origin_x;
    double origin_y;
    double origin_z;
    long long blockID;
  };
  std::vector<BlockID> MyBlockIDs;
  std::vector<std::vector<int>> MySegments;

  struct MPI_Obstacle
  {
    double d [FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX*3 
           + (FluidBlock::sizeZ+2)*(FluidBlock::sizeY+2)*(FluidBlock::sizeX+2)];
    int     i[FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX];
  };
  MPI_Datatype MPI_BLOCKID;
  MPI_Datatype MPI_OBSTACLE;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Fish_h
