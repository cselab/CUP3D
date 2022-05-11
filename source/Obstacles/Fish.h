//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

class FishMidlineData;
struct VolumeSegment_OBB;

class Fish: public Obstacle
{
 protected:
  void integrateMidline();

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
  virtual void create() override;
  FishMidlineData * myFish = nullptr;

  struct BlockID
  {
    Real h;
    Real origin_x;
    Real origin_y;
    Real origin_z;
    long long blockID;
  };
  std::vector<BlockID> MyBlockIDs;
  std::vector<std::vector<int>> MySegments;

  #if 1
  //MPI stuff, for ObstaclesCreate
  struct MPI_Obstacle
  {
    Real d [FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX*3 
           + (FluidBlock::sizeZ+2)*(FluidBlock::sizeY+2)*(FluidBlock::sizeX+2)];
    int     i[FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX];
  };
  MPI_Datatype MPI_BLOCKID;
  MPI_Datatype MPI_OBSTACLE;
  #endif
};

CubismUP_3D_NAMESPACE_END
