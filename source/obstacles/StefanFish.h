//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Fish.h"

CubismUP_3D_NAMESPACE_BEGIN

class StefanFish: public Fish
{
public:
  StefanFish(SimulationData&s, cubism::ArgumentParser&p);

  //Used for PID controller (which tries to maintain the initial fish position)
  bool bCorrectTrajectory, bCorrectPosition;
  double origC[3];   //initial location
  double origAng = 0;//initial planar angle (in xy plane)

  void create() override;
  void save(std::string filename = std::string()) override;
  void restart(std::string filename) override;

  //Reinforcement Learning functions
  void act(const Real lTact, const std::vector<double>& a) const; //get vector of actions that will be taken at time lTact
  std::vector<double> state() const;                              //return vector of state
  double getLearnTPeriod() const; //take actions once every 0.5*getLearnTPeriod() time
  //// Helpers for state function
  //ssize_t holdingBlockID(const std::array<Real,3> pos, const std::vector<cubism::BlockInfo>& velInfo) const;
  //std::array<int, 3> safeIdInBlock(const std::array<Real,3> pos, const std::array<Real,3> org, const Real invh ) const;
  //std::array<Real, 2> getShear(const std::array<Real,3> pSurf, const std::array<Real,3> normSurf, const std::vector<cubism::BlockInfo>& velInfo) const;
  //#endif
};

CubismUP_3D_NAMESPACE_END
