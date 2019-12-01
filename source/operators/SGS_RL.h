//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_SGS_RL_h
#define CubismUP_3D_SGS_RL_h

#define SGSRL_STATE_INVARIANTS

#include "Operator.h"

namespace smarties {
class Communicator;
}

CubismUP_3D_NAMESPACE_BEGIN

struct HITstatistics;

class SGS_RL : public Operator
{
  smarties::Communicator * const commPtr;
  const int nAgentsPerBlock;
  std::vector<int> agentsIDX;
  std::vector<int> agentsIDY;
  std::vector<int> agentsIDZ;
  std::vector<double> localRewards;

public:
  SGS_RL(SimulationData&s, smarties::Communicator*_comm, const int nAgentsPB);

  void run(const double dt, const bool RLinit, const bool RLover,
           const HITstatistics& stats, const Real collectiveReward);
  void operator()(const double dt) override {}

  std::string getName() override { return "SGS_RL"; }
};

CubismUP_3D_NAMESPACE_END
#endif
