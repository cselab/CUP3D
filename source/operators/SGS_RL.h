//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_SGS_RL_h
#define CubismUP_3D_SGS_RL_h

//#define SGSRL_STATE_INVARIANTS

#include "operators/Operator.h"

namespace smarties {
class Communicator;
}

CubismUP_3D_NAMESPACE_BEGIN

class SGS_RL : public Operator
{
  smarties::Communicator * const comm;
  const int step;
  const bool timeOut;
  const double reward;
  const double scaleGrads;
  const int nAgentsPerBlock;

public:
  SGS_RL(SimulationData& s, smarties::Communicator* _comm, const int _step,
         const bool _timeOut, const double _reward, const double _scaleGrads,
         const int nAgentsPerBlock);

  ~SGS_RL() { }

  void run() { return (*this)(0.0); /* call next fun with dummy dt */ };
  void operator()(const double dt);

  std::string getName() { return "SGS_RL"; }
};

CubismUP_3D_NAMESPACE_END
#endif
