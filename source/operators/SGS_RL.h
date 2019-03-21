//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_AdvectionDiffusion_h
#define CubismUP_3D_AdvectionDiffusion_h

#include "operators/Operator.h"

class Communicator;

CubismUP_3D_NAMESPACE_BEGIN

class SGS_RL : public Operator
{
  Communicator * const comm;
  const bool timeOut;
  const double reward;

public:
  SGS_RL(SimulationData& s, Communicator* c, bool timeOut, double rew);

  ~SGS_RL() { }

  void run() { return (*this)(0.0); /* call next fun with dummy dt */ };
  void operator()(const double dt);

  std::string getName() { return "SGS_RL"; }
};

CubismUP_3D_NAMESPACE_END
#endif
