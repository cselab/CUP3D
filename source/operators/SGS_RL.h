//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_SGS_RL_h
#define CubismUP_3D_SGS_RL_h

#include "operators/Operator.h"

class Communicator;

CubismUP_3D_NAMESPACE_BEGIN

class SGS_RL : public Operator
{
  Communicator * const comm;
  const int step;
  const bool timeOut;
  const bool evalStep;
  const double reward;
  const int nAgentsPerBlock;

public:
  SGS_RL(SimulationData& s, Communicator* _comm, const int _step,
         const bool _timeOut,const bool _evalStep, const double _reward, const int nAgentsPerBlock);

  ~SGS_RL() { }

  void run() { return (*this)(0.0); /* call next fun with dummy dt */ };
  void operator()(const double dt);

  template <typename Kernel>
  void runKernel(const Kernel& kernel);

  std::string getName() { return "SGS_RL"; }
};

CubismUP_3D_NAMESPACE_END
#endif
