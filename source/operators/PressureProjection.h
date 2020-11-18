//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_PressureProjection_h
#define CubismUP_3D_PressureProjection_h

#include "Operator.h"
#include "../poisson/PoissonSolverAMR.h"


CubismUP_3D_NAMESPACE_BEGIN

class PressureProjection : public Operator
{
 protected:
  PoissonSolverAMR * pressureSolver;

 public:
  PressureProjection(SimulationData & s); 
  ~PressureProjection()
  {
    if (sim.pressureSolver != nullptr)
    {
        delete sim.pressureSolver;
        sim.pressureSolver = nullptr;
    }
  };

  void operator()(const double dt);

  std::string getName() { return "PressureProjection"; }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_PressureProjection_h
