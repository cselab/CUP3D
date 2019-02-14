//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef __IncompressibleFluids3D__IF3D_NacaOperator__
#define __IncompressibleFluids3D__IF3D_NacaOperator__

#include "IF3D_FishOperator.h"

class IF3D_NacaOperator: public IF3D_FishOperator
{
  double Apitch, Fpitch, Ppitch, Mpitch, Fheave, Aheave, tOld = 0;
  bool bCreated;
 public:
  IF3D_NacaOperator(SimulationData&s, ArgumentParser&p);
  void update(const int stepID, const double t, const double dt, const Real* Uinf) override;
  void computeVelocities(const Real* Uinf) override;
  void writeSDFOnBlocks(const mapBlock2Segs& segmentsPerBlock) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
