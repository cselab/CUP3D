//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifndef CubismUP_3D_SpectralManipACCFFT_h
#define CubismUP_3D_SpectralManipACCFFT_h

#include "SpectralManip.h"

CubismUP_3D_NAMESPACE_BEGIN

class SpectralManipACC : public SpectralManip
{
  MPI_Comm acc_comm;
  // the local pencil size and the allocation size
  int isize[3], osize[3], istart[3], ostart[3];
  size_t alloc_max;
  Real * phi_hat;
  void * plan;

public:

  SpectralManipACC(SimulationData & s);
  ~SpectralManipACC();

  void prepareFwd() override;
  void prepareBwd() override;

  void runFwd() const override;
  void runBwd() const override;

  void _compute_largeModesForcing() override;
  void _compute_analysis() override;
  void _compute_IC(const std::vector<Real> &K,
                   const std::vector<Real> &E) override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralManipACCFFT_h
