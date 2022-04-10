//
//  CubismUP_3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  
//

#pragma once

#include <memory>

#include <Cubism/BlockInfo.h>
#include "../Definitions.h"
#include "../SimulationData.h"
#include "PoissonSolverBase.h"
#include "LocalSpMatDnVec.h"

namespace cubismup3d {

class PoissonSolverExp : public PoissonSolverBase
{
 public:
  PoissonSolverExp(SimulationData&s);
  PoissonSolverExp(const PoissonSolverExp& c) = delete; 
  void solve() override;

 protected:
  SimulationData & sim;
  
  int rank_;
  MPI_Comm m_comm_;
  int comm_size_;

  static constexpr int nx_ = BlockType::sizeX;
  static constexpr int ny_ = BlockType::sizeY;
  static constexpr int nz_ = BlockType::sizeZ;

  // Returns element of preconditioner negative K_{i,j}
  double getA_local(const int& i, const int& j);

  // Distributed linear system with local indexing
  std::unique_ptr<LocalSpMatDnVec> LocalLS_;

};

}//namespace cubismup3d
