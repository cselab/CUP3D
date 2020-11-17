//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include "PoissonSolverAMR.h"

#define PRECOND

#define pVector  tmpU
#define rVector  tmpV
#define xVector  p //current solution estimate
#define sVector  u
#define AxVector v
#define zVector tmpW //preconditioner

namespace cubismup3d {

class PoissonSolverAMR_CG : public PoissonSolverAMR
{
  //will need flux corrections!
  void Get_LHS(bool useX);
 public:
 typedef typename FluidGridMPI::BlockType BlockType;

  void solve() override;

  PoissonSolverAMR_CG(SimulationData& s);

  std::string getName() {
    return "PoissonSolverAMR_CG";
  }

  ~PoissonSolverAMR_CG();

  #ifdef PRECOND
  double getA_local(int I1,int I2);
  void FindZ();
  std::vector<std::vector<double>> Ld;
  std::vector <  std::vector <std::vector< std::pair<int,double> > > >L_row;
  std::vector <  std::vector <std::vector< std::pair<int,double> > > >L_col;
  #endif
};

}//namespace cubismup3d