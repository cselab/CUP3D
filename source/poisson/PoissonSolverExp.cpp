//
//  CubismUP_3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  
//

#include <iostream>
#include <vector>

#include "PoissonSolverExp.h"
#include "BiCGSTAB.cuh"

namespace cubismup3d {

double PoissonSolverExp::getA_local(const int& i, const int& j)
{
  if (i==j)
    return 6.;

  static constexpr int nxny = nx_ * ny_;
  const int ix = i % nx_;
  const int iy = i / nx_;
  const int iz = i / nxny;

  const int jx = j % nx_;
  const int jy = j / nx_;
  const int jz = j / nxny;

  if ((abs(ix-jx) + abs(iy-jy) + abs(iz-jz)) == 1)
    return -1.;

  return 0.;
}

PoissonSolverExp::PoissonSolverExp(SimulationData& s)
  : sim(s), m_comm_(sim.app_comm)
{
  // MPI
  MPI_Comm_rank(m_comm_, &rank_);
  MPI_Comm_size(m_comm_, &comm_size_);

//  Nblocks_xcumsum_.resize(comm_size_ + 1);
//  Nrows_xcumsum_.resize(comm_size_ + 1);

  std::vector<std::vector<double>> L; // lower triangular matrix of Cholesky decomposition
  std::vector<std::vector<double>> L_inv; // inverse of L

  static constexpr int nlen = nx_ * ny_ * nz_;
  L.resize(nlen);
  L_inv.resize(nlen);
  for (int i(0); i<nlen ; i++)
  {
    L[i].resize(i+1);
    L_inv[i].resize(i+1);
    // L_inv will act as right block in GJ algorithm, init as identity
    for (int j(0); j<=i; j++){
      L_inv[i][j] = (i == j) ? 1. : 0.;
    }
  }

  // compute the Cholesky decomposition of the preconditioner with Cholesky-Crout
  for (int i(0); i<nlen ; i++)
  {
    double s1 = 0;
    for (int k(0); k<=i-1; k++)
      s1 += L[i][k]*L[i][k];
    L[i][i] = sqrt(getA_local(i,i) - s1);
    for (int j(i+1); j<nlen; j++)
    {
      double s2 = 0;
      for (int k(0); k<=i-1; k++)
        s2 += L[i][k]*L[j][k];
      L[j][i] = (getA_local(j,i)-s2) / L[i][i];
    }
  }

  /* Compute the inverse of the Cholesky decomposition L using Gauss-Jordan elimination.
     L will act as the left block (it does not need to be modified in the process), 
     L_inv will act as the right block and at the end of the algo will contain the inverse */
  for (int br(0); br<nlen; br++)
  { // 'br' - base row in which all columns up to L_lb[br][br] are already zero
    const double bsf = 1. / L[br][br];
    for (int c(0); c<=br; c++)
      L_inv[br][c] *= bsf;

    for (int wr(br+1); wr<nlen; wr++)
    { // 'wr' - working row where elements below L_lb[br][br] will be set to zero
      const double wsf = L[wr][br];
      for (int c(0); c<=br; c++)
        L_inv[wr][c] -= (wsf * L_inv[br][c]);
    }
  }

  // P_inv_ holds inverse preconditionner in row major order!
  std::vector<double> P_inv(nlen * nlen);
  for (int i(0); i<nlen; i++)
  for (int j(0); j<nlen; j++)
  {
    double aux = 0.;
    for (int k(0); k<nlen; k++) // P_inv_ = (L^T)^{-1} L^{-1}
      aux += (i <= k && j <=k) ? L_inv[k][i] * L_inv[k][j] : 0.;

    P_inv[i*nlen+j] = -aux; // Up to now Cholesky of negative P to avoid complex numbers
  }

  // Create Linear system and backend solver objects
  LocalLS_ = std::make_unique<LocalSpMatDnVec>(m_comm_, nlen, P_inv);
}
 
void PoissonSolverExp::solve() 
{
  std::cerr << "Hello PoissonSolverExp!\n";
  throw;
}
}//namespace cubismup3d
