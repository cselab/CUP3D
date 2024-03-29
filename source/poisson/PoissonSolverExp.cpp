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

double PoissonSolverExp::getA_local(const int i, const int j)
{
  if (i==j)
    return 6.;

  static constexpr int nxny = nx_ * ny_;
  const int ix = i % nx_;
  const int iy = (i / nx_) % ny_;
  const int iz = i / nxny;

  const int jx = j % nx_;
  const int jy = (j / nx_) % ny_;
  const int jz = j / nxny;

  if ((abs(ix-jx) + abs(iy-jy) + abs(iz-jz)) == 1)
    return -1.;

  return 0.;
}

PoissonSolverExp::PoissonSolverExp(SimulationData& s)
  : sim(s), m_comm_(sim.comm), GenericCell(*this), 
    XminCell(*this), XmaxCell(*this), YminCell(*this), YmaxCell(*this), ZminCell(*this), ZmaxCell(*this),
    faceIndexers {&XminCell, &XmaxCell, &YminCell, &YmaxCell, &ZminCell, &ZmaxCell}
{
  // MPI
  MPI_Comm_rank(m_comm_, &rank_);
  MPI_Comm_size(m_comm_, &comm_size_);

  Nblocks_xcumsum_.resize(comm_size_ + 1);
  Nrows_xcumsum_.resize(comm_size_ + 1);

  std::vector<std::vector<double>> L; // lower triangular matrix of Cholesky decomposition
  std::vector<std::vector<double>> L_inv; // inverse of L

  L.resize(nxyz_);
  L_inv.resize(nxyz_);
  for (int i(0); i<nxyz_ ; i++)
  {
    L[i].resize(i+1);
    L_inv[i].resize(i+1);
    // L_inv will act as right block in GJ algorithm, init as identity
    for (int j(0); j<=i; j++){
      L_inv[i][j] = (i == j) ? 1. : 0.;
    }
  }

  // compute the Cholesky decomposition of the preconditioner with Cholesky-Crout
  for (int i(0); i<nxyz_ ; i++)
  {
    double s1 = 0;
    for (int k(0); k<=i-1; k++)
      s1 += L[i][k]*L[i][k];
    L[i][i] = sqrt(getA_local(i,i) - s1);
    for (int j(i+1); j<nxyz_; j++)
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
  for (int br(0); br<nxyz_; br++)
  { // 'br' - base row in which all columns up to L_lb[br][br] are already zero
    const double bsf = 1. / L[br][br];
    for (int c(0); c<=br; c++)
      L_inv[br][c] *= bsf;

    for (int wr(br+1); wr<nxyz_; wr++)
    { // 'wr' - working row where elements below L_lb[br][br] will be set to zero
      const double wsf = L[wr][br];
      for (int c(0); c<=br; c++)
        L_inv[wr][c] -= (wsf * L_inv[br][c]);
    }
  }

  // P_inv_ holds inverse preconditionner in row major order!
  std::vector<double> P_inv(nxyz_ * nxyz_);
  for (int i(0); i<nxyz_; i++)
  for (int j(0); j<nxyz_; j++)
  {
    double aux = 0.;
    for (int k(0); k<nxyz_; k++) // P_inv_ = (L^T)^{-1} L^{-1}
      aux += (i <= k && j <=k) ? L_inv[k][i] * L_inv[k][j] : 0.;

    P_inv[i*nxyz_+j] = -aux; // Up to now Cholesky of negative P to avoid complex numbers
  }

  // Create Linear system and backend solver objects
  LocalLS_ = std::make_unique<LocalSpMatDnVec>(m_comm_, nxyz_, sim.bMeanConstraint, P_inv);
}
void PoissonSolverExp::interpolate(
    const cubism::BlockInfo &info_c, const int ix_c, const int iy_c, const int iz_c,
    const cubism::BlockInfo &info_f, const long long fine_close_idx, const long long fine_far_idx,
    const double sign, const double sign_ds1, const double sign_ds2,
    const FaceCellIndexer &indexer, SpRowInfo& row) const
{
  const int rank_c = sim.lhs->Tree(info_c).rank();
  const int rank_f = sim.lhs->Tree(info_f).rank();
  const double h = info_f.h;

  // 2./3.*p_fine_close_idx - 1./5.*p_fine_far_idx
  row.mapColVal(rank_f, fine_close_idx, sign * h * 2./3.);
  row.mapColVal(rank_f, fine_far_idx,  -sign * h * 1./5.);

  // 8./15 * p_T, constant term
  const double tf = sign * h * 8./15.; // common factor for all terms of Taylor expansion
  row.mapColVal(rank_c, indexer.This(info_c, ix_c, iy_c, iz_c), tf);

  // first derivative terms
  std::array<std::pair<long long, double>, 3> D;

  D = D1(info_c, indexer, ix_c, iy_c, iz_c, 1);
  for (int i(0); i < 3; i++)
    row.mapColVal(rank_c, D[i].first, sign_ds1 * tf * D[i].second);

  D = D1(info_c, indexer, ix_c, iy_c, iz_c, 2);
  for (int i(0); i < 3; i++)
    row.mapColVal(rank_c, D[i].first, sign_ds2 * tf * D[i].second);
  
  // second derivative terms
  D = D2(info_c, indexer, ix_c, iy_c, iz_c, 1);
  for (int i(0); i < 3; i++)
    row.mapColVal(rank_c, D[i].first, tf * D[i].second);

  D = D2(info_c, indexer, ix_c, iy_c, iz_c, 2);
  for (int i(0); i < 3; i++)
    row.mapColVal(rank_c, D[i].first, tf * D[i].second);

  // mixed derivative terms
  std::array<std::pair<long long, double>, 4> Dm = Dmixed(info_c, indexer, ix_c, iy_c, iz_c);
  for (int i(0); i < 4; i++)
    row.mapColVal(rank_c, Dm[i].first, sign_ds1 * sign_ds2 * tf * Dm[i].second);

}
 
void PoissonSolverExp::makeFlux(
    const BlockInfo& rhs_info,
    const int ix, const int iy, const int iz,
    const bool &isBoundary,
    const BlockInfo &rhsNei,
    const FaceCellIndexer& indexer,
    SpRowInfo &row) const
{
  const long long sfc_idx = indexer.This(rhs_info, ix, iy, iz);

  if (this->sim.lhs->Tree(rhsNei).Exists())
  { 
    const int nei_rank = sim.lhs->Tree(rhsNei).rank();
    const long long nei_idx = indexer.neiUnif(rhsNei, ix, iy, iz);
    const double h = rhsNei.h;

    // Map flux associated to out-of-block edges at the same level of refinement
    row.mapColVal(nei_rank, nei_idx, h);
    row.mapColVal(sfc_idx, -h);
  }
  else if (this->sim.lhs->Tree(rhsNei).CheckCoarser())
  {
    const BlockInfo &rhsNei_c = sim.lhs->getBlockInfoAll(rhs_info.level - 1 , rhsNei.Zparent);
    const int ix_c = indexer.ix_c(rhs_info, ix);
    const int iy_c = indexer.iy_c(rhs_info, iy);
    const int iz_c = indexer.iz_c(rhs_info, iz);
    const long long inward_idx = indexer.neiInward(rhs_info, ix, iy, iz);
    const double sign_ds1 = indexer.taylorSign(ix, iy, iz, 1);
    const double sign_ds2 = indexer.taylorSign(ix, iy, iz, 2);
    const double h = rhs_info.h;
      
    interpolate(rhsNei_c, ix_c, iy_c, iz_c, rhs_info, sfc_idx, inward_idx, 1., sign_ds1, sign_ds2, indexer, row);
    row.mapColVal(sfc_idx, -h);
  }
  else if (this->sim.lhs->Tree(rhsNei).CheckFiner())
  {
    const BlockInfo &rhsNei_f = this->sim.lhs->getBlockInfoAll(rhs_info.level + 1, indexer.Zchild(rhsNei, ix, iy, iz));
    const int nei_rank = sim.lhs->Tree(rhsNei_f).rank();
    const double h = rhsNei_f.h;
  
    // F1
    long long fine_close_idx = indexer.neiFine1(rhsNei_f, ix, iy, iz, 0);
    long long fine_far_idx   = indexer.neiFine1(rhsNei_f, ix, iy, iz, 1);
    row.mapColVal(nei_rank, fine_close_idx, h);
    interpolate(rhs_info, ix, iy, iz, rhsNei_f, fine_close_idx, fine_far_idx, -1., -1., -1., indexer, row);
    // F2
    fine_close_idx = indexer.neiFine2(rhsNei_f, ix, iy, iz, 0);
    fine_far_idx   = indexer.neiFine2(rhsNei_f, ix, iy, iz, 1);
    row.mapColVal(nei_rank, fine_close_idx, h);
    interpolate(rhs_info, ix, iy, iz, rhsNei_f, fine_close_idx, fine_far_idx, -1., 1., -1., indexer, row);
    // F3
    fine_close_idx = indexer.neiFine3(rhsNei_f, ix, iy, iz, 0);
    fine_far_idx   = indexer.neiFine3(rhsNei_f, ix, iy, iz, 1);
    row.mapColVal(nei_rank, fine_close_idx, h);
    interpolate(rhs_info, ix, iy, iz, rhsNei_f, fine_close_idx, fine_far_idx, -1., -1., 1., indexer, row);
    // F4
    fine_close_idx = indexer.neiFine4(rhsNei_f, ix, iy, iz, 0);
    fine_far_idx   = indexer.neiFine4(rhsNei_f, ix, iy, iz, 1);
    row.mapColVal(nei_rank, fine_close_idx, h);
    interpolate(rhs_info, ix, iy, iz, rhsNei_f, fine_close_idx, fine_far_idx, -1., 1., 1., indexer, row);
  }
  else { throw std::runtime_error("Neighbour doesn't exist, isn't coarser, nor finer..."); }
}

void PoissonSolverExp::getMat()
{
  //sim.startProfiler("Make LHS");

  // Update blockID's for blocks from other ranks
  sim.lhs->UpdateBlockInfoAll_States(true); 
  // This returns an array with the blocks that the coarsest possible 
  // mesh would have (i.e. all blocks are at level 0)
  std::array<int, 3> blocksPerDim = sim.lhs->getMaxBlocks();

  //Get a vector of all BlockInfos of the grid we're interested in
  std::vector<cubism::BlockInfo>&  RhsInfo = sim.lhsInfo();
  const int Nblocks = RhsInfo.size();
  const int N = Nblocks * nxyz_;

  // Reserve sufficient memory for LS proper to the rank
  LocalLS_->reserve(N);

  // Calculate cumulative sums for blocks and rows for correct global indexing
  const long long Nblocks_long = Nblocks;
  MPI_Allgather(&Nblocks_long, 1, MPI_LONG_LONG, Nblocks_xcumsum_.data(), 1, MPI_LONG_LONG, m_comm_);
  for (int i(Nblocks_xcumsum_.size()-1); i > 0; i--)
    Nblocks_xcumsum_[i] = Nblocks_xcumsum_[i-1]; // shift to right for rank 'i+1' to have cumsum of rank 'i'
  
  // Set cumsum for rank 0 to zero
  Nblocks_xcumsum_[0] = 0;
  Nrows_xcumsum_[0] = 0;

  // Perform cumulative sum
  for (size_t i(1); i < Nblocks_xcumsum_.size(); i++)
  {
    Nblocks_xcumsum_[i] += Nblocks_xcumsum_[i-1];
    Nrows_xcumsum_[i] = Nblocks_xcumsum_[i] * nxyz_;
  }

  // No parallel for to ensure COO are ordered at construction
  for(int i=0; i<Nblocks; i++)
  {    
    const BlockInfo &rhs_info = RhsInfo[i];

    //1.Check if this is a boundary block
    const int aux = 1 << rhs_info.level; // = 2^level
    const int MAX_X_BLOCKS = blocksPerDim[0]*aux - 1; //this means that if level 0 has blocksPerDim[0] blocks in the x-direction, level rhs.level will have this many blocks
    const int MAX_Y_BLOCKS = blocksPerDim[1]*aux - 1;
    const int MAX_Z_BLOCKS = blocksPerDim[2]*aux - 1;

    std::array<bool, 6> isBoundary;
    isBoundary[0] = (rhs_info.index[0] == 0           ); // Xm, same order as faceIndexers made in constructor!
    isBoundary[1] = (rhs_info.index[0] == MAX_X_BLOCKS); // Xp
    isBoundary[2] = (rhs_info.index[1] == 0           ); // Ym
    isBoundary[3] = (rhs_info.index[1] == MAX_Y_BLOCKS); // Yp
    isBoundary[4] = (rhs_info.index[2] == 0           ); // Zm
    isBoundary[5] = (rhs_info.index[2] == MAX_Z_BLOCKS); // Zp

    std::array<bool, 3> isPeriodic; // same dimension ordering as isBoundary
    isPeriodic[0] = (sim.BCx_flag == periodic);
    isPeriodic[1] = (sim.BCy_flag == periodic);
    isPeriodic[2] = (sim.BCz_flag == periodic);

    //2.Access the block's neighbors (for the Poisson solve in two dimensions we care about four neighbors in total)
    std::array<long long, 6> Z;
    Z[0] = rhs_info.Znei[1-1][1][1]; // Xm
    Z[1] = rhs_info.Znei[1+1][1][1]; // Xp
    Z[2] = rhs_info.Znei[1][1-1][1]; // Ym
    Z[3] = rhs_info.Znei[1][1+1][1]; // Yp
    Z[4] = rhs_info.Znei[1][1][1-1]; // Zm
    Z[5] = rhs_info.Znei[1][1][1+1]; // Zp
    //rhs.Z == rhs.Znei[1][1][1] is true always

    std::array<const BlockInfo*, 6> rhsNei;
    rhsNei[0] = &(this->sim.lhs->getBlockInfoAll(rhs_info.level, Z[0]));
    rhsNei[1] = &(this->sim.lhs->getBlockInfoAll(rhs_info.level, Z[1]));
    rhsNei[2] = &(this->sim.lhs->getBlockInfoAll(rhs_info.level, Z[2]));
    rhsNei[3] = &(this->sim.lhs->getBlockInfoAll(rhs_info.level, Z[3]));
    rhsNei[4] = &(this->sim.lhs->getBlockInfoAll(rhs_info.level, Z[4]));
    rhsNei[5] = &(this->sim.lhs->getBlockInfoAll(rhs_info.level, Z[5]));

    // Record local index of row which is to be modified with bMeanConstraint == 1
    if (sim.bMeanConstraint == 1 &&
        rhs_info.index[0] == 0 && 
        rhs_info.index[1] == 0 && 
        rhs_info.index[2] == 0)
      LocalLS_->set_bMeanRow(GenericCell.This(rhs_info, 0, 0, 0) - Nrows_xcumsum_[rank_]);

    for (int iz(0); iz<nz_; iz++)
    for (int iy(0); iy<ny_; iy++)
    for (int ix(0); ix<nx_; ix++)
    { // Logic needs to be in 'for' loop to consruct cooRows in order
      const long long sfc_idx = GenericCell.This(rhs_info, ix, iy, iz);  
      const double h = rhs_info.h;
      if (sim.bMeanConstraint == 3 &&
          rhs_info.index[0] == 0 && 
          rhs_info.index[1] == 0 && 
          rhs_info.index[2] == 0 &&
          ix == 0 && iy == 0 && iz == 0)
      {
        LocalLS_->cooPushBackVal(h, sfc_idx, sfc_idx);
      }
      else if ((ix > 0 && ix<nx_-1) && (iy > 0 && iy<ny_-1) && (iz > 0 && iz<nz_-1))
      { // Inner cells

        // Push back in ascending order for column index
        LocalLS_->cooPushBackVal(h, sfc_idx, GenericCell.This(rhs_info, ix, iy, iz-1));
        LocalLS_->cooPushBackVal(h, sfc_idx, GenericCell.This(rhs_info, ix, iy-1, iz));
        LocalLS_->cooPushBackVal(h, sfc_idx, GenericCell.This(rhs_info, ix-1, iy, iz));
        LocalLS_->cooPushBackVal(-6.*h, sfc_idx, sfc_idx);
        LocalLS_->cooPushBackVal(h, sfc_idx, GenericCell.This(rhs_info, ix+1, iy, iz));
        LocalLS_->cooPushBackVal(h, sfc_idx, GenericCell.This(rhs_info, ix, iy+1, iz));
        LocalLS_->cooPushBackVal(h, sfc_idx, GenericCell.This(rhs_info, ix, iy, iz+1));
      }
      else
      {
        // See if face is shared with a cell from different block
        std::array<bool, 6> validNei;
        validNei[0] = GenericCell.validXm(ix, iy, iz); 
        validNei[1] = GenericCell.validXp(ix, iy, iz); 
        validNei[2] = GenericCell.validYm(ix, iy, iz); 
        validNei[3] = GenericCell.validYp(ix, iy, iz); 
        validNei[4] = GenericCell.validZm(ix, iy, iz); 
        validNei[5] = GenericCell.validZp(ix, iy, iz); 

        // Get index of cell accross the face
        std::array<long long, 6> idxNei;
        idxNei[0] = GenericCell.This(rhs_info, ix-1, iy, iz);
        idxNei[1] = GenericCell.This(rhs_info, ix+1, iy, iz);
        idxNei[2] = GenericCell.This(rhs_info, ix, iy-1, iz);
        idxNei[3] = GenericCell.This(rhs_info, ix, iy+1, iz);
        idxNei[4] = GenericCell.This(rhs_info, ix, iy, iz-1);
        idxNei[5] = GenericCell.This(rhs_info, ix, iy, iz+1);

        SpRowInfo row(sim.lhs->Tree(rhs_info).rank(), sfc_idx, 16);

        for (int j(0); j < 6; j++)
        { // Iterate over each face of cell
          if (validNei[j])
          { // This face is 'inner' wrt to the block
            row.mapColVal(idxNei[j], h);
            row.mapColVal(sfc_idx, -h);  // diagonal element
          }
          // Outer face, do nothing if non-periodic BC
          else if (!isBoundary[j] || (isBoundary[j] && isPeriodic[j/2]))
            this->makeFlux(rhs_info, ix, iy, iz, isBoundary[j], *rhsNei[j], *faceIndexers[j], row);
        }

        LocalLS_->cooPushBackRow(row);
      } // if else
    } // for ix iy iz
  } // for(int i=0; i< Nblocks; i++)

  LocalLS_->make(Nrows_xcumsum_);

  //sim.stopProfiler();
}

void PoissonSolverExp::getVec()
{
  //Get a vector of all BlockInfos of the grid we're interested in
  std::vector<cubism::BlockInfo>&  RhsInfo = sim.lhsInfo();
  std::vector<cubism::BlockInfo>&  zInfo = sim.presInfo();
  const int Nblocks = RhsInfo.size();

  std::vector<double>& x = LocalLS_->get_x();
  std::vector<double>& b = LocalLS_->get_b();
  std::vector<double>& h3 = LocalLS_->get_h3();
  std::vector<double>& invh = LocalLS_->get_invh();
  const long long shift = -Nrows_xcumsum_[rank_];

  // Copy RHS and LHS vec initial guess, if LS was updated getMat reallocates sufficient memory
  #pragma omp parallel for
  for (int i=0; i<Nblocks; i++)
  {
    const ScalarBlock & __restrict__ rhs = *(ScalarBlock*)RhsInfo[i].ptrBlock;
    const ScalarBlock & __restrict__ p = *(ScalarBlock*)zInfo[i].ptrBlock;

    h3[i]   = RhsInfo[i].h * RhsInfo[i].h * RhsInfo[i].h;
    invh[i] = 1. / RhsInfo[i].h;
    for (int iz(0); iz<nz_; iz++)
    for (int iy(0); iy<ny_; iy++)
    for (int ix(0); ix<nx_; ix++)
    {
      const long long sfc_loc = GenericCell.This(RhsInfo[i], ix, iy, iz) + shift;
      if ((sim.bMeanConstraint == 1 || sim.bMeanConstraint == 3) &&
           RhsInfo[i].index[0] == 0 && 
           RhsInfo[i].index[1] == 0 && 
           RhsInfo[i].index[2] == 0 &&
           ix == 0 && iy == 0 && iz == 0)
      {
        b[sfc_loc] = 0.;
        x[sfc_loc] = 0.;
      }
      else
      {
        b[sfc_loc] = rhs(ix, iy, iz).s;
        x[sfc_loc] = p(ix, iy, iz).s;
      }
    }
  }
}

void PoissonSolverExp::solve() 
{
  if (rank_ == 0)
  {
    if (sim.verbose)
      std::cout << "--------------------- Calling on ExpAMRSolver.solve() ------------------------ \n";
    else
      std::cout << '\n';
  }

  //const double max_error = this->sim.step < 10 ? 0.0 : sim.PoissonErrorTol;
  //const double max_rel_error = this->sim.step < 10 ? 0.0 : sim.PoissonErrorTolRel;
  //const int max_restarts = this->sim.step < 10 ? 100 : sim.maxPoissonRestarts;
  const double max_error = sim.PoissonErrorTol;
  const double max_rel_error = sim.PoissonErrorTolRel;
  const int max_restarts = 100;

  if (sim.pres->UpdateFluxCorrection)
  {
    sim.pres->UpdateFluxCorrection = false;
    this->getMat();
    this->getVec();
    LocalLS_->solveWithUpdate(max_error, max_rel_error, max_restarts);
  }
  else
  {
    this->getVec();
    LocalLS_->solveNoUpdate(max_error, max_rel_error, max_restarts);
  }

  //Now that we found the solution, we just substract the mean to get a zero-mean solution. 
  //This can be done because the solver only cares about grad(P) = grad(P-mean(P))
  std::vector<cubism::BlockInfo>&  zInfo = sim.presInfo();
  const int Nblocks = zInfo.size();

  const std::vector<double>& x = LocalLS_->get_x();
  const long long shift = -Nrows_xcumsum_[rank_];

  #pragma omp parallel for
  for (int i=0; i < Nblocks; i++)
  {
    const int m = zInfo[i].level;
    const long long n = zInfo[i].Z;
    const BlockInfo & info = sim.pres->getBlockInfoAll(m,n);
    ScalarBlock & __restrict__ b  = *(ScalarBlock*) info.ptrBlock;
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++)
    {
      const long long sfc_loc = GenericCell.This(zInfo[i], ix, iy, iz) + shift;
      b(ix,iy,iz).s = x[sfc_loc];
    }
  }
}
}//namespace cubismup3d
