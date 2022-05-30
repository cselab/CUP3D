#include <unordered_map>
#include <algorithm> // std::copy
#include <iostream>

#include "LocalSpMatDnVec.h"
#include "BiCGSTAB.cuh"

LocalSpMatDnVec::LocalSpMatDnVec(MPI_Comm m_comm, const int BLEN, const int bMeanConstraint, const std::vector<double>& P_inv) 
  : m_comm_(m_comm), BLEN_(BLEN)
{
  // MPI
  MPI_Comm_rank(m_comm_, &rank_);
  MPI_Comm_size(m_comm_, &comm_size_);

  bd_recv_set_.resize(comm_size_);
  bd_recv_vec_.resize(comm_size_);
  recv_ranks_.reserve(comm_size_); 
  recv_offset_.reserve(comm_size_); 
  recv_sz_.reserve(comm_size_);
  send_ranks_.reserve(comm_size_); 
  send_offset_.reserve(comm_size_); 
  send_sz_.reserve(comm_size_); 

  solver_ = std::make_unique<BiCGSTABSolver>(m_comm, *this, BLEN, bMeanConstraint, P_inv);
}

LocalSpMatDnVec::~LocalSpMatDnVec() {}

void LocalSpMatDnVec::reserve(const int N)
{
  m_ = N;
  bMeanRow_ = -1; // default value

  // Clear previous contents and reserve excess memory
  for (size_t i(0); i < bd_recv_set_.size(); i++)
    bd_recv_set_[i].clear();

  loc_cooValA_.clear(); loc_cooValA_.reserve(10*N);
  loc_cooRowA_long_.clear(); loc_cooRowA_long_.reserve(10*N);
  loc_cooColA_long_.clear(); loc_cooColA_long_.reserve(10*N);
  bd_cooValA_.clear(); bd_cooValA_.reserve(N);
  bd_cooRowA_long_.clear(); bd_cooRowA_long_.reserve(N);
  bd_cooColA_long_.clear(); bd_cooColA_long_.reserve(N);

  x_.resize(N);
  b_.resize(N);
  h3_.resize(N/BLEN_);    // grid-spacing per block
  invh_.resize(N/BLEN_); 
}

void LocalSpMatDnVec::cooPushBackVal(const double val, const long long row, const long long col)
{
  loc_cooValA_.push_back(val);  
  loc_cooRowA_long_.push_back(row);
  loc_cooColA_long_.push_back(col);
}

void LocalSpMatDnVec::cooPushBackRow(const SpRowInfo &row)
{
  for (const auto &[col_idx, val] : row.loc_colval_)
  {
    loc_cooValA_.push_back(val);  
    loc_cooRowA_long_.push_back(row.idx_);
    loc_cooColA_long_.push_back(col_idx);
  }
  if (!row.neirank_cols_.empty())
  { 
    for (const auto &[col_idx, val] : row.bd_colval_)
    {
      bd_cooValA_.push_back(val);  
      bd_cooRowA_long_.push_back(row.idx_);
      bd_cooColA_long_.push_back(col_idx);
    }
    // Update recv set
    for (const auto &[rank, col_idx] : row.neirank_cols_)
    {
      bd_recv_set_[rank].insert(col_idx);
    }
  }
}

void LocalSpMatDnVec::make(const std::vector<long long> &Nrows_xcumsum)
{
  loc_nnz_ = loc_cooValA_.size();
  bd_nnz_  = bd_cooValA_.size();
  
  halo_ = 0;

  std::vector<int> send_sz_allranks(comm_size_);
  std::vector<int> recv_sz_allranks(comm_size_);

  for (int r(0); r < comm_size_; r++)
    recv_sz_allranks[r] = bd_recv_set_[r].size();

  // Exchange message sizes between all ranks
  MPI_Alltoall(recv_sz_allranks.data(), 1, MPI_INT, send_sz_allranks.data(), 1, MPI_INT, m_comm_);

  // Set receiving rules into halo
  recv_ranks_.clear();
  recv_offset_.clear();
  recv_sz_.clear();
  int offset = 0;
  for (int r(0); r < comm_size_; r++)
  {
    if (r != rank_ && recv_sz_allranks[r] > 0)
    {
      recv_ranks_.push_back(r); 
      recv_offset_.push_back(offset);
      recv_sz_.push_back(recv_sz_allranks[r]);
      offset += recv_sz_allranks[r];
    }
  }
  halo_ = offset;

  // Set sending rules from a 'send' buffer
  send_ranks_.clear();
  send_offset_.clear();
  send_sz_.clear();
  offset = 0;
  for (int r(0); r < comm_size_; r++)
  {
    if (r !=rank_ && send_sz_allranks[r] > 0)
    {
      send_ranks_.push_back(r);
      send_offset_.push_back(offset);
      send_sz_.push_back(send_sz_allranks[r]);
      offset += send_sz_allranks[r];
    }
  }
  std::vector<long long> send_pack_idx_long(offset);
  send_pack_idx_.resize(offset);


  // Post receives for column indices from other ranks required for SpMV
  std::vector<MPI_Request> recv_requests(send_ranks_.size());
  for (size_t i(0); i < send_ranks_.size(); i++)
    MPI_Irecv(&send_pack_idx_long[send_offset_[i]], send_sz_[i], MPI_LONG_LONG, send_ranks_[i], 546, m_comm_, &recv_requests[i]);


  // Create sends to inform other ranks what they will need to send here
  std::vector<long long> recv_idx_list(halo_);
  std::vector<MPI_Request> send_requests(recv_ranks_.size());
  for (size_t i(0); i < recv_ranks_.size(); i++)
  {
    std::copy(bd_recv_set_[recv_ranks_[i]].begin(), bd_recv_set_[recv_ranks_[i]].end(), &recv_idx_list[recv_offset_[i]]);
    MPI_Isend(&recv_idx_list[recv_offset_[i]], recv_sz_[i], MPI_LONG_LONG, recv_ranks_[i], 546, m_comm_, &send_requests[i]);
  }

  // Now re-index the linear system from global to local indexing
  const long long shift = -Nrows_xcumsum[rank_];
  loc_cooRowA_int_.resize(loc_nnz_);
  loc_cooColA_int_.resize(loc_nnz_);
  bd_cooRowA_int_.resize(bd_nnz_);
  #pragma omp parallel
  {
    // Shift rows and columns to local indexing
    #pragma omp for
    for (int i=0; i < loc_nnz_; i++)
      loc_cooRowA_int_[i] = (int)(loc_cooRowA_long_[i] + shift);
    #pragma omp for
    for (int i=0; i < loc_nnz_; i++)
      loc_cooColA_int_[i] = (int)(loc_cooColA_long_[i] + shift);
    #pragma omp for
    for (int i=0; i < bd_nnz_; i++)
      bd_cooRowA_int_[i] = (int)(bd_cooRowA_long_[i] + shift);
  }

  // Make sure recv_idx_list is safe to use (and not prone to deallocation
  // because it will no longer be in use)
  MPI_Waitall(send_ranks_.size(), recv_requests.data(), MPI_STATUS_IGNORE);

  // Map indices of columns from other ranks to the halo
  std::unordered_map<long long, int> bd_reindex_map; 
  bd_reindex_map.reserve(halo_);
  for (int i(0); i < halo_; i++)
    bd_reindex_map[recv_idx_list[i]] = m_ + i;

  bd_cooColA_int_.resize(bd_nnz_);
  for (int i=0; i < bd_nnz_; i++)
    bd_cooColA_int_[i] = bd_reindex_map[bd_cooColA_long_[i]];

  MPI_Waitall(recv_ranks_.size(), send_requests.data(), MPI_STATUS_IGNORE);

  #pragma omp parallel for
  for (size_t i=0; i < send_pack_idx_.size(); i++)
    send_pack_idx_[i] = (int)(send_pack_idx_long[i] + shift);

//  if (rank_ == 0)
//    std::cerr << "  [LocalLS]: Rank: " << rank_ << ", m: " << m_ << ", halo: " << halo_ << std::endl;
}

// Solve method with update to LHS matrix
void LocalSpMatDnVec::solveWithUpdate(
  const double max_error,
  const double max_rel_error,
  const int max_restarts)
{
  solver_->solveWithUpdate(max_error, max_rel_error, max_restarts);
}
// Solve method without update to LHS matrix
void LocalSpMatDnVec::solveNoUpdate(
  const double max_error,
  const double max_rel_error,
  const int max_restarts)
{
  solver_->solveNoUpdate(max_error, max_rel_error, max_restarts);
}
