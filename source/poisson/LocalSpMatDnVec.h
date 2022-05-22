#pragma once

#include <map>
#include <set>
#include <vector>
#include <memory>
#include <mpi.h>

class SpRowInfo
{
  public:
    const int rank_;
    const long long idx_; // global row index
    std::map<long long, double> loc_colval_; // col_idx->val map
    std::map<long long, double> bd_colval_;
    // neirank_cols_[i] holds {rank of non-local col, idx of non-local col}
    std::vector<std::pair<int,long long>> neirank_cols_;

    SpRowInfo(const int rank, const long long row_idx, const int neirank_max) : rank_(rank), idx_(row_idx) 
    { 
      neirank_cols_.reserve(neirank_max); 
    }
    ~SpRowInfo() = default;

    void mapColVal(const long long col_idx, const double val) 
    { 
      loc_colval_[col_idx] += val; 
    }
    void mapColVal(const int rank, const long long col_idx, const double val) 
    {
      if (rank == rank_)
        mapColVal(col_idx, val);
      else
      {
        bd_colval_[col_idx] += val;
        neirank_cols_.push_back({rank, col_idx});
      }
    }
};

// Forward declaration of BiCGSTABSolver class
class BiCGSTABSolver;

class LocalSpMatDnVec 
{
  public:
    LocalSpMatDnVec(MPI_Comm m_comm, const int BLEN, const std::vector<double>& P_inv); 
    ~LocalSpMatDnVec();

    // Reserve space for linear system
    void reserve(const int N);
    // Push back value to COO matrix, up to user to insure ordering of row and column elements
    void cooPushBackVal(const double val, const long long row, const long long col);
    // Push back row to COO matrix, up to user to ensure ordering of rows
    void cooPushBackRow(const SpRowInfo &row);
    // Make the distributed linear system for solver
    void make(const std::vector<long long>& Nrows_xcumsum);

    // Solve method with update to LHS matrix
    void solveWithUpdate(
      const double max_error,
      const double max_rel_error,
      const int max_restarts); 

    // Solve method without update to LHS matrix
    void solveNoUpdate(
      const double max_error,
      const double max_rel_error,
      const int max_restarts); 

    // Modifiable references for x and b for setting and getting initial conditions/solution
    std::vector<double>& get_x() { return x_; }
    std::vector<double>& get_b() { return b_; }
    std::vector<double>& get_pScale() { return pScale_; }

    // Expose private variables to friendly solver
    friend class BiCGSTABSolver;

  private:
    int rank_;
    MPI_Comm m_comm_;
    int comm_size_; 
    const int BLEN_; // number of cells in a block

    int m_;
    int halo_;
    int loc_nnz_;
    int bd_nnz_;

    // Local rows of linear system + dense vecs
    std::vector<double> loc_cooValA_;
    std::vector<long long> loc_cooRowA_long_;
    std::vector<long long> loc_cooColA_long_;
    std::vector<double> x_;
    std::vector<double> b_;
    std::vector<double> pScale_; // post-preconditioning scaling factor for each block

    // Non-local rows with columns belonging to halo using rank-local indexing
    std::vector<double> bd_cooValA_;
    std::vector<long long> bd_cooRowA_long_;
    std::vector<long long> bd_cooColA_long_;

    // Identical to above, but with integers
    std::vector<int> loc_cooRowA_int_;
    std::vector<int> loc_cooColA_int_;
    std::vector<int> bd_cooRowA_int_;
    std::vector<int> bd_cooColA_int_;

    // bd_recv_set_[r] contains columns that need to be received from rank 'r'
    std::vector<std::set<long long>> bd_recv_set_;
    std::vector<std::vector<long long>> bd_recv_vec_;

    // Vectors that contain rules for sending and receiving
    std::vector<int> recv_ranks_;
    std::vector<int> recv_offset_;
    std::vector<int> recv_sz_;

    std::vector<int> send_ranks_;
    std::vector<int> send_offset_;
    std::vector<int> send_sz_;
    std::vector<int> send_pack_idx_;

    std::unique_ptr<BiCGSTABSolver> solver_;
};
