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
  static constexpr int nxyz_ = nx_ * ny_ * nz_;

  // Returns element of preconditioner negative K_{i,j}
  double getA_local(const int& i, const int& j);

  // Method construct flux accross block boundaries
  class FaceCellIndexer; // forward declaration
  void makeFlux(
      const cubism::BlockInfo &rhs_info,
      const int &ix, const int &iy, const int &iz,
      const bool &isBoundary,
      const cubism::BlockInfo &rhsNei,
      const FaceCellIndexer &indexer,
      SpRowInfo &row) const;

  // Method to compute A and b for the current mesh
  void getMat(); // update LHS and RHS after refinement
  void getVec(); // update initial guess and RHS vecs only

  // Distributed linear system with local indexing
  std::unique_ptr<LocalSpMatDnVec> LocalLS_;

  std::vector<long long> Nblocks_xcumsum_;
  std::vector<long long> Nrows_xcumsum_;

  class CellIndexer
  {
    public:
      CellIndexer(const PoissonSolverExp& pSolver) : ps(pSolver) {}
      ~CellIndexer() = default;

      long long This(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz) const
      { return blockOffset(info) + (long long)(iz*ny*nx + iy*nx + ix); }

      // Return indices of neighbouring cells
      long long neiXp(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny*nx + iy*nx + ix+dist); }
      long long neiXm(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny*nx + iy*nx + ix-dist); }
      long long neiYp(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny*nx + (iy+dist)*nx + ix); }
      long long neiYm(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny*nx + (iy-dist)*nx + ix); }
      long long neiZp(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)((iz+dist)*ny*nx + iy*nx + ix); }
      long long neiZm(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)((iz-dist)*ny*nx + iy*nx + ix); }

      // Check if neighbouring cells are in this block
      static bool validXp(const int &ix, const int &iy, const int &iz)
      { return ix < nx_ - 1; }
      static bool validXm(const int &ix, const int &iy, const int &iz)
      { return ix > 0; }
      static bool validYp(const int &ix, const int &iy, const int &iz)
      { return iy < ny_ - 1; }
      static bool validYm(const int &ix, const int &iy, const int &iz)
      { return iy > 0; }
      static bool validZp(const int &ix, const int &iy, const int &iz)
      { return iz < nz_ - 1; }
      static bool validZm(const int &ix, const int &iy, const int &iz)
      { return iz > 0; }

      long long Xmax(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny*nx + iy*nx + (nx-1-offset)); }
      long long Xmin(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny*nx + iy*nx + offset); }
      long long Ymax(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny*nx + (ny-1-offset)*nx + ix); }
      long long Ymin(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny*nx + offset*nx + ix); }
      long long Zmax(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)((nz-1-offset)*ny*nx + iy*nx + ix); }
      long long Zmin(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(offset*ny*nx + iy*nx + ix); }

    protected:
      long long blockOffset(const cubism::BlockInfo &info) const
      { return (info.blockID + ps.Nblocks_xcumsum_[ps.sim.lhs->Tree(info).rank()])*nlen; }

      const PoissonSolverExp &ps;
      static constexpr int nx = BlockType::sizeX;
      static constexpr int ny = BlockType::sizeY;
      static constexpr int nz = BlockType::sizeZ;
      static constexpr long long nlen = nx * ny * nz;
  };

  /*
    Class to help with indexing of neighbours for cell located at a face of a block.
    To provide a generic API for flux constructors at block faces, static polymorphism determines how indexing 
    will be carried out with methods overloaded in [XYZ]{max,min}Indexer.

      - neiUnif gives the index of a cell across the face assuming the neighbouring block is at the same refinement level

    Define [XYZ]{max,min}Indexer using an inheritence model which needs 10 class definitions instead of 6, 
    but method reuse makes it less error prone
  */

  // Abstract base class
  class FaceCellIndexer : public CellIndexer
  {
    public:
      FaceCellIndexer(const PoissonSolverExp& pSolver) : CellIndexer(pSolver) {}

      virtual long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const = 0;

  };

  // ------------------------------------------- Faces perpendicular to x-axis -----------------------------------------------
  class XbaseIndexer : public FaceCellIndexer
  {
    public:
      XbaseIndexer(const PoissonSolverExp& pSolver) : FaceCellIndexer(pSolver) {}
  };

  class XminIndexer : public XbaseIndexer
  {
    public:
      XminIndexer(const PoissonSolverExp& pSolver) : XbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmax(nei_info, ix, iy, iz, offset); }
  };

  class XmaxIndexer : public XbaseIndexer
  {
    public:
      XmaxIndexer(const PoissonSolverExp& pSolver) : XbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmin(nei_info, ix, iy, iz, offset); }
  };

  // ------------------------------------------- Faces perpendicular to y-axis -----------------------------------------------
  class YbaseIndexer : public FaceCellIndexer
  {
    public:
      YbaseIndexer(const PoissonSolverExp& pSolver) : FaceCellIndexer(pSolver) {}
  };

  class YminIndexer : public YbaseIndexer
  {
    public:
      YminIndexer(const PoissonSolverExp& pSolver) : YbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymax(nei_info, ix, iy, iz, offset); }
  };

  class YmaxIndexer : public YbaseIndexer
  {
    public:
      YmaxIndexer(const PoissonSolverExp& pSolver) : YbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymin(nei_info, ix, iy, iz, offset); }
  };

  // ------------------------------------------- Faces perpendicular to z-axis -----------------------------------------------
  class ZbaseIndexer : public FaceCellIndexer
  {
    public:
      ZbaseIndexer(const PoissonSolverExp& pSolver) : FaceCellIndexer(pSolver) {}
  };

  class ZminIndexer : public ZbaseIndexer
  {
    public:
      ZminIndexer(const PoissonSolverExp& pSolver) : ZbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmax(nei_info, ix, iy, iz, offset); }
  };

  class ZmaxIndexer : public ZbaseIndexer
  {
    public:
      ZmaxIndexer(const PoissonSolverExp& pSolver) : ZbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmin(nei_info, ix, iy, iz, offset); }
  };

  // Cell indexers
  CellIndexer GenericCell;
  XminIndexer XminCell;
  XmaxIndexer XmaxCell;
  YminIndexer YminCell;
  YmaxIndexer YmaxCell;
  ZminIndexer ZminCell;
  ZmaxIndexer ZmaxCell;
  // Array of the indexers above upcast to their parent class for polymorphism in makeFlux
  std::array<const FaceCellIndexer*, 6> faceIndexers;
};

}//namespace cubismup3d
