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
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + iy*nx_ + ix); }

      // Return indices of neighbouring cells
      long long neiXp(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + iy*nx_ + ix+dist); }
      long long neiXm(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + iy*nx_ + ix-dist); }
      long long neiYp(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + (iy+dist)*nx_ + ix); }
      long long neiYm(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + (iy-dist)*nx_ + ix); }
      long long neiZp(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)((iz+dist)*ny_*nx_ + iy*nx_ + ix); }
      long long neiZm(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const
      { return blockOffset(info) + (long long)((iz-dist)*ny_*nx_ + iy*nx_ + ix); }

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
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + iy*nx_ + (nx_-1-offset)); }
      long long Xmin(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + iy*nx_ + offset); }
      long long Ymax(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + (ny_-1-offset)*nx_ + ix); }
      long long Ymin(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(iz*ny_*nx_ + offset*nx_ + ix); }
      long long Zmax(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)((nz_-1-offset)*ny_*nx_ + iy*nx_ + ix); }
      long long Zmin(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int offset = 0) const
      { return blockOffset(info) + (long long)(offset*ny_*nx_ + iy*nx_ + ix); }

    protected:
      long long blockOffset(const cubism::BlockInfo &info) const
      { return (info.blockID + ps.Nblocks_xcumsum_[ps.sim.lhs->Tree(info).rank()])*nxyz_; }
      static int ix_f(const int &ix) { return (ix % (nx_/2)) * 2; }
      static int iy_f(const int &iy) { return (iy % (ny_/2)) * 2; }
      static int iz_f(const int &iz) { return (iz % (nz_/2)) * 2; }

      const PoissonSolverExp &ps;
  };

  /*
    Class to help with indexing of neighbours for cell located at a face of a block.
    To provide a generic API for flux constructors at block faces, static polymorphism determines how indexing 
    will be carried out with methods overloaded in [XYZ]{max,min}Indexer.

      - neiUnif gives the index of a cell across the face assuming the neighbouring block is at the same refinement level

    Define [XYZ]{max,min}Indexer using an inheritence model which needs 10 class definitions instead of 6, 
    but avoids two-fold reimplementation of some methods
  */

  // Abstract base class
  class FaceCellIndexer : public CellIndexer
  {
    public:
      FaceCellIndexer(const PoissonSolverExp& pSolver) : CellIndexer(pSolver) {}

      // When I am uniform with the neighbouring block
      virtual long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const = 0;

      // When I am finer than neighbouring block
      virtual long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const = 0;
      virtual double sign_ds1(const int &ix, const int &iy, const int &iz) const = 0;
      virtual double sign_ds2(const int &ix, const int &iy, const int &iz) const = 0;

      virtual int ix_c(const cubism::BlockInfo &info, const int &ix) const
      { return (info.index[0] % 2 == 0) ? (ix/2) : (ix/2 + nx_/2); }
      virtual int iy_c(const cubism::BlockInfo &info, const int &iy) const
      { return (info.index[1] % 2 == 0) ? (iy/2) : (iy/2 + ny_/2); }
      virtual int iz_c(const cubism::BlockInfo &info, const int &iz) const
      { return (info.index[2] % 2 == 0) ? (iz/2) : (iz/2 + nz_/2); }

      // When I am coarser than neiboughring block
      // Order of neiFine{1,2,3,4} must be consistent with [sign_ds1, sign_ds2] of corresponding dimension and have order [-1,-1], [1,-1], [-1,1], [1,1]
      virtual long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const = 0;
      virtual long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const = 0;
      virtual long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const = 0;
      virtual long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const = 0;

      // Indexing aids for derivatives in taylor approximation, info, ix, iy, iz of a coarse cell
      virtual bool isBD_ds1(const int &ix, const int &iy, const int &iz) const = 0;
      virtual bool isFD_ds1(const int &ix, const int &iy, const int &iz) const = 0;
      virtual bool isBD_ds2(const int &ix, const int &iy, const int &iz) const = 0;
      virtual bool isFD_ds2(const int &ix, const int &iy, const int &iz) const = 0;
      virtual long long backward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const = 0;
      virtual long long forward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const = 0;
      virtual long long backward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const = 0;
      virtual long long forward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const = 0;

      // When I am coarses and need to determine which Zchild I'm next to
      virtual long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const = 0;

  };

  // ------------------------------------------- Faces perpendicular to x-axis -----------------------------------------------
  // ds1 = y, ds2 = z
  class XbaseIndexer : public FaceCellIndexer
  {
    public:
      XbaseIndexer(const PoissonSolverExp& pSolver) : FaceCellIndexer(pSolver) {}

      double sign_ds1(const int &ix, const int &iy, const int &iz) const override
      { return iy % 2 == 0 ? -1. : 1.; }
      double sign_ds2(const int &ix, const int &iy, const int &iz) const override
      { return iz % 2 == 0 ? -1. : 1.; }

      bool isBD_ds1(const int &ix, const int &iy, const int &iz) const override
      { return iy == ny_-1; }
      bool isFD_ds1(const int &ix, const int &iy, const int &iz) const override
      { return iy == 0; }
      bool isBD_ds2(const int &ix, const int &iy, const int &iz) const override
      { return iz == nz_-1; }
      bool isFD_ds2(const int &ix, const int &iy, const int &iz) const override
      { return iz == 0; }
      long long backward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiYm(info, ix, iy, iz, dist); }
      long long forward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiYp(info, ix, iy, iz, dist); }
      long long backward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiZm(info, ix, iy, iz, dist); }
      long long forward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiZp(info, ix, iy, iz, dist); }

  };

  class XminIndexer : public XbaseIndexer
  {
    public:
      XminIndexer(const PoissonSolverExp& pSolver) : XbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmax(nei_info, ix, iy, iz, offset); }

      long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiXp(info, ix, iy, iz, dist); }

      int ix_c(const cubism::BlockInfo &info, const int &ix) const override
      { return nx_ - 1; }

      long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmax(nei_info, ix_f(ix), iy_f(iy), iz_f(iz), offset); }
      long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmax(nei_info, ix_f(ix), iy_f(iy)+1, iz_f(iz), offset); }
      long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmax(nei_info, ix_f(ix), iy_f(iy), iz_f(iz)+1, offset); }
      long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmax(nei_info, ix_f(ix), iy_f(iy)+1, iz_f(iz)+1, offset); }

      long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const override
      { return nei_info.Zchild[1][int(iy >= ny_/2)][int(iz >= nz_/2)];}
  };

  class XmaxIndexer : public XbaseIndexer
  {
    public:
      XmaxIndexer(const PoissonSolverExp& pSolver) : XbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmin(nei_info, ix, iy, iz, offset); }

      long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiXm(info, ix, iy, iz, dist); }

      int ix_c(const cubism::BlockInfo &info, const int &ix) const override
      { return 0; }

      long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmin(nei_info, ix_f(ix), iy_f(iy), iz_f(iz), offset); }
      long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmin(nei_info, ix_f(ix), iy_f(iy)+1, iz_f(iz), offset); }
      long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmin(nei_info, ix_f(ix), iy_f(iy), iz_f(iz)+1, offset); }
      long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Xmin(nei_info, ix_f(ix), iy_f(iy)+1, iz_f(iz)+1, offset); }

      long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const override
      { return nei_info.Zchild[0][int(iy >= ny_/2)][int(iz >= nz_/2)];}
  };

  // ------------------------------------------- Faces perpendicular to y-axis -----------------------------------------------
  // ds1 = z, ds2 = y
  class YbaseIndexer : public FaceCellIndexer
  {
    public:
      YbaseIndexer(const PoissonSolverExp& pSolver) : FaceCellIndexer(pSolver) {}

      double sign_ds1(const int &ix, const int &iy, const int &iz) const override
      { return iz % 2 == 0 ? -1. : 1.; }
      double sign_ds2(const int &ix, const int &iy, const int &iz) const override
      { return ix % 2 == 0 ? -1. : 1.; }

      bool isBD_ds1(const int &ix, const int &iy, const int &iz) const override
      { return iz == nz_-1; }
      bool isFD_ds1(const int &ix, const int &iy, const int &iz) const override
      { return iz == 0; }
      bool isBD_ds2(const int &ix, const int &iy, const int &iz) const override
      { return ix == nx_-1; }
      bool isFD_ds2(const int &ix, const int &iy, const int &iz) const override
      { return ix == 0; }
      long long backward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiZm(info, ix, iy, iz, dist); }
      long long forward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiZp(info, ix, iy, iz, dist); }
      long long backward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiXm(info, ix, iy, iz, dist); }
      long long forward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiXp(info, ix, iy, iz, dist); }

  };

  class YminIndexer : public YbaseIndexer
  {
    public:
      YminIndexer(const PoissonSolverExp& pSolver) : YbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymax(nei_info, ix, iy, iz, offset); }

      long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiYp(info, ix, iy, iz, dist); }

      int iy_c(const cubism::BlockInfo &info, const int &iy) const override
      { return ny_ - 1; }

      long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymax(nei_info, ix_f(ix), iy_f(iy), iz_f(iz), offset); }
      long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymax(nei_info, ix_f(ix), iy_f(iy), iz_f(iz)+1, offset); } // z is ds1, modulate it first
      long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymax(nei_info, ix_f(ix)+1, iy_f(iy), iz_f(iz), offset); }
      long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymax(nei_info, ix_f(ix)+1, iy_f(iy), iz_f(iz)+1, offset); }

      long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const override
      { return nei_info.Zchild[int(ix >= nx_/2)][1][int(iz >= nz_/2)];}
  };

  class YmaxIndexer : public YbaseIndexer
  {
    public:
      YmaxIndexer(const PoissonSolverExp& pSolver) : YbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymin(nei_info, ix, iy, iz, offset); }

      long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiYm(info, ix, iy, iz, dist); }

      int iy_c(const cubism::BlockInfo &info, const int &iy) const override
      { return 0; }

      long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymin(nei_info, ix_f(ix), iy_f(iy), iz_f(iz), offset); }
      long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymin(nei_info, ix_f(ix), iy_f(iy), iz_f(iz)+1, offset); } // z is ds1, modulate it first
      long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymin(nei_info, ix_f(ix)+1, iy_f(iy), iz_f(iz), offset); }
      long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Ymin(nei_info, ix_f(ix)+1, iy_f(iy), iz_f(iz)+1, offset); }

      long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const override
      { return nei_info.Zchild[int(ix >= nx_/2)][0][int(iz >= nz_/2)];}
  };

  // ------------------------------------------- Faces perpendicular to z-axis -----------------------------------------------
  // ds1 = x, ds2 = y
  class ZbaseIndexer : public FaceCellIndexer
  {
    public:
      ZbaseIndexer(const PoissonSolverExp& pSolver) : FaceCellIndexer(pSolver) {}

      double sign_ds1(const int &ix, const int &iy, const int &iz) const override
      { return ix % 2 == 0 ? -1. : 1.; }
      double sign_ds2(const int &ix, const int &iy, const int &iz) const override
      { return iy % 2 == 0 ? -1. : 1.; }

      bool isBD_ds1(const int &ix, const int &iy, const int &iz) const override
      { return ix == nx_-1; }
      bool isFD_ds1(const int &ix, const int &iy, const int &iz) const override
      { return ix == 0; }
      bool isBD_ds2(const int &ix, const int &iy, const int &iz) const override
      { return iy == ny_-1; }
      bool isFD_ds2(const int &ix, const int &iy, const int &iz) const override
      { return iy == 0; }
      long long backward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiXm(info, ix, iy, iz, dist); }
      long long forward_ds1(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiXp(info, ix, iy, iz, dist); }
      long long backward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiYm(info, ix, iy, iz, dist); }
      long long forward_ds2(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiYp(info, ix, iy, iz, dist); }

  };

  class ZminIndexer : public ZbaseIndexer
  {
    public:
      ZminIndexer(const PoissonSolverExp& pSolver) : ZbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmax(nei_info, ix, iy, iz, offset); }

      long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiZp(info, ix, iy, iz, dist); }

      int iz_c(const cubism::BlockInfo &info, const int &iz) const override
      { return nz_ - 1; }

      long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmax(nei_info, ix_f(ix), iy_f(iy), iz_f(iz), offset); }
      long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmax(nei_info, ix_f(ix)+1, iy_f(iy), iz_f(iz), offset); }
      long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmax(nei_info, ix_f(ix), iy_f(iy)+1, iz_f(iz), offset); }
      long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmax(nei_info, ix_f(ix)+1, iy_f(iy)+1, iz_f(iz), offset); }

      long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const override
      { return nei_info.Zchild[int(ix >= nx_/2)][int(iy >= ny_/2)][1];}
  };

  class ZmaxIndexer : public ZbaseIndexer
  {
    public:
      ZmaxIndexer(const PoissonSolverExp& pSolver) : ZbaseIndexer(pSolver) {}

      long long neiUnif(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmin(nei_info, ix, iy, iz, offset); }

      long long neiInward(const cubism::BlockInfo &info, const int &ix, const int &iy, const int &iz, const int dist = 1) const override
      { return neiZm(info, ix, iy, iz, dist); }

      int iz_c(const cubism::BlockInfo &info, const int &iz) const override
      { return 0; }

      long long neiFine1(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmin(nei_info, ix_f(ix), iy_f(iy), iz_f(iz), offset); }
      long long neiFine2(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmin(nei_info, ix_f(ix)+1, iy_f(iy), iz_f(iz), offset); }
      long long neiFine3(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmin(nei_info, ix_f(ix), iy_f(iy)+1, iz_f(iz), offset); }
      long long neiFine4(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz, const int offset = 0) const override
      { return Zmin(nei_info, ix_f(ix)+1, iy_f(iy)+1, iz_f(iz), offset); }

      long long Zchild(const cubism::BlockInfo &nei_info, const int &ix, const int &iy, const int &iz) const override
      { return nei_info.Zchild[int(ix >= nx_/2)][int(iy >= ny_/2)][0];}
  };

  // Cell indexers
  CellIndexer GenericCell;
  XminIndexer XminCell;
  XmaxIndexer XmaxCell;
  YminIndexer YminCell;
  YmaxIndexer YmaxCell;
  ZminIndexer ZminCell;
  ZmaxIndexer ZmaxCell;
  // Array of the indexers above pointed to by a base pointer for polymorphism in makeFlux
  std::array<const FaceCellIndexer*, 6> faceIndexers;

  std::array<std::pair<long long, double>, 3> D1_ds1(const cubism::BlockInfo &info, const FaceCellIndexer &indexer, const int ix, const int iy, const int iz) const
  {
    // Scale D1 by h^l/4
    if (indexer.isBD_ds1(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.backward_ds1(info, ix, iy, iz, 2),  1./8.), 
              std::make_pair<long long, double>(indexer.backward_ds1(info, ix, iy, iz, 1), -1./2.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),             3./8.)};
    else if (indexer.isFD_ds1(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.forward_ds1(info, ix, iy, iz, 2), -1./8.), 
              std::make_pair<long long, double>(indexer.forward_ds1(info, ix, iy, iz, 1),  1./2.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),           -3./8.)};

    return {std::make_pair<long long, double>(indexer.backward_ds1(info, ix, iy, iz, 1), -1./8.), 
            std::make_pair<long long, double>(indexer.forward_ds1(info, ix, iy, iz, 1),   1./8.), 
            std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),             0.)};
  }

  std::array<std::pair<long long, double>, 3> D1_ds2(const cubism::BlockInfo &info, const FaceCellIndexer &indexer, const int ix, const int iy, const int iz) const
  {
    // Scale D1 by h^l/4
    if (indexer.isBD_ds2(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.backward_ds2(info, ix, iy, iz, 2),  1./8.), 
              std::make_pair<long long, double>(indexer.backward_ds2(info, ix, iy, iz, 1), -1./2.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),             3./8.)};
    else if (indexer.isFD_ds2(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.forward_ds2(info, ix, iy, iz, 2), -1./8.), 
              std::make_pair<long long, double>(indexer.forward_ds2(info, ix, iy, iz, 1),  1./2.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),           -3./8.)};

    return {std::make_pair<long long, double>(indexer.backward_ds2(info, ix, iy, iz, 1), -1./8.), 
            std::make_pair<long long, double>(indexer.forward_ds2(info, ix, iy, iz, 1),   1./8.), 
            std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),             0.)};
  }

  std::array<std::pair<long long, double>, 3> D2_ds1(const cubism::BlockInfo &info, const FaceCellIndexer &indexer, const int ix, const int iy, const int iz) const
  {
    // Scale D2 by 0.5*(h^l/4)^2
    if (indexer.isBD_ds1(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.backward_ds1(info, ix, iy, iz, 2),  1./32.), 
              std::make_pair<long long, double>(indexer.backward_ds1(info, ix, iy, iz, 1), -1./16.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),             1./32.)};
    else if (indexer.isFD_ds1(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.forward_ds1(info, ix, iy, iz, 2),  1./32.), 
              std::make_pair<long long, double>(indexer.forward_ds1(info, ix, iy, iz, 1), -1./16.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),            1./32.)};

    return {std::make_pair<long long, double>(indexer.backward_ds1(info, ix, iy, iz, 1),  1./32.), 
            std::make_pair<long long, double>(indexer.forward_ds1(info, ix, iy, iz, 1),   1./32.), 
            std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),            -1./16.)};
  }

  std::array<std::pair<long long, double>, 3> D2_ds2(const cubism::BlockInfo &info, const FaceCellIndexer &indexer, const int ix, const int iy, const int iz) const
  {
    // Scale D2 by 0.5*(h^l/4)^2
    if (indexer.isBD_ds2(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.backward_ds2(info, ix, iy, iz, 2),  1./32.), 
              std::make_pair<long long, double>(indexer.backward_ds2(info, ix, iy, iz, 1), -1./16.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),             1./32.)};
    else if (indexer.isFD_ds2(ix, iy, iz)) 
      return {std::make_pair<long long, double>(indexer.forward_ds2(info, ix, iy, iz, 2),  1./32.), 
              std::make_pair<long long, double>(indexer.forward_ds2(info, ix, iy, iz, 1), -1./16.), 
              std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),            1./32.)};

    return {std::make_pair<long long, double>(indexer.backward_ds2(info, ix, iy, iz, 1),  1./32.), 
            std::make_pair<long long, double>(indexer.forward_ds2(info, ix, iy, iz, 1),   1./32.), 
            std::make_pair<long long, double>(indexer.This(info, ix, iy, iz),            -1./16.)};
  }

  void interpolate(
      const cubism::BlockInfo &info_c, const int ix_c, const int iy_c, const int iz_c,
      const cubism::BlockInfo &info_f, const long long fine_close_idx, const long long fine_far_idx,
      const double sign, const double sign_ds1, const double sign_ds2,
      const FaceCellIndexer &indexer, SpRowInfo& row) const;
};

}//namespace cubismup3d
