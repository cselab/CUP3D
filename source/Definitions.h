//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

// https://github.com/open-mpi/ompi/issues/5157#issuecomment-388495496
#define OMPI_SKIP_MPICXX 1 //silence annoying openmpi warnings

#include "Base.h"

#include "Utils/AlignedAllocator.h"

#include <Cubism/Grid.h>
#include <Cubism/GridMPI.h>
#include <Cubism/BlockInfo.h>
#include <Cubism/BlockLab.h>
#include <Cubism/BlockLabMPI.h>
#include <Cubism/Definitions.h>
#include <Cubism/AMR_MeshAdaptationMPI.h>

using namespace cubism;

#ifndef CUP_BLOCK_SIZEX
#define CUP_BLOCK_SIZEX 8
#define CUP_BLOCK_SIZEY 8
#define CUP_BLOCK_SIZEZ 8
#endif
#include <array>
#include <cassert>
#include <iosfwd>
#include <string>

CubismUP_3D_NAMESPACE_BEGIN

enum BCflag {freespace, periodic, wall};
inline BCflag string2BCflag(const std::string &strFlag)
{
  if      (strFlag == "periodic" ) return periodic;
  else if (strFlag == "wall"     ) return wall;
  else if (strFlag == "freespace") return freespace;
  else
  {
    fprintf(stderr,"BC not recognized %s\n",strFlag.c_str());
    fflush(0);abort();
    return periodic; // dummy
  }
}

//CAREFUL THESE ARE GLOBAL VARIABLES!
extern BCflag cubismBCX;
extern BCflag cubismBCY;
extern BCflag cubismBCZ;

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBC: public cubism::BlockLab<BlockType,allocator>
{
  typedef typename BlockType::ElementType ElementTypeBlock;
  static constexpr int sizeX = BlockType::sizeX;
  static constexpr int sizeY = BlockType::sizeY;
  static constexpr int sizeZ = BlockType::sizeZ;

  // Used for Boundary Conditions:
  // Apply bc on face of direction dir and side side (0 or 1):
  template<int dir, int side> void applyBCfaceOpen(const bool coarse = false)
  {
    //First we apply zero Neumann BCs to all variables. Then, if we are dealing with a 
    //VectorElement (which has DIM=3), we set the normal velocity to zero.
    //Second order Neumann   BCs mean than u_{i} = u_{i+1}
    //Second order Dirichlet BCs mean than u_{i} = -u_{i+1}
    if (!coarse)
    {
      auto * const cb = this->m_cacheBlock;
  
      int s[3] = {0,0,0}, e[3] = {0,0,0};
      const int* const stenBeg = this->m_stencilStart;
      const int* const stenEnd = this->m_stencilEnd;
      s[0] =  dir==0 ? (side==0 ? stenBeg[0] : sizeX ) : 0;
      s[1] =  dir==1 ? (side==0 ? stenBeg[1] : sizeY ) : 0;
      s[2] =  dir==2 ? (side==0 ? stenBeg[2] : sizeZ ) : 0;
      e[0] =  dir==0 ? (side==0 ? 0 : sizeX + stenEnd[0]-1 ) : sizeX;
      e[1] =  dir==1 ? (side==0 ? 0 : sizeY + stenEnd[1]-1 ) : sizeY;
      e[2] =  dir==2 ? (side==0 ? 0 : sizeZ + stenEnd[2]-1 ) : sizeZ;
      for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
      for(int ix=s[0]; ix<e[0]; ix++)
      {
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) = cb->Access
          (
            ( dir==0 ? (side==0 ? 0 : sizeX-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ-1 ) : iz ) - stenBeg[2]
          );      
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(dir) = (-1.)*cb->Access
          (
            ( dir==0 ? (side==0 ? 0 : sizeX-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ-1 ) : iz ) - stenBeg[2]
          ).member(dir);      
      }

      //tensorial edges and corners also filled
      const int aux = coarse ? 2:1;
      const int bsize[3] = {sizeX/aux, sizeY/aux, sizeZ/aux};
      int s_[3], e_[3];
      s_[dir] = stenBeg[dir]*(1-side) + bsize[dir]*side;
      e_[dir] = (bsize[dir]-1+stenEnd[dir])*side;
      const int d1 = (dir + 1) % 3;
      const int d2 = (dir + 2) % 3;
      for(int b=0; b<2; ++b)
      for(int a=0; a<2; ++a)
      {
        s_[d1] = stenBeg[d1] + a*b*(bsize[d1] - stenBeg[d1]);
        s_[d2] = stenBeg[d2] + (a-a*b)*(bsize[d2] - stenBeg[d2]);
        e_[d1] = (1-b+a*b)*(bsize[d1] - 1 + stenEnd[d1]);
        e_[d2] = (a+b-a*b)*(bsize[d2] - 1 + stenEnd[d2]);
        for(int iz=s_[2]; iz<e_[2]; iz++)
        for(int iy=s_[1]; iy<e_[1]; iy++)
        for(int ix=s_[0]; ix<e_[0]; ix++)
        {
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) = dir==0?
          cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) : (dir==1?
          cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]) :
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]));
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(dir) = dir==0?
          (-1.)*cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(dir) : (dir==1?
          (-1.)*cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(dir) :
          (-1.)*cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(dir));
        }
      }
    }
    else
    {
      auto * const cb = this->m_CoarsenedBlock;
  
      int s[3] = {0,0,0}, e[3] = {0,0,0};
      const int eI[3] = {(this->m_stencilEnd[0])/2 + 1 + this->m_InterpStencilEnd[0] -1,
                         (this->m_stencilEnd[1])/2 + 1 + this->m_InterpStencilEnd[1] -1,
                         (this->m_stencilEnd[2])/2 + 1 + this->m_InterpStencilEnd[2] -1};
      const int sI[3] = {(this->m_stencilStart[0]-1)/2+  this->m_InterpStencilStart[0],
                         (this->m_stencilStart[1]-1)/2+  this->m_InterpStencilStart[1],
                         (this->m_stencilStart[2]-1)/2+  this->m_InterpStencilStart[2]};

      const int* const stenBeg = sI;
      const int* const stenEnd = eI;
      s[0] =  dir==0 ? (side==0 ? stenBeg[0] : sizeX/2 ) : 0;
      s[1] =  dir==1 ? (side==0 ? stenBeg[1] : sizeY/2 ) : 0;
      s[2] =  dir==2 ? (side==0 ? stenBeg[2] : sizeZ/2 ) : 0;
      e[0] =  dir==0 ? (side==0 ? 0 : sizeX/2 + stenEnd[0]-1 ) : sizeX/2;
      e[1] =  dir==1 ? (side==0 ? 0 : sizeY/2 + stenEnd[1]-1 ) : sizeY/2;
      e[2] =  dir==2 ? (side==0 ? 0 : sizeZ/2 + stenEnd[2]-1 ) : sizeZ/2;
      for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
      for(int ix=s[0]; ix<e[0]; ix++)
      {
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) = cb->Access
          ( ( dir==0 ? (side==0 ? 0 : sizeX/2-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY/2-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ/2-1 ) : iz ) - stenBeg[2]);
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(dir) = (-1.)*cb->Access
          ( ( dir==0 ? (side==0 ? 0 : sizeX/2-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY/2-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ/2-1 ) : iz ) - stenBeg[2]).member(dir);
      }

      //tensorial edges and corners also filled (this is necessary for the coarse block!)
      const int aux = coarse ? 2:1;
      const int bsize[3] = {sizeX/aux, sizeY/aux, sizeZ/aux};
      int s_[3], e_[3];
      s_[dir] = stenBeg[dir]*(1-side) + bsize[dir]*side;
      e_[dir] = (bsize[dir]-1+stenEnd[dir])*side;
      const int d1 = (dir + 1) % 3;
      const int d2 = (dir + 2) % 3;
      for(int b=0; b<2; ++b)
      for(int a=0; a<2; ++a)
      {
        s_[d1] = stenBeg[d1] + a*b*(bsize[d1] - stenBeg[d1]);
        s_[d2] = stenBeg[d2] + (a-a*b)*(bsize[d2] - stenBeg[d2]);
        e_[d1] = (1-b+a*b)*(bsize[d1] - 1 + stenEnd[d1]);
        e_[d2] = (a+b-a*b)*(bsize[d2] - 1 + stenEnd[d2]);
        for(int iz=s_[2]; iz<e_[2]; iz++)
        for(int iy=s_[1]; iy<e_[1]; iy++)
        for(int ix=s_[0]; ix<e_[0]; ix++)
        {
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) = dir==0?
          cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) : (dir==1?
          cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]) :
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]));
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(dir) = dir==0?
          (-1.)*cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(dir) : (dir==1?
          (-1.)*cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(dir) :
          (-1.)*cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(dir));
        }
      }
    }
  }
  template<int dir, int side> void applyBCfaceWall(const bool coarse=false)
  {
    if (!coarse)
    {
      auto * const cb = this->m_cacheBlock;
      int s[3] = {0,0,0}, e[3] = {0,0,0};
      const int* const stenBeg = this->m_stencilStart;
      const int* const stenEnd = this->m_stencilEnd;
      s[0] =  dir==0 ? (side==0 ? stenBeg[0] : sizeX ) : 0;
      s[1] =  dir==1 ? (side==0 ? stenBeg[1] : sizeY ) : 0;
      s[2] =  dir==2 ? (side==0 ? stenBeg[2] : sizeZ ) : 0;
      e[0] =  dir==0 ? (side==0 ? 0 : sizeX + stenEnd[0]-1 ) : sizeX;
      e[1] =  dir==1 ? (side==0 ? 0 : sizeY + stenEnd[1]-1 ) : sizeY;
      e[2] =  dir==2 ? (side==0 ? 0 : sizeZ + stenEnd[2]-1 ) : sizeZ;
      for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
      for(int ix=s[0]; ix<e[0]; ix++)
      for(int k=0; k<ElementTypeBlock::DIM; k++)
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = (-1.0)*cb->Access
          (
            ( dir==0 ? (side==0 ? 0 : sizeX-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ-1 ) : iz ) - stenBeg[2]
          ).member(k);

      //tensorial edges and corners also filled
      const int aux = coarse ? 2:1;
      const int bsize[3] = {sizeX/aux, sizeY/aux, sizeZ/aux};
      int s_[3], e_[3];
      s_[dir] = stenBeg[dir]*(1-side) + bsize[dir]*side;
      e_[dir] = (bsize[dir]-1+stenEnd[dir])*side;
      const int d1 = (dir + 1) % 3;
      const int d2 = (dir + 2) % 3;
      for(int b=0; b<2; ++b)
      for(int a=0; a<2; ++a)
      {
        s_[d1] = stenBeg[d1] + a*b*(bsize[d1] - stenBeg[d1]);
        s_[d2] = stenBeg[d2] + (a-a*b)*(bsize[d2] - stenBeg[d2]);
        e_[d1] = (1-b+a*b)*(bsize[d1] - 1 + stenEnd[d1]);
        e_[d2] = (a+b-a*b)*(bsize[d2] - 1 + stenEnd[d2]);
        for(int iz=s_[2]; iz<e_[2]; iz++)
        for(int iy=s_[1]; iy<e_[1]; iy++)
        for(int ix=s_[0]; ix<e_[0]; ix++)
        for(int k=0; k<ElementTypeBlock::DIM; k++)
        {
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = (-1.0) * (dir==0?
          cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) : (dir==1?
          cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(k):
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(k)));
        }
      }
    }
    else
    {
      auto * const cb = this->m_CoarsenedBlock;
  
      int s[3] = {0,0,0}, e[3] = {0,0,0};
      const int eI[3] = {(this->m_stencilEnd[0])/2 + 1 + this->m_InterpStencilEnd[0] -1,
                         (this->m_stencilEnd[1])/2 + 1 + this->m_InterpStencilEnd[1] -1,
                         (this->m_stencilEnd[2])/2 + 1 + this->m_InterpStencilEnd[2] -1};
      const int sI[3] = {(this->m_stencilStart[0]-1)/2+  this->m_InterpStencilStart[0],
                         (this->m_stencilStart[1]-1)/2+  this->m_InterpStencilStart[1],
                         (this->m_stencilStart[2]-1)/2+  this->m_InterpStencilStart[2]};

      const int* const stenBeg = sI;
      const int* const stenEnd = eI;
      s[0] =  dir==0 ? (side==0 ? stenBeg[0] : sizeX/2 ) : 0;
      s[1] =  dir==1 ? (side==0 ? stenBeg[1] : sizeY/2 ) : 0;
      s[2] =  dir==2 ? (side==0 ? stenBeg[2] : sizeZ/2 ) : 0;
      e[0] =  dir==0 ? (side==0 ? 0 : sizeX/2 + stenEnd[0]-1 ) : sizeX/2;
      e[1] =  dir==1 ? (side==0 ? 0 : sizeY/2 + stenEnd[1]-1 ) : sizeY/2;
      e[2] =  dir==2 ? (side==0 ? 0 : sizeZ/2 + stenEnd[2]-1 ) : sizeZ/2;
      for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
      for(int ix=s[0]; ix<e[0]; ix++)
      for(int k=0; k<ElementTypeBlock::DIM; k++)
      {
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = (-1.0)*cb->Access
          ( ( dir==0 ? (side==0 ? 0 : sizeX/2-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY/2-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ/2-1 ) : iz ) - stenBeg[2]).member(k);
      }

      //tensorial edges and corners also filled (this is necessary for the coarse block!)
      const int aux = coarse ? 2:1;
      const int bsize[3] = {sizeX/aux, sizeY/aux, sizeZ/aux};
      int s_[3], e_[3];
      s_[dir] = stenBeg[dir]*(1-side) + bsize[dir]*side;
      e_[dir] = (bsize[dir]-1+stenEnd[dir])*side;
      const int d1 = (dir + 1) % 3;
      const int d2 = (dir + 2) % 3;
      for(int b=0; b<2; ++b)
      for(int a=0; a<2; ++a)
      {
        s_[d1] = stenBeg[d1] + a*b*(bsize[d1] - stenBeg[d1]);
        s_[d2] = stenBeg[d2] + (a-a*b)*(bsize[d2] - stenBeg[d2]);
        e_[d1] = (1-b+a*b)*(bsize[d1] - 1 + stenEnd[d1]);
        e_[d2] = (a+b-a*b)*(bsize[d2] - 1 + stenEnd[d2]);
        for(int iz=s_[2]; iz<e_[2]; iz++)
        for(int iy=s_[1]; iy<e_[1]; iy++)
        for(int ix=s_[0]; ix<e_[0]; ix++)
        for(int k=0; k<ElementTypeBlock::DIM; k++)
        {
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = (-1.0) * (dir==0?
          cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) : (dir==1?
          cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(k) :
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(k)));
        }
      }
    }
  }

 public:

  typedef typename BlockType::ElementType ElementType;

  virtual bool is_xperiodic() override { return cubismBCX == periodic; }
  virtual bool is_yperiodic() override { return cubismBCY == periodic; }
  virtual bool is_zperiodic() override { return cubismBCZ == periodic; }

  BlockLabBC() = default;
  BlockLabBC(const BlockLabBC&) = delete;
  BlockLabBC& operator=(const BlockLabBC&) = delete;

  // Called by Cubism:
  // wall means u=v=w=0 and dp/dn=0
  // freespace means normal velocity=0, dp/dn=0 and d(tangential velocities)/dn=0
  void _apply_bc(const cubism::BlockInfo& info, const Real t=0, const bool coarse = false)
  {
    const BCflag BCX = cubismBCX;
    const BCflag BCY = cubismBCY;
    const BCflag BCZ = cubismBCZ;
    if (BCX == wall)
    {
      if(info.index[0]==0 )          this->template applyBCfaceWall<0,0>(coarse);
      if(info.index[0]==this->NX-1 ) this->template applyBCfaceWall<0,1>(coarse);
    }
    else if (BCX != periodic)
    {
      //then either we have freespace or wall. If we have wall, we still apply 
      //freespace to the PoissonElement, as pressure has zero Neumann BCS
      if(info.index[0]==0 )          this->template applyBCfaceOpen<0,0>(coarse);
      if(info.index[0]==this->NX-1 ) this->template applyBCfaceOpen<0,1>(coarse);
    }
    if (BCY == wall)
    {
      if(info.index[1]==0 )          this->template applyBCfaceWall<1,0>(coarse);
      if(info.index[1]==this->NY-1 ) this->template applyBCfaceWall<1,1>(coarse);
    }
    else if (BCY != periodic)
    {
      //then either we have freespace or wall. If we have wall, we still apply 
      //freespace to the PoissonElement, as pressure has zero Neumann BCS
      if(info.index[1]==0 )          this->template applyBCfaceOpen<1,0>(coarse);
      if(info.index[1]==this->NY-1 ) this->template applyBCfaceOpen<1,1>(coarse);
    }
    if (BCZ == wall)
    {
      if(info.index[2]==0 )          this->template applyBCfaceWall<2,0>(coarse);
      if(info.index[2]==this->NZ-1 ) this->template applyBCfaceWall<2,1>(coarse);
    }
    else if (BCZ != periodic)
    {
      //then either we have freespace or wall. If we have wall, we still apply 
      //freespace to the PoissonElement, as pressure has zero Neumann BCS
      if(info.index[2]==0 )          this->template applyBCfaceOpen<2,0>(coarse);
      if(info.index[2]==this->NZ-1 ) this->template applyBCfaceOpen<2,1>(coarse);
    }
  }
};


template<typename BlockType, template<typename X> class allocator = std::allocator>
class BlockLabNeumann3D: public cubism::BlockLabNeumann<BlockType, 3, allocator>
{
 public:
   using cubismLab = cubism::BlockLabNeumann<BlockType, 3, allocator>;
   virtual bool is_xperiodic() override{ return cubismBCX == periodic; }
   virtual bool is_yperiodic() override{ return cubismBCY == periodic; }
   virtual bool is_zperiodic() override{ return cubismBCZ == periodic; }

   // Called by Cubism:
   void _apply_bc(const cubism::BlockInfo& info, const Real t = 0, const bool coarse = false) override
   {
       if (is_xperiodic() == false)
       {
        if(info.index[0]==0 )          cubismLab::template Neumann3D<0,0>(coarse);
        if(info.index[0]==this->NX-1 ) cubismLab::template Neumann3D<0,1>(coarse);
       }
       if (is_yperiodic() == false)
       {
        if(info.index[1]==0 )          cubismLab::template Neumann3D<1,0>(coarse);
        if(info.index[1]==this->NY-1 ) cubismLab::template Neumann3D<1,1>(coarse);
       }
       if (is_zperiodic() == false)
       {
        if(info.index[2]==0 )          cubismLab::template Neumann3D<2,0>(coarse);
        if(info.index[2]==this->NZ-1 ) cubismLab::template Neumann3D<2,1>(coarse);
       }
   }
};

struct StreamerVectorX
{
  static constexpr int NCHANNELS = 1;
  template <typename TBlock, typename T>
  static void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
  {
      output[0] = b(ix,iy,iz).u[0];
  }
  static std::string prefix() { return std::string(""); }
  static const char * getAttributeName() { return "VectorX"; }
};
struct StreamerVectorY
{
  static constexpr int NCHANNELS = 1;
  template <typename TBlock, typename T>
  static void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
  {
      output[0] = b(ix,iy,iz).u[1];
  }
  static std::string prefix() { return std::string(""); }
  static const char * getAttributeName() { return "VectorY"; }
};
struct StreamerVectorZ
{
  static constexpr int NCHANNELS = 1;
  template <typename TBlock, typename T>
  static void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
  {
      output[0] = b(ix,iy,iz).u[2];
  }
  static std::string prefix() { return std::string(""); }
  static const char * getAttributeName() { return "VectorZ"; }
};

// The alignmnet of 32 is sufficient for AVX-256, but we put 64 to cover AVX-512.
static constexpr int kBlockAlignment = 64;
template <typename T>
using aligned_block_allocator = aligned_allocator<T, kBlockAlignment>;

using ScalarElement = cubism::ScalarElement<Real>;
using ScalarBlock   = cubism::GridBlock<CUP_BLOCK_SIZEX,3,ScalarElement>;
using ScalarGrid    = cubism::GridMPI<cubism::Grid<ScalarBlock, aligned_block_allocator>>;
using ScalarLab     = cubism::BlockLabMPI<BlockLabNeumann3D <ScalarBlock, aligned_block_allocator>,ScalarGrid>;

using VectorElement = cubism::VectorElement<3,Real>;
using VectorBlock   = cubism::GridBlock<CUP_BLOCK_SIZEX,3,VectorElement>;
using VectorGrid    = cubism::GridMPI<cubism::Grid<VectorBlock, aligned_block_allocator>>;
using VectorLab     = cubism::BlockLabMPI<BlockLabBC<VectorBlock, aligned_block_allocator>, VectorGrid>;

using ScalarAMR     = cubism::MeshAdaptationMPI<ScalarGrid,ScalarLab,ScalarGrid>;
using VectorAMR     = cubism::MeshAdaptationMPI<VectorGrid,VectorLab,ScalarGrid>;

CubismUP_3D_NAMESPACE_END
