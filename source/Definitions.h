//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Base.h"

#include "utils/AlignedAllocator.h"

#include <Cubism/Grid.h>
#include <Cubism/GridMPI.h>
#include <Cubism/BlockInfo.h>
#include <Cubism/BlockLab.h>
#include <Cubism/BlockLabMPI.h>

#ifndef CUP_BLOCK_SIZEX
#define CUP_BLOCK_SIZEX 8
#define CUP_BLOCK_SIZEY 8
#define CUP_BLOCK_SIZEZ 8
#endif

#include <array>
#include <cassert>
#include <iosfwd>
#include <string>

#include "MeshAdaptation_CUP.h"

CubismUP_3D_NAMESPACE_BEGIN

struct PoissonElement
{
  using RealType = Real;
  Real s = 0;
  Real lhs = 0;
  inline void clear() { s = 0; lhs=0;}
  inline void set(const Real v) { s = v; }
  inline void copy(const PoissonElement& c) { s = c.s; }

  PoissonElement& operator=(const PoissonElement& c) = default;

  PoissonElement &operator*=(const Real a)
  {
    this->s*=a;
    this->lhs*=a;
    return *this;
  }
  PoissonElement &operator+=(const PoissonElement &rhs)
  {
    this->s+=rhs.s;
    this->lhs+=rhs.s;
    return *this;
  }
  PoissonElement &operator-=(const PoissonElement &rhs)
  {
    this->s-=rhs.s;
    this->lhs-=rhs.s;
    return *this;
  }
  PoissonElement &operator/=(const PoissonElement &rhs)
  {
    this->s/=rhs.s;
    this->lhs/=rhs.s;
    return *this;
  }
  friend PoissonElement operator*(const Real a, PoissonElement el)
  {
      return (el *= a);
  }
  friend PoissonElement operator+(PoissonElement lhs1, const PoissonElement &rhs)
  {
      return (lhs1 += rhs);
  }
  friend PoissonElement operator-(PoissonElement lhs1, const PoissonElement &rhs)
  {
      return (lhs1 -= rhs);
  }
  friend PoissonElement operator/(PoissonElement lhs1, const PoissonElement &rhs)
  {
      return (lhs1 /= rhs);
  }
  bool operator<(const PoissonElement &other) const
  {
     return (s < other.s);
  }
  bool operator>(const PoissonElement &other) const
  {
     return (s > other.s);
  }
  bool operator<=(const PoissonElement &other) const
  {
     return (s <= other.s);
  }
  bool operator>=(const PoissonElement &other) const
  {
     return (s >= other.s);
  }
  double magnitude()
  {
    return s;
  }

  Real & member(int i)
  {
    return (i==0) ? s : lhs;
  }
  static constexpr int DIM = 2;
};

enum { FE_CHI = 0, FE_U, FE_V, FE_W, FE_P, FE_TMPU, FE_TMPV, FE_TMPW };
struct FluidElement
{
  static constexpr int DIM = 8;
  typedef Real RealType;
  Real chi=0, u=0, v=0, w=0, p=0, tmpU=0, tmpV=0, tmpW=0;
  
  void clear() { chi =0; u =0; v =0; w =0; p =0; tmpU =0; tmpV =0; tmpW =0; }

  ~FluidElement() {}

  FluidElement& operator=(const FluidElement& c) = default;

  Real & member(int i)
  {
    Real * tmp = & this->chi;
    return *(tmp + i);
  }
  Real magnitude()//used in TagLoadedBlock, to adapt the mesh
  {
      return chi;
  }
  FluidElement &operator*=(const Real a)
  {
    this->chi  *= a;
    this->u    *= a;
    this->v    *= a;
    this->w    *= a;
    this->p    *= a;
    this->tmpU *= a;
    this->tmpV *= a;
    this->tmpW *= a;
    return *this;
  }
  FluidElement &operator+=(const FluidElement &rhs)
  {
    this->chi  += rhs.chi ;
    this->u    += rhs.u   ;
    this->v    += rhs.v   ;
    this->w    += rhs.w   ;
    this->p    += rhs.p   ;
    this->tmpU += rhs.tmpU;
    this->tmpV += rhs.tmpV;
    this->tmpW += rhs.tmpW;
    return *this;
  }
  FluidElement &operator-=(const FluidElement &rhs)
  {
    this->chi  -= rhs.chi ;
    this->u    -= rhs.u   ;
    this->v    -= rhs.v   ;
    this->w    -= rhs.w   ;
    this->p    -= rhs.p   ;
    this->tmpU -= rhs.tmpU;
    this->tmpV -= rhs.tmpV;
    this->tmpW -= rhs.tmpW;
    return *this;
  }
  //only for debug
  FluidElement &operator=(const double a)
  {
    this->chi  = a;
    this->u    = a;
    this->v    = a;
    this->w    = a;
    this->p    = a;
    this->tmpU = a;
    this->tmpV = a;
    this->tmpW = a;
    return *this;
  }
  friend FluidElement operator*(const Real a, FluidElement el)
  {
    return (el *= a);
  }
  friend FluidElement operator+(FluidElement lhs, const FluidElement &rhs)
  {
    return (lhs += rhs);
  }
  friend FluidElement operator-(FluidElement lhs, const FluidElement &rhs)
  {
    return (lhs -= rhs);
  }
};

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

template <typename TElement>
struct BaseBlock
{
  //these identifiers are required by cubism!
  static constexpr int sizeX = CUP_BLOCK_SIZEX;
  static constexpr int sizeY = CUP_BLOCK_SIZEY;
  static constexpr int sizeZ = CUP_BLOCK_SIZEZ;
  typedef TElement ElementType;
  typedef Real   RealType;
  //__attribute__((aligned(32)))
  TElement data[sizeZ][sizeY][sizeX];
  Real  dataOld[sizeZ][sizeY][sizeX][4];//contains (u,v,w,p) from the previous timestep

  //required from Grid.h
  void clear()
  {
      TElement * entry = &data[0][0][0];
      const int N = sizeX*sizeY*sizeZ;
      for(int i=0; i<N; ++i) entry[i].clear();
      Real * entry1 = &dataOld[0][0][0][0];
      for(int i=0; i<4*N; ++i) entry1[i] = 0.0;
  }

  TElement& operator()(int ix, int iy=0, int iz=0)
  {
    assert(ix>=0); assert(ix<sizeX);
    assert(iy>=0); assert(iy<sizeY);
    assert(iz>=0); assert(iz<sizeZ);
    return data[iz][iy][ix];
  }

  const TElement& operator()(int ix, int iy = 0, int iz = 0) const
  {
    assert(ix>=0); assert(ix<sizeX);
    assert(iy>=0); assert(iy<sizeY);
    assert(iz>=0); assert(iz<sizeZ);
    return data[iz][iy][ix];
  }

  template <typename Streamer>
  inline void Write(std::ofstream& output, Streamer streamer) const
  {
    for(int iz=0; iz<sizeZ; iz++)
      for(int iy=0; iy<sizeY; iy++)
        for(int ix=0; ix<sizeX; ix++)
          streamer.operate(data[iz][iy][ix], output);
  }

  template <typename Streamer>
  inline void Read(std::ifstream& input, Streamer streamer)
  {
    for(int iz=0; iz<sizeZ; iz++)
      for(int iy=0; iy<sizeY; iy++)
        for(int ix=0; ix<sizeX; ix++)
          streamer.operate(input, data[iz][iy][ix]);
  }
};

template <typename TElement>
struct BaseBlockPoisson
{
  static constexpr int sizeX = CUP_BLOCK_SIZEX;
  static constexpr int sizeY = CUP_BLOCK_SIZEY;
  static constexpr int sizeZ = CUP_BLOCK_SIZEZ;
  typedef TElement ElementType;
  typedef Real RealType;
  TElement data[sizeZ][sizeY][sizeX];

  void clear()
  {
      TElement * entry = &data[0][0][0];
      const int N = sizeX*sizeY*sizeZ;
      for(int i=0; i<N; ++i) entry[i].clear();
  }

  TElement& operator()(int ix, int iy=0, int iz=0)
  {
    assert(ix>=0); assert(ix<sizeX);
    assert(iy>=0); assert(iy<sizeY);
    assert(iz>=0); assert(iz<sizeZ);
    return data[iz][iy][ix];
  }

  const TElement& operator()(int ix, int iy = 0, int iz = 0) const
  {
    assert(ix>=0); assert(ix<sizeX);
    assert(iy>=0); assert(iy<sizeY);
    assert(iz>=0); assert(iz<sizeZ);
    return data[iz][iy][ix];
  }
};

struct StreamerChi
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).chi;
    }
    static std::string prefix() {return std::string("chi_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerVelocityVector
{
    static const int NCHANNELS = 3;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
        output[0] = b(ix,iy,iz).u; output[1] = b(ix,iy,iz).v; output[2] = b(ix,iy,iz).w;
    }
    static std::string prefix(){return std::string("vel_");}
    static const char * getAttributeName() { return "Vector"; }
};

struct StreamerTmpVector
{
    static const int NCHANNELS = 3;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
        output[0] = b(ix,iy,iz).tmpU; output[1] = b(ix,iy,iz).tmpV; output[2] = b(ix,iy,iz).tmpW;
    }
    static std::string prefix(){return std::string("tmp_");}
    static const char * getAttributeName() { return "Vector"; }
};

struct StreamerPressure
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).p;
    }
    static std::string prefix(){return std::string("pres_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerTmpVectorX
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).tmpU;
    }
    static std::string prefix(){return std::string("tmpU_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerTmpVectorY
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).tmpV;
    }
    static std::string prefix(){return std::string("tmpV_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerTmpVectorZ
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).tmpW;
    }
    static std::string prefix(){return std::string("tmpW_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerVelVectorX
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).u;
    }
    static std::string prefix(){return std::string("u_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerVelVectorY
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).v;
    }
    static std::string prefix(){return std::string("v_");}
    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerVelVectorZ
{
    static const int NCHANNELS = 1;
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).w;
    }
    static std::string prefix(){return std::string("w_");}
    static const char * getAttributeName() { return "Scalar"; }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBC: public cubism::BlockLab<BlockType,allocator>
{
  typedef typename BlockType::ElementType ElementTypeBlock;
  static constexpr int sizeX = BlockType::sizeX;
  static constexpr int sizeY = BlockType::sizeY;
  static constexpr int sizeZ = BlockType::sizeZ;

  BCflag BCX = freespace;
  BCflag BCY = freespace;
  BCflag BCZ = freespace;

  // Used for Boundary Conditions:
  // Apply bc on face of direction dir and side side (0 or 1):
  template<int dir, int side> void applyBCfaceOpen(const bool coarse = false)
  {
    //First we apply zero Neumann BCs to all variables. Then, if we are dealing with a 
    //FluidElement (which has DIM=8), we set the normal velocity to zero.
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
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]) = cb->Access
          (
            ( dir==0 ? (side==0 ? 0 : sizeX-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ-1 ) : iz ) - stenBeg[2]
          );      
      if (ElementTypeBlock::DIM==8) // (u,v,w)*(nx,ny,nz) = 0
      for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
      for(int ix=s[0]; ix<e[0]; ix++)
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(1+dir) = (-1.)*cb->Access
          (
            ( dir==0 ? (side==0 ? 0 : sizeX-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ-1 ) : iz ) - stenBeg[2]
          ).member(1+dir);      

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
        }
      }

      if (ElementTypeBlock::DIM==8) // (u,v,w)*(nx,ny,nz) = 0
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
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(1+dir) = dir==0?
          (-1.)*cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(1+dir) : (dir==1?
          (-1.)*cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(1+dir) :
          (-1.)*cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(1+dir));
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
      }

      if (ElementTypeBlock::DIM==8) // (u,v,w)*(nx,ny,nz) = 0
      for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
      for(int ix=s[0]; ix<e[0]; ix++)
      {
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(1+dir) = (-1.)*cb->Access
          ( ( dir==0 ? (side==0 ? 0 : sizeX/2-1 ) : ix ) - stenBeg[0],
            ( dir==1 ? (side==0 ? 0 : sizeY/2-1 ) : iy ) - stenBeg[1],
            ( dir==2 ? (side==0 ? 0 : sizeZ/2-1 ) : iz ) - stenBeg[2]).member(1+dir);
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
        }
      }

      if (ElementTypeBlock::DIM==8) // (u,v,w)*(nx,ny,nz) = 0
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
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(1+dir) = dir==0?
          (-1.)*cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(1+dir) : (dir==1?
          (-1.)*cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(1+dir) :
          (-1.)*cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(1+dir));
        }
      }

    }
  }
  template<int dir, int side> void applyBCfaceWall(const bool coarse=false)
  {
    assert (ElementTypeBlock::DIM==FluidElement::DIM);
    const int mask [8] = {-1,-1,-1,-1,+1,-1,-1,-1};//apply zero Dirichlet to all variables except pressure
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
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = mask[k]*cb->Access
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
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = mask[k] * dir==0?
          cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) : (dir==1?
          cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(k):
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(k));
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
        cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = mask[k]*cb->Access
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
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) = mask[k] * dir==0?
          cb->Access(side*(bsize[0]-1)-stenBeg[0], iy-stenBeg[1], iz-stenBeg[2]).member(k) : (dir==1?
          cb->Access(ix-stenBeg[0], side*(bsize[1]-1)-stenBeg[1], iz-stenBeg[2]).member(k) :
          cb->Access(ix-stenBeg[0], iy-stenBeg[1], side*(bsize[2]-1)-stenBeg[2]).member(k));
        }
      }
    }
  }

 public:

  typedef typename BlockType::ElementType ElementType;
  void setBC(const BCflag _BCX, const BCflag _BCY, const BCflag _BCZ) {
    BCX=_BCX; BCY=_BCY; BCZ=_BCZ;
  }
  bool is_xperiodic() { return BCX == periodic; }
  bool is_yperiodic() { return BCY == periodic; }
  bool is_zperiodic() { return BCZ == periodic; }

  BlockLabBC() = default;
  BlockLabBC(const BlockLabBC&) = delete;
  BlockLabBC& operator=(const BlockLabBC&) = delete;

  // Called by Cubism:
  // wall means u=v=w=0 and dp/dn=0
  // freespace means normal velocity=0, dp/dn=0 and d(tangential velocities)/dn=0
  void _apply_bc(const cubism::BlockInfo& info, const Real t=0, const bool coarse = false)
  {
    if (BCX == wall && ElementTypeBlock::DIM==FluidElement::DIM)//wall only makes sense for FluidElement
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
    if (BCY == wall && ElementTypeBlock::DIM==FluidElement::DIM)//wall only makes sense for FluidElement
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
    if (BCZ == wall && ElementTypeBlock::DIM==FluidElement::DIM)//wall only makes sense for FluidElement
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

using FluidBlock = BaseBlock<FluidElement>;
using FluidGrid    = cubism::Grid<FluidBlock, aligned_allocator>;
using FluidGridMPI = cubism::GridMPI<FluidGrid>;
using Lab          = BlockLabBC<FluidBlock, aligned_allocator>;
using LabMPI       = cubism::BlockLabMPI<Lab,FluidGridMPI>;

using AMR = MeshAdaptation_CUP<FluidGridMPI,LabMPI>;

using FluidBlockPoisson = BaseBlockPoisson<PoissonElement>;
using FluidGridPoisson  = cubism::Grid<FluidBlockPoisson, aligned_allocator>;
using FluidGridMPIPoisson = cubism::GridMPI<FluidGridPoisson>;
using LabPoisson          = BlockLabBC<FluidBlockPoisson, aligned_allocator>;
using LabMPIPoisson       = cubism::BlockLabMPI<LabPoisson,FluidGridMPIPoisson>;
using AMR2 = MeshAdaptationMPI<FluidGridMPIPoisson,LabMPIPoisson,FluidGridMPI>;

CubismUP_3D_NAMESPACE_END
