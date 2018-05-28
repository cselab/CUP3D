//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Christian Conti
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_DataStructures_h
#define CubismUP_3D_DataStructures_h

//#include <cassert>
//#define __2Leads_
//#define __DumpWakeStefan 9
#define __useSkin_
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <math.h>
#include <string>
#include <vector>
#include <array>
using namespace std;

#define DUMPGRID

#include <mpi.h>
#include <omp.h>

#ifndef _FLOAT_PRECISION_
typedef double Real;
#else // _FLOAT_PRECISION_
typedef float Real;
#endif // _FLOAT_PRECISION_

//this is all cubism file we need
#include <ArgumentParser.h>
#include <AlignedAllocator.h>
#include <Grid.h>
#include <GridMPI.h>
#include <BlockInfo.h>
#ifdef _VTK_
#include <SerializerIO_ImageVTK.h>
#endif
//#include <HDF5Dumper_MPI.h>
//#include <ZBinDumper_MPI.h>
#include <BlockLab.h>
#include <BlockLabMPI.h>
#include <Profiler.h>
#include <StencilInfo.h>
#include "Timer.h"
#include "BoundaryConditions.h"
#include "ObstacleBlock.h"
#include "StateRewardData.h"
#include "Slice.h"


struct FluidElement
{
  //If you modify these (adding, moving, shuffling, whatever) you kill the code
  Real chi, u, v, w, p, tmpU, tmpV, tmpW;
  FluidElement(): chi(0),u(0),v(0),w(0),p(0),tmpU(0),tmpV(0),tmpW(0) {}
  void clear() { chi =0; u =0; v =0; w =0; p =0; tmpU =0; tmpV =0; tmpW =0; }
};

struct DumpElement {
    float u, v, w, chi, p;
    DumpElement() : u(0), v(0), w(0), chi(0), p(0) {}
    void clear() { u = v = w = chi = p = 0; }
};

struct StreamerDiv
{
  static const int channels = 1;
  static void operate(const FluidElement& input, Real output[1])
  { output[0] = input.p; }

  static void operate(const Real input[1], FluidElement& output)
  { output.p = input[0]; }
};

template <typename TElement>
struct BaseBlock
{
  //these identifiers are required by cubism!
  static constexpr int sizeX = _BLOCKSIZEX_;
  static constexpr int sizeY = _BLOCKSIZEY_;
  static constexpr int sizeZ = _BLOCKSIZEZ_;
  static constexpr int sizeArray[3] = {_BLOCKSIZEX_, _BLOCKSIZEY_, _BLOCKSIZEZ_};
  typedef TElement ElementType;
  typedef TElement element_type;
  __attribute__((aligned(32))) TElement data[sizeZ][sizeY][sizeX];

  //required from Grid.h
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

  template <typename Streamer>
  inline void Write(ofstream& output, Streamer streamer) const
  {
    for(int iz=0; iz<sizeZ; iz++)
      for(int iy=0; iy<sizeY; iy++)
        for(int ix=0; ix<sizeX; ix++)
          streamer.operate(data[iz][iy][ix], output);
  }

  template <typename Streamer>
  inline void Read(ifstream& input, Streamer streamer)
  {
    for(int iz=0; iz<sizeZ; iz++)
      for(int iy=0; iy<sizeY; iy++)
        for(int ix=0; ix<sizeX; ix++)
          streamer.operate(input, data[iz][iy][ix]);
  }
};

struct StreamerChi
{
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename T>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
      output[0] = b(ix,iy,iz).chi;
    }
    static std::string postfix()
    {
      return std::string("-chi");
    }
    static std::string descriptor()
    {
      return std::string("Chi field");
    }

    static const char * getAttributeName() { return "Scalar"; }
};

struct StreamerVelocityVector
{
    static const int NCHANNELS = 3;
    static const int CLASS = 0;

    // Write
    template <typename TBlock, typename T>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
        output[0] = b(ix,iy,iz).u;
        output[1] = b(ix,iy,iz).v;
        output[2] = b(ix,iy,iz).w;
    }

    // Read
    template <typename TBlock, typename T>
    static inline void operate(TBlock& b, const T input[NCHANNELS], const int ix, const int iy, const int iz)
    {
        b(ix,iy,iz).u = input[0];
        b(ix,iy,iz).v = input[1];
        b(ix,iy,iz).w = input[2];
    }

    static std::string postfix() {
      return std::string("-vel");
    }
    static std::string descriptor() {
      return std::string("Velocity vector");
    }

    static const char * getAttributeName() { return "Vector"; }
};

struct StreamerPressure
{
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename T>
    inline void operate(const TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS]) const
    {
      output[0] = b(ix,iy,iz).p;
    }
    static std::string postfix() {
      return std::string("-pre");
    }
    static std::string descriptor() {
      return std::string("Pressure field");
    }

    static const char * getAttributeName() { return "Scalar"; }
};


template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabOpen: public BlockLab<BlockType,allocator>
{
  typedef typename BlockType::ElementType ElementTypeBlock;
 public:
  //virtual inline std::string name() const { return "BlockLabOpen"; }
  bool is_xperiodic() {return false;}
  bool is_yperiodic() {return false;}
  bool is_zperiodic() {return false;}
  BlockLabOpen(): BlockLab<BlockType,allocator>(){}
  void _apply_bc(const BlockInfo& info, const Real t=0)
  {
    BoundaryCondition<BlockType,ElementTypeBlock,allocator>
      bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

    if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
    //if (info.index[0]==0)           bc.template applyBC_dirichlet_inflow<0,0>();
    if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
    if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
    if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
    if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
    if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
  }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabPeriodicZ : public BlockLab<BlockType,allocator>
{
  typedef typename BlockType::ElementType ElementTypeBlock;
 public:
  bool is_xperiodic() {return false;}
  bool is_yperiodic() {return false;}
  bool is_zperiodic() {return true; }
  BlockLabPeriodicZ(): BlockLab<BlockType,allocator>(){}

  void _apply_bc(const BlockInfo& info, const Real t=0)
  {
    BoundaryCondition<BlockType,ElementTypeBlock,allocator>
      bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

    if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
    //if (info.index[0]==0)           bc.template applyBC_dirichlet_inflow<0,0>();
    if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
    if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
    if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
  }
};

using FluidBlock = BaseBlock<FluidElement>;
typedef Grid<FluidBlock, aligned_allocator> FluidGrid;
typedef GridMPI<FluidGrid> FluidGridMPI;

#ifdef DUMPGRID
  using DumpBlock = BaseBlock<DumpElement>;
  typedef GridMPI<Grid<DumpBlock, aligned_allocator>> DumpGridMPI;
  typedef SliceMPI<DumpGridMPI> SliceType;
#else
  typedef SliceMPI<FluidGridMPI> SliceType;
#endif

#if defined(BC_PERIODICZ)
  typedef  BlockLabPeriodicZ<FluidBlock, aligned_allocator> Lab;
#else
  typedef  BlockLabOpen<FluidBlock, aligned_allocator> Lab;
#endif

typedef BlockLabMPI<Lab> LabMPI;

//#ifdef _VORTEX_
//typedef BlockLabVortex<FluidBlock, std::allocator> Lab;
//#endif // _VORTEX_

//#ifdef _PIPE_
//typedef BlockLabPipe<FluidBlock, std::allocator> Lab;
//#endif // _PIPE_

#endif
