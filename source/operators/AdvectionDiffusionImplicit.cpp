//
//  Cubism3D
//  Copyright (c) 2021 CSE-Lab, ETH Zurich, Switzerland.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "AdvectionDiffusionImplicit.h"

CubismUP_3D_NAMESPACE_BEGIN

//#define WENO
#ifdef PRESERVE_SYMMETRY
#define DISABLE_OPTIMIZATIONS __attribute__((optimize("-O1")))
#else
#define DISABLE_OPTIMIZATIONS
#endif

struct KernelDiffusionRHS
{
  SimulationData & sim;
  StencilInfo stencil  = StencilInfo(-1,-1,-1, 2,2,2, false, {0,1,2});
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  const std::vector<BlockInfo> &  tmpVInfo = sim.tmpVInfo();

  KernelDiffusionRHS(SimulationData& s) :sim(s){}

  void operator()(const VectorLab & lab, const BlockInfo& info) const
  {
    VectorBlock & __restrict__ TMPV  = (*sim.tmpV)(info.blockID);
    const Real facD = info.h;

    for(int z=0; z<Nz; ++z)
    for(int y=0; y<Ny; ++y)
    for(int x=0; x<Nx; ++x)
    {
      const Real duD =( (lab(x+1,y,z).u[0] + lab(x-1,y,z).u[0]) + ((lab(x,y+1,z).u[0] + lab(x,y-1,z).u[0]) + (lab(x,y,z+1).u[0] + lab(x,y,z-1).u[0])) )- 6 * lab(x,y,z).u[0];
      const Real dvD =( (lab(x,y+1,z).u[1] + lab(x,y-1,z).u[1]) + ((lab(x,y,z+1).u[1] + lab(x,y,z-1).u[1]) + (lab(x+1,y,z).u[1] + lab(x-1,y,z).u[1])) )- 6 * lab(x,y,z).u[1];
      const Real dwD =( (lab(x,y,z+1).u[2] + lab(x,y,z-1).u[2]) + ((lab(x+1,y,z).u[2] + lab(x-1,y,z).u[2]) + (lab(x,y+1,z).u[2] + lab(x,y-1,z).u[2])) )- 6 * lab(x,y,z).u[2];
      TMPV(x,y,z).u[0] = facD * duD;
      TMPV(x,y,z).u[1] = facD * dvD;
      TMPV(x,y,z).u[2] = facD * dwD;
    }

        BlockCase<VectorBlock> * tempCase = (BlockCase<VectorBlock> *)(tmpVInfo[info.blockID].auxiliary);

        if (tempCase == nullptr) return; //no flux corrections needed for this block

        VectorElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
        VectorElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
        VectorElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
        VectorElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
        VectorElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
        VectorElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
        if (faceXm != nullptr)
        {
          const int x = 0;
          for(int z=0; z<Nz; ++z)
          for(int y=0; y<Ny; ++y)
          {
            faceXm[y + Ny * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x-1,y,z).u[0]);
            faceXm[y + Ny * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x-1,y,z).u[1]);
            faceXm[y + Ny * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x-1,y,z).u[2]);
          }
        }
        if (faceXp != nullptr)
        {
          const int x = Nx-1;
          for(int z=0; z<Nz; ++z)
          for(int y=0; y<Ny; ++y)
          {
            faceXp[y + Ny * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x+1,y,z).u[0]);
            faceXp[y + Ny * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x+1,y,z).u[1]);
            faceXp[y + Ny * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x+1,y,z).u[2]);
          }
        }
        if (faceYm != nullptr)
        {
          const int y = 0;
          for(int z=0; z<Nz; ++z)
          for(int x=0; x<Nx; ++x)
          {
            faceYm[x + Nx * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y-1,z).u[0]);
            faceYm[x + Nx * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y-1,z).u[1]);
            faceYm[x + Nx * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y-1,z).u[2]);
          }
        }
        if (faceYp != nullptr)
        {
          const int y = Ny-1;
          for(int z=0; z<Nz; ++z)
          for(int x=0; x<Nx; ++x)
          {
            faceYp[x + Nx * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y+1,z).u[0]);
            faceYp[x + Nx * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y+1,z).u[1]);
            faceYp[x + Nx * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y+1,z).u[2]);
          }
        }
        if (faceZm != nullptr)
        {
          const int z = 0;
          for(int y=0; y<Ny; ++y)
          for(int x=0; x<Nx; ++x)
          {
            faceZm[x + Nx * y].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y,z-1).u[0]);
            faceZm[x + Nx * y].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y,z-1).u[1]);
            faceZm[x + Nx * y].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y,z-1).u[2]);
          }
        }
        if (faceZp != nullptr)
        {
          const int z = Nz-1;
          for(int y=0; y<Ny; ++y)
          for(int x=0; x<Nx; ++x)
          {
            faceZp[x + Nx * y].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y,z+1).u[0]);
            faceZp[x + Nx * y].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y,z+1).u[1]);
            faceZp[x + Nx * y].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y,z+1).u[2]);
          }
        }

  }
};

struct KernelAdvect
{
    KernelAdvect(const SimulationData&s, const Real _dt) : sim(s), dt(_dt){}
    const SimulationData & sim;
    const Real dt;
    const Real mu = sim.nu;
    const std::array<Real, 3>& uInf = sim.uinf;
    const std::vector<BlockInfo> &tmpVInfo = sim.tmpVInfo();
    const StencilInfo stencil{-3,-3,-3,4,4,4,false, {0,1,2}};
    const int Nx = VectorBlock::sizeX;
    const int Ny = VectorBlock::sizeY;
    const int Nz = VectorBlock::sizeZ;
  #ifdef WENO
    DISABLE_OPTIMIZATIONS
    inline Real weno5_plus(const Real & um2, const Real & um1, const Real & u, const Real & up1, const Real & up2) const
    {
      const Real exponent = 2;
      const Real e = 1e-6;
      const Real b1 = 13.0/12.0*pow((um2+u)-2*um1,2)+0.25*pow((um2+3*u)-4*um1,2);
      const Real b2 = 13.0/12.0*pow((um1+up1)-2*u,2)+0.25*pow(um1-up1,2);
      const Real b3 = 13.0/12.0*pow((u+up2)-2*up1,2)+0.25*pow((3*u+up2)-4*up1,2);
      const Real g1 = 0.1;
      const Real g2 = 0.6;
      const Real g3 = 0.3;
      const Real what1 = g1/pow(b1+e,exponent);
      const Real what2 = g2/pow(b2+e,exponent);
      const Real what3 = g3/pow(b3+e,exponent);
      const Real aux = 1.0/((what1+what3)+what2);
      const Real w1 = what1*aux;
      const Real w2 = what2*aux;
      const Real w3 = what3*aux;
      const Real f1 = (11.0/6.0)*u + ( ( 1.0/3.0)*um2- (7.0/6.0)*um1);
      const Real f2 = (5.0 /6.0)*u + ( (-1.0/6.0)*um1+ (1.0/3.0)*up1);
      const Real f3 = (1.0 /3.0)*u + ( (+5.0/6.0)*up1- (1.0/6.0)*up2);
      return (w1*f1+w3*f3)+w2*f2;
    }
    DISABLE_OPTIMIZATIONS
    inline Real weno5_minus(const Real & um2, const Real & um1, const Real & u, const Real & up1, const Real & up2) const
    {
      const Real exponent = 2;
      const Real e = 1e-6;
      const Real b1 = 13.0/12.0*pow((um2+u)-2*um1,2)+0.25*pow((um2+3*u)-4*um1,2);
      const Real b2 = 13.0/12.0*pow((um1+up1)-2*u,2)+0.25*pow(um1-up1,2);
      const Real b3 = 13.0/12.0*pow((u+up2)-2*up1,2)+0.25*pow((3*u+up2)-4*up1,2);
      const Real g1 = 0.3;
      const Real g2 = 0.6;
      const Real g3 = 0.1;
      const Real what1 = g1/pow(b1+e,exponent);
      const Real what2 = g2/pow(b2+e,exponent);
      const Real what3 = g3/pow(b3+e,exponent);
      const Real aux = 1.0/((what1+what3)+what2);
      const Real w1 = what1*aux;
      const Real w2 = what2*aux;
      const Real w3 = what3*aux;
      const Real f1 = ( 1.0/3.0)*u + ( (-1.0/6.0)*um2+ (5.0/6.0)*um1);
      const Real f2 = ( 5.0/6.0)*u + ( ( 1.0/3.0)*um1- (1.0/6.0)*up1);
      const Real f3 = (11.0/6.0)*u + ( (-7.0/6.0)*up1+ (1.0/3.0)*up2);
      return (w1*f1+w3*f3)+w2*f2;
    }
    inline Real derivative(const Real & U, const Real & um3, const Real & um2, const Real & um1,
                           const Real & u, const Real & up1, const Real & up2, const Real & up3) const
    {
      Real fp = 0.0;
      Real fm = 0.0;
      if (U > 0)
      {
        fp = weno5_plus (um2,um1,u,up1,up2);
        fm = weno5_plus (um3,um2,um1,u,up1);
      }
      else
      {
        fp = weno5_minus(um1,u,up1,up2,up3);
        fm = weno5_minus(um2,um1,u,up1,up2);
      }
      return (fp-fm);
    }
  #else
    inline Real derivative(const Real & U, const Real & um3, const Real & um2, const Real & um1,
                           const Real & u, const Real & up1, const Real & up2, const Real & up3) const
    {
      if (U > 0) return(-2*um3+15*um2-60*um1+20*u+30*up1-3*up2)/60.;
      else       return( 2*up3-15*up2+60*up1-20*u-30*um1+3*um2)/60.;
    }
  #endif

    void operator()(const VectorLab & lab, const BlockInfo& info) const
    {
        VectorBlock & o = (*sim.tmpV)(info.blockID);
        VectorBlock & v = (*sim.vel )(info.blockID);

        const Real h3   = info.h*info.h*info.h;
        const Real facA = -dt/info.h * h3;
        const Real facD = (mu/info.h)*(dt/info.h) * h3;

        for (int z=0; z<Nz; ++z)
        for (int y=0; y<Ny; ++y)
        for (int x=0; x<Nx; ++x)
        {
            const Real uAbs[3] = { lab(x,y,z).u[0] + uInf[0], lab(x,y,z).u[1] + uInf[1], lab(x,y,z).u[2] + uInf[2] };
            const Real dudx = derivative(uAbs[0],lab(x-3,y,z).u[0],lab(x-2,y,z).u[0],lab(x-1,y,z).u[0],lab(x,y,z).u[0],lab(x+1,y,z).u[0],lab(x+2,y,z).u[0],lab(x+3,y,z).u[0]);
            const Real dvdx = derivative(uAbs[0],lab(x-3,y,z).u[1],lab(x-2,y,z).u[1],lab(x-1,y,z).u[1],lab(x,y,z).u[1],lab(x+1,y,z).u[1],lab(x+2,y,z).u[1],lab(x+3,y,z).u[1]);
            const Real dwdx = derivative(uAbs[0],lab(x-3,y,z).u[2],lab(x-2,y,z).u[2],lab(x-1,y,z).u[2],lab(x,y,z).u[2],lab(x+1,y,z).u[2],lab(x+2,y,z).u[2],lab(x+3,y,z).u[2]);
            const Real dudy = derivative(uAbs[1],lab(x,y-3,z).u[0],lab(x,y-2,z).u[0],lab(x,y-1,z).u[0],lab(x,y,z).u[0],lab(x,y+1,z).u[0],lab(x,y+2,z).u[0],lab(x,y+3,z).u[0]);
            const Real dvdy = derivative(uAbs[1],lab(x,y-3,z).u[1],lab(x,y-2,z).u[1],lab(x,y-1,z).u[1],lab(x,y,z).u[1],lab(x,y+1,z).u[1],lab(x,y+2,z).u[1],lab(x,y+3,z).u[1]);
            const Real dwdy = derivative(uAbs[1],lab(x,y-3,z).u[2],lab(x,y-2,z).u[2],lab(x,y-1,z).u[2],lab(x,y,z).u[2],lab(x,y+1,z).u[2],lab(x,y+2,z).u[2],lab(x,y+3,z).u[2]);
            const Real dudz = derivative(uAbs[2],lab(x,y,z-3).u[0],lab(x,y,z-2).u[0],lab(x,y,z-1).u[0],lab(x,y,z).u[0],lab(x,y,z+1).u[0],lab(x,y,z+2).u[0],lab(x,y,z+3).u[0]);
            const Real dvdz = derivative(uAbs[2],lab(x,y,z-3).u[1],lab(x,y,z-2).u[1],lab(x,y,z-1).u[1],lab(x,y,z).u[1],lab(x,y,z+1).u[1],lab(x,y,z+2).u[1],lab(x,y,z+3).u[1]);
            const Real dwdz = derivative(uAbs[2],lab(x,y,z-3).u[2],lab(x,y,z-2).u[2],lab(x,y,z-1).u[2],lab(x,y,z).u[2],lab(x,y,z+1).u[2],lab(x,y,z+2).u[2],lab(x,y,z+3).u[2]);
            const Real duD =( (lab(x+1,y,z).u[0] + lab(x-1,y,z).u[0]) + ((lab(x,y+1,z).u[0] + lab(x,y-1,z).u[0]) + (lab(x,y,z+1).u[0] + lab(x,y,z-1).u[0])) )- 6 * lab(x,y,z).u[0];
            const Real dvD =( (lab(x,y+1,z).u[1] + lab(x,y-1,z).u[1]) + ((lab(x,y,z+1).u[1] + lab(x,y,z-1).u[1]) + (lab(x+1,y,z).u[1] + lab(x-1,y,z).u[1])) )- 6 * lab(x,y,z).u[1];
            const Real dwD =( (lab(x,y,z+1).u[2] + lab(x,y,z-1).u[2]) + ((lab(x+1,y,z).u[2] + lab(x-1,y,z).u[2]) + (lab(x,y+1,z).u[2] + lab(x,y-1,z).u[2])) )- 6 * lab(x,y,z).u[2];
            const Real duA = uAbs[0] * dudx + (uAbs[1] * dudy + uAbs[2] * dudz);
            const Real dvA = uAbs[1] * dvdy + (uAbs[2] * dvdz + uAbs[0] * dvdx);
            const Real dwA = uAbs[2] * dwdz + (uAbs[0] * dwdx + uAbs[1] * dwdy);
            o(x,y,z).u[0] =  facD*duD ;
            o(x,y,z).u[1] =  facD*dvD ;
            o(x,y,z).u[2] =  facD*dwD ;
            v(x,y,z).u[0] += facA*duA/h3;
            v(x,y,z).u[1] += facA*dvA/h3;
            v(x,y,z).u[2] += facA*dwA/h3;
        }

        BlockCase<VectorBlock> * tempCase = (BlockCase<VectorBlock> *)(tmpVInfo[info.blockID].auxiliary);

        if (tempCase == nullptr) return; //no flux corrections needed for this block

        VectorElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
        VectorElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
        VectorElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
        VectorElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
        VectorElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
        VectorElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
        if (faceXm != nullptr)
        {
          const int x = 0;
          for(int z=0; z<Nz; ++z)
          for(int y=0; y<Ny; ++y)
          {
            faceXm[y + Ny * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x-1,y,z).u[0]);
            faceXm[y + Ny * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x-1,y,z).u[1]);
            faceXm[y + Ny * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x-1,y,z).u[2]);
          }
        }
        if (faceXp != nullptr)
        {
          const int x = Nx-1;
          for(int z=0; z<Nz; ++z)
          for(int y=0; y<Ny; ++y)
          {
            faceXp[y + Ny * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x+1,y,z).u[0]);
            faceXp[y + Ny * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x+1,y,z).u[1]);
            faceXp[y + Ny * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x+1,y,z).u[2]);
          }
        }
        if (faceYm != nullptr)
        {
          const int y = 0;
          for(int z=0; z<Nz; ++z)
          for(int x=0; x<Nx; ++x)
          {
            faceYm[x + Nx * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y-1,z).u[0]);
            faceYm[x + Nx * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y-1,z).u[1]);
            faceYm[x + Nx * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y-1,z).u[2]);
          }
        }
        if (faceYp != nullptr)
        {
          const int y = Ny-1;
          for(int z=0; z<Nz; ++z)
          for(int x=0; x<Nx; ++x)
          {
            faceYp[x + Nx * z].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y+1,z).u[0]);
            faceYp[x + Nx * z].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y+1,z).u[1]);
            faceYp[x + Nx * z].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y+1,z).u[2]);
          }
        }
        if (faceZm != nullptr)
        {
          const int z = 0;
          for(int y=0; y<Ny; ++y)
          for(int x=0; x<Nx; ++x)
          {
            faceZm[x + Nx * y].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y,z-1).u[0]);
            faceZm[x + Nx * y].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y,z-1).u[1]);
            faceZm[x + Nx * y].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y,z-1).u[2]);
          }
        }
        if (faceZp != nullptr)
        {
          const int z = Nz-1;
          for(int y=0; y<Ny; ++y)
          for(int x=0; x<Nx; ++x)
          {
            faceZp[x + Nx * y].u[0] = facD*(lab(x,y,z).u[0] - lab(x,y,z+1).u[0]);
            faceZp[x + Nx * y].u[1] = facD*(lab(x,y,z).u[1] - lab(x,y,z+1).u[1]);
            faceZp[x + Nx * y].u[2] = facD*(lab(x,y,z).u[2] - lab(x,y,z+1).u[2]);
          }
        }
    }
};

void AdvectionDiffusionImplicit::euler(const Real dt)
{
  const std::vector<BlockInfo> &  velInfo = sim.velInfo();
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  const size_t Nblocks = velInfo.size();

  pressure.resize(Nblocks*Nx*Ny*Nz);
  velocity.resize(Nblocks*Nx*Ny*Nz*3);

  //Explicit Euler timestep for advection terms. We also store pressure field.
  compute<VectorLab>(KernelAdvect(sim,dt),sim.vel,sim.tmpV);
  #pragma omp parallel for
  for(size_t i=0; i<Nblocks; i++)
  {
    const VectorBlock & TMPV = (*sim.tmpV)(i);
    const ScalarBlock & P    = (*sim.pres)(i);
    VectorBlock & V          = (*sim.vel )(i);
    const Real ih3 = 1.0/(velInfo[i].h*velInfo[i].h*velInfo[i].h);
    for (int z=0; z<Nz; ++z)
    for (int y=0; y<Ny; ++y)
    for (int x=0; x<Nx; ++x)
    {
      const int idx = i*Nx*Ny*Nz+z*Ny*Nx+y*Nx+x;
      pressure[idx] = P(x,y,z).s;
      velocity[3*idx+0] = V(x,y,z).u[0];
      velocity[3*idx+1] = V(x,y,z).u[1];
      velocity[3*idx+2] = V(x,y,z).u[2];
      V(x,y,z).u[0] = TMPV(x,y,z).u[0]*ih3 + V(x,y,z).u[0];
      V(x,y,z).u[1] = TMPV(x,y,z).u[1]*ih3 + V(x,y,z).u[1];
      V(x,y,z).u[2] = TMPV(x,y,z).u[2]*ih3 + V(x,y,z).u[2];
    }
  }

  compute<VectorLab>(KernelDiffusionRHS(sim),sim.vel,sim.tmpV);

  #pragma omp parallel for
  for(size_t i=0; i<Nblocks; i++)
  {
    VectorBlock & V    = (*sim.vel )(i);
    VectorBlock & TMPV = (*sim.tmpV)(i);
    const Real ih3 = 1.0/(velInfo[i].h*velInfo[i].h*velInfo[i].h);
    for (int z=0; z<Nz; ++z)
    for (int y=0; y<Ny; ++y)
    for (int x=0; x<Nx; ++x)
    {
      const int idx = i*Nx*Ny*Nz+z*Ny*Nx+y*Nx+x;
      TMPV(x,y,z).u[0] = -TMPV(x,y,z).u[0]*ih3 + (V(x,y,z).u[0]- velocity[3*idx+0])/(dt*sim.nu);
      TMPV(x,y,z).u[1] = -TMPV(x,y,z).u[1]*ih3 + (V(x,y,z).u[1]- velocity[3*idx+1])/(dt*sim.nu);
      TMPV(x,y,z).u[2] = -TMPV(x,y,z).u[2]*ih3 + (V(x,y,z).u[2]- velocity[3*idx+2])/(dt*sim.nu);
    }
  }


  DiffusionSolver mmysolver(sim);
  for (int index = 0 ; index < 3; index ++)
  {
    mmysolver.mydirection = index;
    mmysolver.dt = dt;
    #pragma omp parallel for
    for(size_t i=0; i<Nblocks; i++)//Set u^{n+1/2} = u^{n} + 0.5*dt*RHS(u^{n})
    {
        const Real h3 = (velInfo[i].h*velInfo[i].h*velInfo[i].h);
        ScalarBlock & RHS = (*sim.lhs)(i);
        ScalarBlock & P   = (*sim.pres)(i);
        const VectorBlock & TMPV = (*sim.tmpV)(i);
        for (int z=0; z<Nz; ++z)
        for (int y=0; y<Ny; ++y)
        for (int x=0; x<Nx; ++x)
        {
          P(x,y,z).s = 0;
          RHS(x,y,z).s = h3*TMPV(x,y,z).u[index];
        }
    }
    mmysolver.solve();
    #pragma omp parallel for
    for(size_t i=0; i<Nblocks; i++)//Set u^{n+1/2} = u^{n} + 0.5*dt*RHS(u^{n})
    {
      ScalarBlock & P   = (*sim.pres)(i);
      VectorBlock & V = (*sim.vel )(i);
      for (int z=0; z<Nz; ++z)
      for (int y=0; y<Ny; ++y)
      for (int x=0; x<Nx; ++x)
      {
        V(x,y,z).u[index] += P(x,y,z).s;
      }
    }  
  }

  #pragma omp parallel for
  for(size_t i=0; i<Nblocks; i++)
  {
    ScalarBlock & P    = (*sim.pres)(i);
    for (int z=0; z<Nz; ++z)
    for (int y=0; y<Ny; ++y)
    for (int x=0; x<Nx; ++x)
    {
      const int idx = i*Nx*Ny*Nz+z*Ny*Nx+y*Nx+x;
      P(x,y,z).s = pressure[idx];
    }
  }
}

void AdvectionDiffusionImplicit::operator()(const Real dt)
{
  euler(sim.dt);
}

CubismUP_3D_NAMESPACE_END
