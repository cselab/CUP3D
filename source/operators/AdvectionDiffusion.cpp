//
//  Cubism3D
//  Copyright (c) 2021 CSE-Lab, ETH Zurich, Switzerland.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "AdvectionDiffusion.h"
#include "../Obstacles/ObstacleVector.h"
#include "../poisson/PoissonSolverAMR.h"
CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

namespace {

struct KernelAdvectDiffuse
{
    KernelAdvectDiffuse(const SimulationData&s, const Real a_coef) : sim(s), coef(a_coef) {}
    const SimulationData & sim;
    const Real dt = sim.dt;
    const Real mu = sim.nu;
    const Real coef;
    const std::array<Real, 3>& uInf = sim.uinf;
    const StencilInfo stencil{-3,-3,-3,4,4,4,false, {FE_U, FE_V, FE_W}};

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


    void operator()(LabMPI & lab, const BlockInfo& info) const
    {
        FluidBlock& o = *(FluidBlock*)info.ptrBlock;
        const Real h3   = info.h*info.h*info.h;
        const Real facA = -dt/info.h * h3 * coef;
        const Real facD = (mu/info.h)*(dt/info.h) * h3 * coef;

        for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
            const Real uAbs[3] = { lab(ix,iy,iz).u + uInf[0], lab(ix,iy,iz).v + uInf[1], lab(ix,iy,iz).w + uInf[2] };
            const Real dudx = derivative(uAbs[0],lab(ix-3,iy,iz).u,lab(ix-2,iy,iz).u,lab(ix-1,iy,iz).u,lab(ix,iy,iz).u,lab(ix+1,iy,iz).u,lab(ix+2,iy,iz).u,lab(ix+3,iy,iz).u);
            const Real dvdx = derivative(uAbs[0],lab(ix-3,iy,iz).v,lab(ix-2,iy,iz).v,lab(ix-1,iy,iz).v,lab(ix,iy,iz).v,lab(ix+1,iy,iz).v,lab(ix+2,iy,iz).v,lab(ix+3,iy,iz).v);
            const Real dwdx = derivative(uAbs[0],lab(ix-3,iy,iz).w,lab(ix-2,iy,iz).w,lab(ix-1,iy,iz).w,lab(ix,iy,iz).w,lab(ix+1,iy,iz).w,lab(ix+2,iy,iz).w,lab(ix+3,iy,iz).w);
            const Real dudy = derivative(uAbs[1],lab(ix,iy-3,iz).u,lab(ix,iy-2,iz).u,lab(ix,iy-1,iz).u,lab(ix,iy,iz).u,lab(ix,iy+1,iz).u,lab(ix,iy+2,iz).u,lab(ix,iy+3,iz).u);
            const Real dvdy = derivative(uAbs[1],lab(ix,iy-3,iz).v,lab(ix,iy-2,iz).v,lab(ix,iy-1,iz).v,lab(ix,iy,iz).v,lab(ix,iy+1,iz).v,lab(ix,iy+2,iz).v,lab(ix,iy+3,iz).v);
            const Real dwdy = derivative(uAbs[1],lab(ix,iy-3,iz).w,lab(ix,iy-2,iz).w,lab(ix,iy-1,iz).w,lab(ix,iy,iz).w,lab(ix,iy+1,iz).w,lab(ix,iy+2,iz).w,lab(ix,iy+3,iz).w);
            const Real dudz = derivative(uAbs[2],lab(ix,iy,iz-3).u,lab(ix,iy,iz-2).u,lab(ix,iy,iz-1).u,lab(ix,iy,iz).u,lab(ix,iy,iz+1).u,lab(ix,iy,iz+2).u,lab(ix,iy,iz+3).u);
            const Real dvdz = derivative(uAbs[2],lab(ix,iy,iz-3).v,lab(ix,iy,iz-2).v,lab(ix,iy,iz-1).v,lab(ix,iy,iz).v,lab(ix,iy,iz+1).v,lab(ix,iy,iz+2).v,lab(ix,iy,iz+3).v);
            const Real dwdz = derivative(uAbs[2],lab(ix,iy,iz-3).w,lab(ix,iy,iz-2).w,lab(ix,iy,iz-1).w,lab(ix,iy,iz).w,lab(ix,iy,iz+1).w,lab(ix,iy,iz+2).w,lab(ix,iy,iz+3).w);
            const Real duD =  lab(ix+1,iy,iz).u + lab(ix-1,iy,iz).u
                            + lab(ix,iy+1,iz).u + lab(ix,iy-1,iz).u
                            + lab(ix,iy,iz+1).u + lab(ix,iy,iz-1).u - 6 * lab(ix,iy,iz).u;
            const Real dvD =  lab(ix+1,iy,iz).v + lab(ix-1,iy,iz).v
                            + lab(ix,iy+1,iz).v + lab(ix,iy-1,iz).v
                            + lab(ix,iy,iz+1).v + lab(ix,iy,iz-1).v - 6 * lab(ix,iy,iz).v;
            const Real dwD =  lab(ix+1,iy,iz).w + lab(ix-1,iy,iz).w
                            + lab(ix,iy+1,iz).w + lab(ix,iy-1,iz).w
                            + lab(ix,iy,iz+1).w + lab(ix,iy,iz-1).w - 6 * lab(ix,iy,iz).w;
            const Real duA = uAbs[0] * dudx + uAbs[1] * dudy + uAbs[2] * dudz;
            const Real dvA = uAbs[0] * dvdx + uAbs[1] * dvdy + uAbs[2] * dvdz;
            const Real dwA = uAbs[0] * dwdx + uAbs[1] * dwdy + uAbs[2] * dwdz;
            o(ix,iy,iz).tmpU =  facA*duA + facD*duD ;
            o(ix,iy,iz).tmpV =  facA*dvA + facD*dvD ;
            o(ix,iy,iz).tmpW =  facA*dwA + facD*dwD ;
        }

        BlockCase<FluidBlock> * tempCase = (BlockCase<FluidBlock> *)(info.auxiliary);
        typename FluidBlock::ElementType * faceXm = nullptr;
        typename FluidBlock::ElementType * faceXp = nullptr;
        typename FluidBlock::ElementType * faceYm = nullptr;
        typename FluidBlock::ElementType * faceYp = nullptr;
        typename FluidBlock::ElementType * faceZp = nullptr;
        typename FluidBlock::ElementType * faceZm = nullptr;
        if (tempCase != nullptr)
        {
          faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
          faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
          faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
          faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
          faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
          faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
        }
        if (faceXm != nullptr)
        {
          const int ix = 0;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          {
            faceXm[iy + FluidBlock::sizeY * iz].clear();
            faceXm[iy + FluidBlock::sizeY * iz].tmpU = facD*(lab(ix,iy,iz).u - lab(ix-1,iy,iz).u);
            faceXm[iy + FluidBlock::sizeY * iz].tmpV = facD*(lab(ix,iy,iz).v - lab(ix-1,iy,iz).v);
            faceXm[iy + FluidBlock::sizeY * iz].tmpW = facD*(lab(ix,iy,iz).w - lab(ix-1,iy,iz).w);
          }
        }
        if (faceXp != nullptr)
        {
          const int ix = FluidBlock::sizeX-1;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          {
            faceXp[iy + FluidBlock::sizeY * iz].clear();
            faceXp[iy + FluidBlock::sizeY * iz].tmpU = facD*(lab(ix,iy,iz).u - lab(ix+1,iy,iz).u);
            faceXp[iy + FluidBlock::sizeY * iz].tmpV = facD*(lab(ix,iy,iz).v - lab(ix+1,iy,iz).v);
            faceXp[iy + FluidBlock::sizeY * iz].tmpW = facD*(lab(ix,iy,iz).w - lab(ix+1,iy,iz).w);
          }
        }
        if (faceYm != nullptr)
        {
          const int iy = 0;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceYm[ix + FluidBlock::sizeX * iz].clear();
            faceYm[ix + FluidBlock::sizeX * iz].tmpU = facD*(lab(ix,iy,iz).u - lab(ix,iy-1,iz).u);
            faceYm[ix + FluidBlock::sizeX * iz].tmpV = facD*(lab(ix,iy,iz).v - lab(ix,iy-1,iz).v);
            faceYm[ix + FluidBlock::sizeX * iz].tmpW = facD*(lab(ix,iy,iz).w - lab(ix,iy-1,iz).w);
          }
        }
        if (faceYp != nullptr)
        {
          const int iy = FluidBlock::sizeY-1;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceYp[ix + FluidBlock::sizeX * iz].clear();
            faceYp[ix + FluidBlock::sizeX * iz].tmpU = facD*(lab(ix,iy,iz).u - lab(ix,iy+1,iz).u);
            faceYp[ix + FluidBlock::sizeX * iz].tmpV = facD*(lab(ix,iy,iz).v - lab(ix,iy+1,iz).v);
            faceYp[ix + FluidBlock::sizeX * iz].tmpW = facD*(lab(ix,iy,iz).w - lab(ix,iy+1,iz).w);
          }
        }
        if (faceZm != nullptr)
        {
          const int iz = 0;
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceZm[ix + FluidBlock::sizeX * iy].clear();
            faceZm[ix + FluidBlock::sizeX * iy].tmpU = facD*(lab(ix,iy,iz).u - lab(ix,iy,iz-1).u);
            faceZm[ix + FluidBlock::sizeX * iy].tmpV = facD*(lab(ix,iy,iz).v - lab(ix,iy,iz-1).v);
            faceZm[ix + FluidBlock::sizeX * iy].tmpW = facD*(lab(ix,iy,iz).w - lab(ix,iy,iz-1).w);
          }
        }
        if (faceZp != nullptr)
        {
          const int iz = FluidBlock::sizeZ-1;
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceZp[ix + FluidBlock::sizeX * iy].clear();
            faceZp[ix + FluidBlock::sizeX * iy].tmpU = facD*(lab(ix,iy,iz).u - lab(ix,iy,iz+1).u);
            faceZp[ix + FluidBlock::sizeX * iy].tmpV = facD*(lab(ix,iy,iz).v - lab(ix,iy,iz+1).v);
            faceZp[ix + FluidBlock::sizeX * iy].tmpW = facD*(lab(ix,iy,iz).w - lab(ix,iy,iz+1).w);
          }
        }
    }
};

}

void AdvectionDiffusion::operator()(const Real dt)
{
    //Midpoint integration

    const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();

    //1.Save u^{n} to dataOld
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
        FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
        for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
            b.dataOld[iz][iy][ix][0] = b(ix,iy,iz).u;
            b.dataOld[iz][iy][ix][1] = b(ix,iy,iz).v;
            b.dataOld[iz][iy][ix][2] = b(ix,iy,iz).w;
        }
    }

    /********************************************************************/
    // 2. Set u^{n+1/2} = u^{n} + 0.5*dt*RHS(u^{n})

    //   2a) Compute 0.5*dt*RHS(u^{n}) and store it to tmpU,tmpV,tmpW
    const KernelAdvectDiffuse step1(sim,0.5);
    cubism::compute<LabMPI>(step1,sim.grid,sim.grid);

    //   2b) Set u^{n+1/2} = u^{n} + 0.5*dt*RHS(u^{n})
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
        FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
        const Real ih3 = 1.0/(vInfo[i].h*vInfo[i].h*vInfo[i].h);
        for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
            b(ix,iy,iz).u = b.dataOld[iz][iy][ix][0] + b(ix,iy,iz).tmpU*ih3;
            b(ix,iy,iz).v = b.dataOld[iz][iy][ix][1] + b(ix,iy,iz).tmpV*ih3;
            b(ix,iy,iz).w = b.dataOld[iz][iy][ix][2] + b(ix,iy,iz).tmpW*ih3;
        }
    }
    /********************************************************************/


    /********************************************************************/
    // 3. Set u^{n+1} = u^{n} + dt*RHS(u^{n+1/2})
    //   3a) Compute dt*RHS(u^{n+1/2}) and store it to tmpU,tmpV,tmpW
    const KernelAdvectDiffuse step2(sim,1.0);
    cubism::compute<LabMPI>(step2,sim.grid,sim.grid);

    //   3b) Set u^{n+1} = u^{n} + dt*RHS(u^{n+1/2})
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
        FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
        const Real ih3 = 1.0/(vInfo[i].h*vInfo[i].h*vInfo[i].h);
        for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
            b(ix,iy,iz).u = b.dataOld[iz][iy][ix][0] + b(ix,iy,iz).tmpU*ih3;
            b(ix,iy,iz).v = b.dataOld[iz][iy][ix][1] + b(ix,iy,iz).tmpV*ih3;
            b(ix,iy,iz).w = b.dataOld[iz][iy][ix][2] + b(ix,iy,iz).tmpW*ih3;
        }
    }
    /********************************************************************/

}

CubismUP_3D_NAMESPACE_END
