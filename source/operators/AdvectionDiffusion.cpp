//
//  Cubism3D
//  Copyright (c) 2021 CSE-Lab, ETH Zurich, Switzerland.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "AdvectionDiffusion.h"
#include "../obstacles/ObstacleVector.h"
#include "../poisson/PoissonSolverAMR.h"
CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

namespace {

// input : field with ghosts from which finite differences are computed
template<int i> Real& inp(LabMPI& L, const int ix, const int iy, const int iz);

template<typename Discretization>
struct KernelAdvectDiffuse : public Discretization
{
    KernelAdvectDiffuse(const SimulationData&s, const double a_coef) : Discretization(s), sim(s), coef(a_coef) {}
    const SimulationData & sim;
    const Real dt = sim.dt;
    const Real mu = sim.nu;
    const double coef;
    const std::array<Real, 3>& uInf = sim.uinf;
    const int loopBeg = this->getStencilBeg();
    const int loopEndX = CUP_BLOCK_SIZEX-1 + this->getStencilEnd();
    const int loopEndY = CUP_BLOCK_SIZEY-1 + this->getStencilEnd();
    const int loopEndZ = CUP_BLOCK_SIZEZ-1 + this->getStencilEnd();
    const Real norUinf = 1 / std::max({std::fabs(uInf[0]), std::fabs(uInf[1]), std::fabs(uInf[2]), EPS});
    const StencilInfo stencil{this->getStencilBeg(), this->getStencilBeg(),
                              this->getStencilBeg(), this->getStencilEnd(),
                              this->getStencilEnd(), this->getStencilEnd(), false, {FE_U, FE_V, FE_W}};

    void operator()(LabMPI & lab, const BlockInfo& info, FluidBlock& o) const
    {
        const Real h3   = info.h*info.h*info.h;
        const Real facA = this->template advectionCoef(dt, info.h);
        const Real facD = this->template diffusionCoef(dt, info.h, mu);
        for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
            const Real uAbs[3] = { lab(ix,iy,iz).u + uInf[0], lab(ix,iy,iz).v + uInf[1], lab(ix,iy,iz).w + uInf[2] };
            const Real dudx = this->template diffx<0>(lab, o, uAbs, ix, iy, iz);
            const Real dvdx = this->template diffx<1>(lab, o, uAbs, ix, iy, iz);
            const Real dwdx = this->template diffx<2>(lab, o, uAbs, ix, iy, iz);
            const Real dudy = this->template diffy<0>(lab, o, uAbs, ix, iy, iz);
            const Real dvdy = this->template diffy<1>(lab, o, uAbs, ix, iy, iz);
            const Real dwdy = this->template diffy<2>(lab, o, uAbs, ix, iy, iz);
            const Real dudz = this->template diffz<0>(lab, o, uAbs, ix, iy, iz);
            const Real dvdz = this->template diffz<1>(lab, o, uAbs, ix, iy, iz);
            const Real dwdz = this->template diffz<2>(lab, o, uAbs, ix, iy, iz);
            const Real duD = this->template lap<0>(lab, o, ix, iy, iz);
            const Real dvD = this->template lap<1>(lab, o, ix, iy, iz);
            const Real dwD = this->template lap<2>(lab, o, ix, iy, iz);
            const Real duA = uAbs[0] * dudx + uAbs[1] * dudy + uAbs[2] * dudz;
            const Real dvA = uAbs[0] * dvdx + uAbs[1] * dvdy + uAbs[2] * dvdz;
            const Real dwA = uAbs[0] * dwdx + uAbs[1] * dwdy + uAbs[2] * dwdz;
            o(ix,iy,iz).tmpU = coef * h3*( facA*duA + facD*duD );
            o(ix,iy,iz).tmpV = coef * h3*( facA*dvA + facD*dvD );
            o(ix,iy,iz).tmpW = coef * h3*( facA*dwA + facD*dwD );
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
          int ix = 0;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          {
            faceXm[iy + FluidBlock::sizeY * iz].clear();
            faceXm[iy + FluidBlock::sizeY * iz].tmpU = h3*facD*(lab(ix,iy,iz).u - lab(ix-1,iy,iz).u);
          }
        }
        if (faceXp != nullptr)
        {
          int ix = FluidBlock::sizeX-1;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          {
            faceXp[iy + FluidBlock::sizeY * iz].clear();
            faceXp[iy + FluidBlock::sizeY * iz].tmpU = h3*facD*(lab(ix,iy,iz).u - lab(ix+1,iy,iz).u);
          }
        }
        if (faceYm != nullptr)
        {
          int iy = 0;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceYm[ix + FluidBlock::sizeX * iz].clear();
            faceYm[ix + FluidBlock::sizeX * iz].tmpV = h3*facD*(lab(ix,iy,iz).v - lab(ix,iy-1,iz).v);
          }
        }
        if (faceYp != nullptr)
        {
          int iy = FluidBlock::sizeY-1;
          for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceYp[ix + FluidBlock::sizeX * iz].clear();
            faceYp[ix + FluidBlock::sizeX * iz].tmpV = h3*facD*(lab(ix,iy,iz).v - lab(ix,iy+1,iz).v);
          }
        }
        if (faceZm != nullptr)
        {
          int iz = 0;
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceZm[ix + FluidBlock::sizeX * iy].clear();
            faceZm[ix + FluidBlock::sizeX * iy].tmpW = h3*facD*(lab(ix,iy,iz).w - lab(ix,iy,iz-1).w);
          }
        }
        if (faceZp != nullptr)
        {
          int iz = FluidBlock::sizeZ-1;
          for(int iy=0; iy<FluidBlock::sizeY; ++iy)
          for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          {
            faceZp[ix + FluidBlock::sizeX * iy].clear();
            faceZp[ix + FluidBlock::sizeX * iy].tmpW = h3*facD*(lab(ix,iy,iz).w - lab(ix,iy,iz+1).w);
          }
        }
    }
};

template<> inline Real& inp<0>(LabMPI& L, const int ix, const int iy, const int iz) { return L(ix,iy,iz).u; }
template<> inline Real& inp<1>(LabMPI& L, const int ix, const int iy, const int iz) { return L(ix,iy,iz).v; }
template<> inline Real& inp<2>(LabMPI& L, const int ix, const int iy, const int iz) { return L(ix,iy,iz).w; }

struct Upwind3rd
{
    const SimulationData& sim;
    Upwind3rd(const SimulationData& s) : sim(s) {}

    template<int dir>
    inline Real diffx(LabMPI& L, const FluidBlock& o, const Real uAbs[3], const int ix, const int iy, const int iz) const
    {
        const Real ucc = inp<dir>(L,ix,iy,iz);
        const Real um1 = inp<dir>(L,ix-1,iy,iz);
        const Real um2 = inp<dir>(L,ix-2,iy,iz);
        const Real up1 = inp<dir>(L,ix+1,iy,iz);
        const Real up2 = inp<dir>(L,ix+2,iy,iz);
        return uAbs[0]>0? 2*up1 +3*ucc -6*um1 +um2 : -up2 +6*up1 -3*ucc -2*um1;
    }
    template<int dir>
    inline Real diffy(LabMPI& L, const FluidBlock& o, const Real uAbs[3], const int ix, const int iy, const int iz) const
    {
        const Real ucc = inp<dir>(L,ix,iy,iz);
        const Real um1 = inp<dir>(L,ix,iy-1,iz);
        const Real um2 = inp<dir>(L,ix,iy-2,iz);
        const Real up1 = inp<dir>(L,ix,iy+1,iz);
        const Real up2 = inp<dir>(L,ix,iy+2,iz);
        return uAbs[1]>0? 2*up1 +3*ucc -6*um1 +um2 : -up2 +6*up1 -3*ucc -2*um1;
    }
    template<int dir>
    inline Real diffz(LabMPI& L, const FluidBlock& o, const Real uAbs[3], const int ix, const int iy, const int iz) const
    {
        const Real ucc = inp<dir>(L,ix,iy,iz);
        const Real um1 = inp<dir>(L,ix,iy,iz-1);
        const Real um2 = inp<dir>(L,ix,iy,iz-2);
        const Real up1 = inp<dir>(L,ix,iy,iz+1);
        const Real up2 = inp<dir>(L,ix,iy,iz+2);
        return uAbs[2]>0? 2*up1 +3*ucc -6*um1 +um2 : -up2 +6*up1 -3*ucc -2*um1;
    }
    template<int dir>
    inline Real   lap(LabMPI& L, const FluidBlock& o, const int ix, const int iy, const int iz) const
    {
        return  inp<dir>(L,ix+1,iy,iz) + inp<dir>(L,ix-1,iy,iz)
              + inp<dir>(L,ix,iy+1,iz) + inp<dir>(L,ix,iy-1,iz)
              + inp<dir>(L,ix,iy,iz+1) + inp<dir>(L,ix,iy,iz-1) - 6 * inp<dir>(L,ix,iy,iz);
    }
    Real advectionCoef(const Real dt, const Real h) const {return -dt/(6*h);}
    Real diffusionCoef(const Real dt, const Real h, const Real mu) const {return (mu/h) * (dt/h);}
    int getStencilBeg() const { return -2; }
    int getStencilEnd() const { return  3; }
};

}

void AdvectionDiffusion::operator()(const double dt)
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
    const KernelAdvectDiffuse<Upwind3rd> step1(sim,0.5);
    compute(step1,true);

    //   2b) Set u^{n+1/2} = u^{n} + 0.5*dt*RHS(u^{n})
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
        FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
        const double ih3 = 1.0/(vInfo[i].h*vInfo[i].h*vInfo[i].h);
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
    const KernelAdvectDiffuse<Upwind3rd> step2(sim,1.0);
    compute(step2,true);

    //   3b) Set u^{n+1} = u^{n} + dt*RHS(u^{n+1/2})
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
        FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
        const double ih3 = 1.0/(vInfo[i].h*vInfo[i].h*vInfo[i].h);
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

    check("AdvectionDiffusion");
}

CubismUP_3D_NAMESPACE_END
