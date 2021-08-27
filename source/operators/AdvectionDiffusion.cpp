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
    KernelAdvectDiffuse(const SimulationData&s) : Discretization(s), sim(s) {}
    const SimulationData & sim;
    const Real dt = sim.dt;
    const Real mu = sim.nu;
    const std::array<Real, 3>& uInf = sim.uinf;
    const int loopBeg = this->getStencilBeg();
    const int loopEndX = CUP_BLOCK_SIZEX-1 + this->getStencilEnd();
    const int loopEndY = CUP_BLOCK_SIZEY-1 + this->getStencilEnd();
    const int loopEndZ = CUP_BLOCK_SIZEZ-1 + this->getStencilEnd();
    const Real norUinf = 1 / std::max({std::fabs(uInf[0]), std::fabs(uInf[1]), std::fabs(uInf[2]), EPS});
    const StencilInfo stencil{this->getStencilBeg(), this->getStencilBeg(),
                              this->getStencilBeg(), this->getStencilEnd(),
                              this->getStencilEnd(), this->getStencilEnd(), false, {FE_U, FE_V, FE_W, FE_P}};

    void operator()(LabMPI & lab, const BlockInfo& info, FluidBlock& o) const
    {
        const bool Euler = (sim.TimeOrder == 1 || sim.step < sim.step_2nd_start);
        const Real h3   = info.h*info.h*info.h;
        const Real ih   = 0.5/info.h;
        const Real facA = this->template advectionCoef(dt, info.h);
        const Real facD = this->template diffusionCoef(dt, info.h, mu);

        if (Euler)
        {
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
                o.data[iz][iy][ix].tmpU = h3*( facA*duA + facD*duD );
                o.data[iz][iy][ix].tmpV = h3*( facA*dvA + facD*dvD );
                o.data[iz][iy][ix].tmpW = h3*( facA*dwA + facD*dwD );
            }
        }
        else
        {
            const double aux = h3/sim.coefU[0];
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
                o.data[iz][iy][ix].tmpU =  aux * ( facA*duA + facD*duD - dt*(lab(ix+1,iy,iz).p-lab(ix-1,iy,iz).p)*ih );
                o.data[iz][iy][ix].tmpV =  aux * ( facA*dvA + facD*dvD - dt*(lab(ix,iy+1,iz).p-lab(ix,iy-1,iz).p)*ih );
                o.data[iz][iy][ix].tmpW =  aux * ( facA*dwA + facD*dwD - dt*(lab(ix,iy,iz+1).p-lab(ix,iy,iz-1).p)*ih );
            }
        }

        const double aux = h3/sim.coefU[0];
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
        if (Euler)
        {
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
        else
        {
            if (faceXm != nullptr)
            {
              int ix = 0;
              for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
              for(int iy=0; iy<FluidBlock::sizeY; ++iy)
              {
                faceXm[iy + FluidBlock::sizeY * iz].clear();
                faceXm[iy + FluidBlock::sizeY * iz].tmpU = h3*facD*(lab(ix,iy,iz).u - lab(ix-1,iy,iz).u) - aux*dt*ih*(lab(ix,iy,iz).p+lab(ix-1,iy,iz).p);
              }
            }
            if (faceXp != nullptr)
            {
              int ix = FluidBlock::sizeX-1;
              for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
              for(int iy=0; iy<FluidBlock::sizeY; ++iy)
              {
                faceXp[iy + FluidBlock::sizeY * iz].clear();
                faceXp[iy + FluidBlock::sizeY * iz].tmpU = h3*facD*(lab(ix,iy,iz).u - lab(ix+1,iy,iz).u) + aux*dt*ih*(lab(ix,iy,iz).p+lab(ix+1,iy,iz).p);
              }
            }
            if (faceYm != nullptr)
            {
              int iy = 0;
              for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
              for(int ix=0; ix<FluidBlock::sizeX; ++ix)
              {
                faceYm[ix + FluidBlock::sizeX * iz].clear();
                faceYm[ix + FluidBlock::sizeX * iz].tmpV = h3*facD*(lab(ix,iy,iz).v - lab(ix,iy-1,iz).v) - aux*dt*ih*(lab(ix,iy,iz).p+lab(ix,iy-1,iz).p);
              }
            }
            if (faceYp != nullptr)
            {
              int iy = FluidBlock::sizeY-1;
              for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
              for(int ix=0; ix<FluidBlock::sizeX; ++ix)
              {
                faceYp[ix + FluidBlock::sizeX * iz].clear();
                faceYp[ix + FluidBlock::sizeX * iz].tmpV = h3*facD*(lab(ix,iy,iz).v - lab(ix,iy+1,iz).v) + aux*dt*ih*(lab(ix,iy,iz).p+lab(ix,iy+1,iz).p);
              }
            }
            if (faceZm != nullptr)
            {
              int iz = 0;
              for(int iy=0; iy<FluidBlock::sizeY; ++iy)
              for(int ix=0; ix<FluidBlock::sizeX; ++ix)
              {
                faceZm[ix + FluidBlock::sizeX * iy].clear();
                faceZm[ix + FluidBlock::sizeX * iy].tmpW = h3*facD*(lab(ix,iy,iz).w - lab(ix,iy,iz-1).w) - aux*dt*ih*(lab(ix,iy,iz).p+lab(ix,iy,iz-1).p);
              }
            }
            if (faceZp != nullptr)
            {
              int iz = FluidBlock::sizeZ-1;
              for(int iy=0; iy<FluidBlock::sizeY; ++iy)
              for(int ix=0; ix<FluidBlock::sizeX; ++ix)
              {
                faceZp[ix + FluidBlock::sizeX * iy].clear();
                faceZp[ix + FluidBlock::sizeX * iy].tmpW = h3*facD*(lab(ix,iy,iz).w - lab(ix,iy,iz+1).w) + aux*dt*ih*(lab(ix,iy,iz).p+lab(ix,iy,iz+1).p);
              }
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


struct UpdateAndCorrectInflow
{
    SimulationData & sim;
    FluidGridMPI * const grid = sim.grid;
    const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();

    inline bool isW(const BlockInfo&I) const
    {
        if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return false;
        return I.index[0] == 0;
    };
    inline bool isE(const BlockInfo&I) const
    {
        if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return false;
        return I.index[0] == sim.bpdx*(1 << I.level)-1;
    };
    inline bool isS(const BlockInfo&I) const
    {
        if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return false;
        return I.index[1] == 0;
    };
    inline bool isN(const BlockInfo&I) const
    {
        if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return false;
        return I.index[1] == sim.bpdy*(1 << I.level)-1;
    };
    inline bool isF(const BlockInfo&I) const
    {
        if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return false;
        return I.index[2] == 0;
    };
    inline bool isB(const BlockInfo&I) const
    {
        if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return false;
        return I.index[2] == sim.bpdz*(1 << I.level)-1;
    };

    UpdateAndCorrectInflow(SimulationData & s) : sim(s) { }

    void operate() const
    {
        // no need for this correction if appropriate U,p BCs are applied
        //double sumInflow = 0;
        //#pragma omp parallel for schedule(static) reduction(+:sumInflow)
        #pragma omp parallel for
        for(size_t i=0; i<vInfo.size(); i++)
        {
            FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
            const double ih3 = 1.0/(vInfo[i].h*vInfo[i].h*vInfo[i].h);
            if (sim.TimeOrder == 2)
            {
                const double aux = 1.0/sim.coefU[0];
                if (sim.step >= sim.step_2nd_start)
                {
                    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                    {
                        b(ix,iy,iz).tmpU = b(ix,iy,iz).tmpU*ih3 + aux * ( -sim.coefU[1]*b.data[iz][iy][ix].u-sim.coefU[2]*b.dataOld[iz][iy][ix][0] );
                        b(ix,iy,iz).tmpV = b(ix,iy,iz).tmpV*ih3 + aux * ( -sim.coefU[1]*b.data[iz][iy][ix].v-sim.coefU[2]*b.dataOld[iz][iy][ix][1] );
                        b(ix,iy,iz).tmpW = b(ix,iy,iz).tmpW*ih3 + aux * ( -sim.coefU[1]*b.data[iz][iy][ix].w-sim.coefU[2]*b.dataOld[iz][iy][ix][2] );
                        b.dataOld[iz][iy][ix][0] = b(ix,iy,iz).u;
                        b.dataOld[iz][iy][ix][1] = b(ix,iy,iz).v;
                        b.dataOld[iz][iy][ix][2] = b(ix,iy,iz).w;
                    }

                }
                else
                {
                    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                    {
                        b(ix,iy,iz).tmpU = b(ix,iy,iz).u + b(ix,iy,iz).tmpU*ih3;
                        b(ix,iy,iz).tmpV = b(ix,iy,iz).v + b(ix,iy,iz).tmpV*ih3;
                        b(ix,iy,iz).tmpW = b(ix,iy,iz).w + b(ix,iy,iz).tmpW*ih3;
                        b.dataOld[iz][iy][ix][0] = b(ix,iy,iz).u;
                        b.dataOld[iz][iy][ix][1] = b(ix,iy,iz).v;
                        b.dataOld[iz][iy][ix][2] = b(ix,iy,iz).w;
                    }
                }
            }
            else
            {
                for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    b(ix,iy,iz).tmpU = b(ix,iy,iz).u + b(ix,iy,iz).tmpU*ih3;
                    b(ix,iy,iz).tmpV = b(ix,iy,iz).v + b(ix,iy,iz).tmpV*ih3;
                    b(ix,iy,iz).tmpW = b(ix,iy,iz).w + b(ix,iy,iz).tmpW*ih3;
                }               
            }
            for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
            for (int ix=0; ix<FluidBlock::sizeX; ++ix)
            {
                b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
                b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
                b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
            }

            // no need for this correction if appropriate U,p BCs are applied
            /*
            const Real h2 = vInfo[i].h_gridpoint * vInfo[i].h_gridpoint;
            if(isW(vInfo[i]))
                for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                    sumInflow -= h2*b(0,iy,iz).u;
            if(isE(vInfo[i]))
                for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                    sumInflow += h2*b(FluidBlock::sizeX-1,iy,iz).u;
            if(isS(vInfo[i]))
                for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                    sumInflow -= h2*b(ix,0,iz).v;
            if(isN(vInfo[i]))
                for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                    sumInflow += h2*b(ix,FluidBlock::sizeY-1,iz).v;
            if(isF(vInfo[i]))
                for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                    sumInflow -= h2*b(ix,iy,0).w;
            if(isB(vInfo[i]))
                for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                    sumInflow += h2*b(ix,iy,FluidBlock::sizeZ-1).w;
            */
        }

        // no need for this correction if appropriate U,p BCs are applied
        /*
        MPI_Allreduce(MPI_IN_PLACE, &sumInflow, 1, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

        const double nTotX = FluidBlock::sizeX * sim.bpdx * (1<<(sim.levelMax-1))*sim.hmin;
        const double nTotY = FluidBlock::sizeY * sim.bpdy * (1<<(sim.levelMax-1))*sim.hmin;
        const double nTotZ = FluidBlock::sizeZ * sim.bpdz * (1<<(sim.levelMax-1))*sim.hmin;
        const Real corr = sumInflow / (2*(nTotX*nTotY + nTotX*nTotZ + nTotY*nTotZ));

        if(std::fabs(corr) < EPS) return;
        if(sim.verbose) printf("Inflow correction %e\n", corr);

        #pragma omp parallel for schedule(static)
        for(size_t i=0; i<vInfo.size(); i++)
        {
          FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
          if(isW(vInfo[i]))
            for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                b(0,iy,iz).u += corr;
          if(isE(vInfo[i]))
            for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                b(FluidBlock::sizeX-1,iy,iz).u -= corr;
          if(isS(vInfo[i]))
            for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                b(ix,0,iz).v += corr;
          if(isN(vInfo[i]))
            for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                b(ix,FluidBlock::sizeY-1,iz).v -= corr;
          if(isF(vInfo[i]))
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
            for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                b(ix,iy,0).w += corr;
          if(isB(vInfo[i]))
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
            for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                b(ix,iy,FluidBlock::sizeZ-1).w -= corr;
        }
        */
    }
};

}

void AdvectionDiffusion::operator()(const double dt)
{
    sim.startProfiler("AdvDiff Kernel");
    const KernelAdvectDiffuse<Upwind3rd> K(sim);
    compute(K,true);
    const UpdateAndCorrectInflow U(sim);
    U.operate();
    sim.stopProfiler();
    check("AdvectionDiffusion");
}

CubismUP_3D_NAMESPACE_END
