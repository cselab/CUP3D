//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Sid Verma in May 2018.
//

#include "ComputeDissipation.h"
#include "../Utils/BufferedLogger.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

class KernelDissipation
{
  public:
  const Real dt, nu, center[3];
  Real * QOI;
  StencilInfo stencil{-1,-1,-1, 2,2,2, false, {0,1,2}};
  StencilInfo stencil2{-1,-1,-1, 2,2,2, false, {0}};
  SimulationData & sim;

  const std::vector<cubism::BlockInfo>& chiInfo = sim.chiInfo();

  KernelDissipation(Real _dt, const Real ext[3], Real _nu, Real * RDX, SimulationData& s)
  : dt(_dt), nu(_nu), center{ext[0]/2, ext[1]/2, ext[2]/2}, QOI(RDX), sim(s) { }

  void operator()(VectorLab& lab, ScalarLab& pLab, const BlockInfo& info, const BlockInfo& info2) const
  {
    const Real h = info.h;
    const Real hCube = std::pow(h,3), inv2h = .5 / h, invHh = 1/(h*h);
    const ScalarBlock & chiBlock = *(ScalarBlock*)chiInfo[info.blockID].ptrBlock;

    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix)
    {
      const VectorElement &L =lab(ix,iy,iz);
      const VectorElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const VectorElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const VectorElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real X = chiBlock(ix,iy,iz).s;

      Real p[3]; info.pos(p, ix, iy, iz);
      const Real PX = p[0]-center[0], PY = p[1]-center[1], PZ = p[2]-center[2];
      // vorticity
      const Real WX = inv2h * ( (LN.u[2]-LS.u[2]) - (LB.u[1]-LF.u[1]) );
      const Real WY = inv2h * ( (LB.u[0]-LF.u[0]) - (LE.u[2]-LW.u[2]) );
      const Real WZ = inv2h * ( (LE.u[1]-LW.u[1]) - (LN.u[0]-LS.u[0]) );
      //  - \nabla P \cdot \bm{u}
      const Real dPdx = inv2h * (pLab(ix+1,iy,iz).s - pLab(ix-1,iy,iz).s);
      const Real dPdy = inv2h * (pLab(ix,iy+1,iz).s - pLab(ix,iy-1,iz).s);
      const Real dPdz = inv2h * (pLab(ix,iy,iz+1).s - pLab(ix,iy,iz-1).s);
      //  + \mu \bm{u} \cdot \nabla^2 \bm{u}
      const Real lapU = invHh * ( LE.u[0]+LW.u[0]+LN.u[0]+LS.u[0]+LB.u[0]+LF.u[0] -6*L.u[0] );
      const Real lapV = invHh * ( LE.u[1]+LW.u[1]+LN.u[1]+LS.u[1]+LB.u[1]+LF.u[1] -6*L.u[1] );
      const Real lapW = invHh * ( LE.u[2]+LW.u[2]+LN.u[2]+LS.u[2]+LB.u[2]+LF.u[2] -6*L.u[2] );
      const Real V1 = lapU * L.u[0] + lapV * L.u[1] + lapW * L.u[2];
      // + 2 \mu \bm{D} : \bm{D}
      const Real D11 = inv2h*(LE.u[0] - LW.u[0]); // shear stresses
      const Real D22 = inv2h*(LN.u[1] - LS.u[1]); // shear stresses
      const Real D33 = inv2h*(LB.u[2] - LF.u[2]); // shear stresses
      const Real D12 = inv2h*(LN.u[0] - LS.u[0] + LE.u[1] - LW.u[1])/2; // shear stresses
      const Real D13 = inv2h*(LB.u[0] - LF.u[0] + LE.u[2] - LW.u[2])/2; // shear stresses
      const Real D23 = inv2h*(LN.u[2] - LS.u[2] + LB.u[1] - LF.u[1])/2; // shear stresses
      const Real V2 = D11*D11 +D22*D22 +D33*D33 +2*(D12*D12 +D13*D13 +D23*D23);

      #pragma omp critical
      {
        // three linear invariants (conserved in inviscid and viscous flows)
        // conservation of vorticity: int w dx = 0
        QOI[0] += hCube * WX;
        QOI[1] += hCube * WY;
        QOI[2] += hCube * WZ;
        // conservation of linear impulse: int u dx = 0.5 int (x cross w) dx
        QOI[3] += hCube/2 * ( PY * WZ - PZ * WY );
        QOI[4] += hCube/2 * ( PZ * WX - PX * WZ );
        QOI[5] += hCube/2 * ( PX * WY - PY * WX );
        QOI[6] += hCube * L.u[0];
        QOI[7] += hCube * L.u[1];
        QOI[8] += hCube * L.u[2];
        //conserve ang imp.: int (x cross u)dx = 1/3 int (x cross (x cross w) )dx
        // = 1/3 int x (w \cdot x) - w (x \cdot x) ) dx (some terms cancel)
        QOI[9] += hCube/3 * ( PX*(PY*WY + PZ*WZ) - WX*(PY*PY + PZ*PZ) );
        QOI[10] += hCube/3 * ( PY*(PX*WX + PZ*WZ) - WY*(PX*PX + PZ*PZ) );
        QOI[11] += hCube/3 * ( PZ*(PX*WX + PY*WY) - WZ*(PX*PX + PY*PY) );
        QOI[12] += hCube * ( PY * L.u[2] - PZ * L.u[1] );
        QOI[13] += hCube * ( PZ * L.u[0] - PX * L.u[2] );
        QOI[14] += hCube * ( PX * L.u[1] - PY * L.u[0] );

        //presPow
        //viscPow
        QOI[15] -= (1-X) * hCube * ( dPdx * L.u[0] + dPdy * L.u[1] + dPdz * L.u[2]);
        QOI[16] += (1-X) * hCube*nu * ( V1 + 2*V2 );

        // two quadratic invariants: kinetic energy (from solver) and helicity (conserved in inviscid flows)
        //helicity 
        //kineticEn
        //enstrophy
        QOI[17] += hCube * ( WX * L.u[0] + WY * L.u[1] + WZ * L.u[2] );
        QOI[18] += hCube * ( L.u[0]*L.u[0] + L.u[1]*L.u[1] + L.u[2]*L.u[2] )/2;
        QOI[19] += hCube*std::sqrt(WX * WX + WY * WY + WZ * WZ);
      }
    }
  }
};
}

void ComputeDissipation::operator()(const Real dt)
{
  if(sim.freqDiagnostics == 0 || sim.step % sim.freqDiagnostics) return;

  Real RDX[20] = { 0.0 };
  KernelDissipation diss(dt, sim.extents.data(), sim.nu, RDX,sim);
  cubism::compute<KernelDissipation,VectorGrid,VectorLab,ScalarGrid,ScalarLab>(diss,*sim.vel,*sim.pres);

  MPI_Allreduce(MPI_IN_PLACE, RDX, 20,MPI_Real, MPI_SUM,sim.comm);

  size_t loc = sim.velInfo().size();
  size_t tot;
  MPI_Reduce(&loc, &tot, 1, MPI_LONG, MPI_SUM, 0, sim.comm);
  if(sim.rank==0 && sim.muteAll == false)
  {
    std::ofstream outfile;
    outfile.open("diagnostics.dat", std::ios_base::app);
    if(sim.step==0)
      outfile<<"step_id time circ_x circ_y circ_y linImp_x linImp_y linImp_z "
      "linMom_x linMom_y linMom_z angImp_x angImp_y angImp_z angMom_x angMom_y "
      "angMom_z presPow viscPow helicity kineticEn enstrophy blocks"<<std::endl;
    outfile<<sim.step<<" "<<sim.time<<" "<<
    RDX[ 0]<<" "<<RDX[ 1]<<" "<<RDX[ 2]<<" "<<RDX[ 3]<<" "<<RDX[ 4]<<" "<<
    RDX[ 5]<<" "<<RDX[ 6]<<" "<<RDX[ 7]<<" "<<RDX[ 8]<<" "<<RDX[ 9]<<" "<<
    RDX[10]<<" "<<RDX[11]<<" "<<RDX[12]<<" "<<RDX[13]<<" "<<RDX[14]<<" "<<
    RDX[15]<<" "<<RDX[16]<<" "<<RDX[17]<<" "<<RDX[18]<<" "<<RDX[19]<<" "<< tot << std::endl;
    outfile.close();
  }
}

CubismUP_3D_NAMESPACE_END
