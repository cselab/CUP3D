//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/ObstacleManagement.h"
#include "obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

struct VelocityObstacleVisitor : public ObstacleVisitor
{
  const double dt, lambda;
  int * const nSum;
  double * const uSum;

  VelocityObstacleVisitor(SimulationData& s, int*const nS, double*const uS) :
   dt(s.dt), lambda(s.lambda), nSum(nS), uSum(uS) { }

  void visit(Obstacle* const obstacle)
  {
    const auto &bFixFrameOfRef = obstacle->bFixFrameOfRef;
    obstacle->computeVelocities();
    double povU[3];
    obstacle->getTranslationVelocity(povU);

    if (bFixFrameOfRef[0]) { (nSum[0])++; uSum[0] -= povU[0]; }
    if (bFixFrameOfRef[1]) { (nSum[1])++; uSum[1] -= povU[1]; }
    if (bFixFrameOfRef[2]) { (nSum[2])++; uSum[2] -= povU[2]; }
  }
};

class KernelCharacteristicFunction
{
  using v_v_ob = std::vector<std::vector<ObstacleBlock*>*>;
  const v_v_ob & vec_obstacleBlocks;

  public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 5);

  KernelCharacteristicFunction(const v_v_ob& v) : vec_obstacleBlocks(v) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& b) const
  {
    using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
    static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
    const Real h = info.h_gridpoint, inv2h = .5/h, fac1 = .5*h*h;

    for (size_t obst_id = 0; obst_id<vec_obstacleBlocks.size(); obst_id++)
    {
      const auto& obstacleBlocks = * vec_obstacleBlocks[obst_id];
      ObstacleBlock*const o = obstacleBlocks[info.blockID];
      if(o == nullptr) return;
      CHIMAT & __restrict__ CHI = o->chi;
      CHIMAT & __restrict__ SDF = o->sdf;
      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
      for(int iy=0; iy<FluidBlock::sizeY; iy++)
      for(int ix=0; ix<FluidBlock::sizeX; ix++)
      {
        Real p[3]; info.pos(p, ix,iy,iz);
        if (SDF[iz][iy][ix] > +2*h || SDF[iz][iy][ix] < -2*h)
        {
          const Real H = SDF[iz][iy][ix] > 0 ? 1 : 0;
          b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
          CHI[iz][iy][ix] = H;
          o->CoM_x += p[0]*H;
          o->CoM_y += p[1]*H;
          o->CoM_z += p[2]*H;
          o->mass  += H;
          continue;
        }

        const Real distPx = lab(ix+1,iy,iz).tmpU, distMx = lab(ix-1,iy,iz).tmpU;
        const Real distPy = lab(ix,iy+1,iz).tmpU, distMy = lab(ix,iy-1,iz).tmpU;
        const Real distPz = lab(ix,iy,iz+1).tmpU, distMz = lab(ix,iy,iz-1).tmpU;
        // gradU
        const Real gradUX = inv2h*(distPx - distMx);
        const Real gradUY = inv2h*(distPy - distMy);
        const Real gradUZ = inv2h*(distPz - distMz);
        const Real gradUSq = gradUX*gradUX +gradUY*gradUY +gradUZ*gradUZ + EPS;

        const Real IplusX = distPx < 0 ? 0 : distPx;
        const Real IminuX = distMx < 0 ? 0 : distMx;
        const Real IplusY = distPy < 0 ? 0 : distPy;
        const Real IminuY = distMy < 0 ? 0 : distMy;
        const Real IplusZ = distPz < 0 ? 0 : distPz;
        const Real IminuZ = distMz < 0 ? 0 : distMz;
        const Real HplusX = std::fabs(distPx)<EPS ? 0.5 : (distPx < 0 ? 0 : 1);
        const Real HminuX = std::fabs(distMx)<EPS ? 0.5 : (distMx < 0 ? 0 : 1);
        const Real HplusY = std::fabs(distPy)<EPS ? 0.5 : (distPy < 0 ? 0 : 1);
        const Real HminuY = std::fabs(distMy)<EPS ? 0.5 : (distMy < 0 ? 0 : 1);
        const Real HplusZ = std::fabs(distPz)<EPS ? 0.5 : (distPz < 0 ? 0 : 1);
        const Real HminuZ = std::fabs(distMz)<EPS ? 0.5 : (distMz < 0 ? 0 : 1);

        // gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
        const Real gradIX = inv2h*(IplusX - IminuX);
        const Real gradIY = inv2h*(IplusY - IminuY);
        const Real gradIZ = inv2h*(IplusZ - IminuZ);
        const Real gradHX = (HplusX - HminuX);
        const Real gradHY = (HplusY - HminuY);
        const Real gradHZ = (HplusZ - HminuZ);
        const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;
        const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;
        const Real Delta = fac1 * numD/gradUSq; //h^3 * Delta
        const Real H     = numH/gradUSq;

        CHI[iz][iy][ix] = H;
        if (Delta>1e-6) o->write(ix, iy, iz, Delta, gradUX, gradUY, gradUZ);
        b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
        o->CoM_x += p[0]*H;
        o->CoM_y += p[1]*H;
        o->CoM_z += p[2]*H;
        o->mass  += H;
      }
      o->allocate_surface();
    }
  }
};

void CreateObstacles::operator()(const double dt)
{
  sim.startProfiler("Obst. Reset");
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < vInfo.size(); ++i)
  {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      b(ix,iy,iz).chi = 0; b(ix,iy,iz).tmpU = 0;
    }
  }
  sim.stopProfiler();

  sim.startProfiler("Obst. Create");
  sim.obstacle_vector->create();
  {
    auto vecOB = sim.obstacle_vector->getAllObstacleBlocks();
    const KernelCharacteristicFunction K(vecOB);
    compute<KernelCharacteristicFunction>(K);
  }
  sim.obstacle_vector->finalize();
  sim.stopProfiler();
  check("obst. create");
}

void UpdateObstacles::operator()(const double dt)
{
  int nSum[3] = {0,0,0};
  double uSum[3] = {0,0,0};

  sim.startProfiler("Obst. Update");
  ObstacleVisitor* VIS = new VelocityObstacleVisitor(sim, nSum, uSum);
  sim.obstacle_vector->Accept(VIS); // accept you son of a french cow
  if( nSum[0] > 0 ) sim.uinf[0] = uSum[0] / nSum[0];
  if( nSum[1] > 0 ) sim.uinf[1] = uSum[1] / nSum[1];
  if( nSum[2] > 0 ) sim.uinf[2] = uSum[2] / nSum[2];
  //if(rank == 0) if(nSum[0] || nSum[1] || nSum[2])
  //  printf("New Uinf %g %g %g (from %d %d %d)\n",
  //  uInf[0],uInf[1],uInf[2],nSum[0],nSum[1],nSum[2]);
  delete VIS;

  sim.obstacle_vector->update();

  sim.stopProfiler();
  check("obst. update");
}

CubismUP_3D_NAMESPACE_END
