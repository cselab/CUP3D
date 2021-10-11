//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Penalization.h"
#include "../obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];

template<bool implicitPenalization>
struct KernelPenalization : public ObstacleVisitor
{
  const Real dt, invdt = 1.0/dt, lambda;
  ObstacleVector * const obstacle_vector;
  const cubism::BlockInfo * info_ptr = nullptr;

  KernelPenalization(double _dt, double _lambda, ObstacleVector* ov) :
    dt(_dt), lambda(_lambda), obstacle_vector(ov) {}

  void operator()(const cubism::BlockInfo& info)
  {
    // first store the lab and info, then do visitor
    info_ptr = & info;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr;
  }

  void visit(Obstacle* const obstacle)
  {
    const BlockInfo& info = * info_ptr;
    assert(info_ptr not_eq nullptr);
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const std::array<double,3> vel = obstacle->getTranslationVelocity();
    const std::array<double,3> omega = obstacle->getAngularVelocity();
    const Real dv = std::pow(info.h_gridpoint, 3);

    // Obstacle-specific lambda, useful for gradually adding an obstacle to the flow.
    const double rampUp = obstacle->lambda_factor;
    // lambda = 1/dt hardcoded for expl time int, other options are wrong.
    const double lambdaFac = rampUp * (implicitPenalization? lambda : invdt);

    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // What if multiple obstacles share a block? Do not write udef onto
      // grid if CHI stored on the grid is greater than obst's CHI.
      if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
      if(CHI[iz][iy][ix] <= 0) continue; // no need to do anything
      double p[3]; info.pos(p, ix, iy, iz);
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      const double U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };
      const Real X = CHI[iz][iy][ix];
      const Real penalFac = implicitPenalization? X*lambdaFac/(1+X*lambdaFac*dt):X*lambdaFac;
      const Real FPX = penalFac * (U_TOT[0] - b(ix,iy,iz).u);
      const Real FPY = penalFac * (U_TOT[1] - b(ix,iy,iz).v);
      const Real FPZ = penalFac * (U_TOT[2] - b(ix,iy,iz).w);
      // What if two obstacles overlap? Let's plus equal. We will need a
      // repulsion term of the velocity at some point in the code.
      b(ix,iy,iz).u = b(ix,iy,iz).u + dt * FPX;
      b(ix,iy,iz).v = b(ix,iy,iz).v + dt * FPY;
      b(ix,iy,iz).w = b(ix,iy,iz).w + dt * FPZ;

      FX += dv * FPX; FY += dv * FPY; FZ += dv * FPZ;
      TX += dv * ( p[1] * FPZ - p[2] * FPY );
      TY += dv * ( p[2] * FPX - p[0] * FPZ );
      TZ += dv * ( p[0] * FPY - p[1] * FPX );
    }
  }
};

struct KernelFinalizePenalizationForce : public ObstacleVisitor
{
  FluidGridMPI * const grid;

  KernelFinalizePenalizationForce(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obst)
  {
    static constexpr int nQoI = 6;
    double M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); ++i) {
      if(oBlock[i] == nullptr) continue;
      M[0] += oBlock[i]->FX; M[1] += oBlock[i]->FY; M[2] += oBlock[i]->FZ;
      M[3] += oBlock[i]->TX; M[4] += oBlock[i]->TY; M[5] += oBlock[i]->TZ;
    }
    const auto comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_DOUBLE, MPI_SUM, comm);
    obst->force[0]  = M[0]; obst->force[1]  = M[1]; obst->force[2]  = M[2];
    obst->torque[0] = M[3]; obst->torque[1] = M[4]; obst->torque[2] = M[5];
  }
};

#if 0
void ElasticCollision(const double m1,
                      const double m2,
                      const double *I1,
                      const double *I2,
                      const double *v1,
                      const double *v2,
                      const double *o1,
                      const double *o2,
                      double *hv1,
                      double *hv2,
                      double *ho1,
                      double *ho2,
                      const double *C1,
                      const double *C2,
                      const double NX,
                      const double NY,
                      const double NZ,
                      const double CX,
                      const double CY,
                      const double CZ)
{
    double N[3] ={NX,NY,NZ};
    double C[3] ={CX,CY,CZ};

    double R1[3];
    double R2[3];
    //R1 = (rc-c1) x n
    R1[0] = ( C[1] - C1[1] )*N[2] - ( C[2] - C1[2] )*N[1];
    R1[1] = ( C[2] - C1[2] )*N[0] - ( C[0] - C1[0] )*N[2];
    R1[2] = ( C[0] - C1[0] )*N[1] - ( C[1] - C1[1] )*N[0];
    //R2 = -(rc-c2) x n
    R2[0] = - ( ( C[1] - C2[1] )*N[2] - ( C[2] - C2[2] )*N[1] );
    R2[1] = - ( ( C[2] - C2[2] )*N[0] - ( C[0] - C2[0] )*N[2] );
    R2[2] = - ( ( C[0] - C2[0] )*N[1] - ( C[1] - C2[1] )*N[0] );

    const Real m00 = I1[0];
    const Real m01 = I1[3];
    const Real m02 = I1[4];
    const Real m11 = I1[1];
    const Real m12 = I1[5];
    const Real m22 = I1[2];
    Real a00 = m22*m11 - m12*m12;
    Real a01 = m02*m12 - m22*m01;
    Real a02 = m01*m12 - m02*m11;
    Real a11 = m22*m00 - m02*m02;
    Real a12 = m01*m02 - m00*m12;
    Real a22 = m00*m11 - m01*m01;
    const Real determinant =  1.0/((m00 * a00) + (m01 * a01) + (m02 * a02));
    a00 *= determinant;
    a01 *= determinant;
    a02 *= determinant;
    a11 *= determinant;
    a12 *= determinant;
    a22 *= determinant;

    const Real n00 = I2[0];
    const Real n01 = I2[3];
    const Real n02 = I2[4];
    const Real n11 = I2[1];
    const Real n12 = I2[5];
    const Real n22 = I2[2];
    Real b00 = n22*n11 - n12*n12;
    Real b01 = n02*n12 - n22*n01;
    Real b02 = n01*n12 - n02*n11;
    Real b11 = n22*n00 - n02*n02;
    Real b12 = n01*n02 - n00*n12;
    Real b22 = n00*n11 - n01*n01;
    const Real determinant2 =  1.0/((n00 * b00) + (n01 * b01) + (n02 * b02));
    b00 *= determinant2;
    b01 *= determinant2;
    b02 *= determinant2;
    b11 *= determinant2;
    b12 *= determinant2;
    b22 *= determinant2;

    double k1[3];
    double k2[3];
    double J1[3];
    double J2[3];
    k1[0] =  N[0]/m1;
    k1[1] =  N[1]/m1;
    k1[2] =  N[2]/m1;
    k2[0] = -N[0]/m2;
    k2[1] = -N[1]/m2;
    k2[2] = -N[2]/m2;
    J1[0] = a00*R1[0] + a01*R1[1] + a02*R1[2];
    J1[1] = a01*R1[0] + a11*R1[1] + a12*R1[2];
    J1[2] = a02*R1[0] + a12*R1[1] + a22*R1[2];
    J2[0] = b00*R2[0] + b01*R2[1] + b02*R2[2];
    J2[1] = b01*R2[0] + b11*R2[1] + b12*R2[2];
    J2[2] = b02*R2[0] + b12*R2[1] + b22*R2[2];

    double nom   = 0;
    double denom = 0;
    nom = 2*m1*( v1[0]*k1[0] + v1[1]*k1[1] + v1[2]*k1[2])
        + 2*m2*( v2[0]*k2[0] + v2[1]*k2[1] + v2[2]*k2[2]);
    denom = m1*(k1[0]*k1[0]+k1[1]*k1[1]+k1[2]*k1[2])+
            m2*(k2[0]*k2[0]+k2[1]*k2[1]+k2[2]*k2[2]);
    double II1[9] = {m00,m01,m02,m01,m11,m12,m02,m12,m22};
    double II2[9] = {n00,n01,n02,n01,n11,n12,n02,n12,n22};
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    {
      nom += II1[3*i+j]*(o1[i]*J1[j]+o1[j]*J1[i]);
      nom += II2[3*i+j]*(o2[i]*J2[j]+o2[j]*J2[i]);
      denom += II1[3*i+j]*J1[i]*J1[j];
      denom += II2[3*i+j]*J2[i]*J2[j];
    }
    const double impulse = -nom/denom;

    hv1[0] = v1[0] + k1[0]*impulse;
    hv1[1] = v1[1] + k1[1]*impulse;
    hv1[2] = v1[2] + k1[2]*impulse;
    hv2[0] = v2[0] + k2[0]*impulse;
    hv2[1] = v2[1] + k2[1]*impulse;
    hv2[2] = v2[2] + k2[2]*impulse;
    ho1[0] = o1[0] + J1[0]*impulse;
    ho1[1] = o1[1] + J1[1]*impulse;
    ho1[2] = o1[2] + J1[2]*impulse;
    ho2[0] = o2[0] + J2[0]*impulse;
    ho2[1] = o2[1] + J2[1]*impulse;
    ho2[2] = o2[2] + J2[2]*impulse;

    double vp1_x = v1[0] + o1[1]*(C[2]-C1[2]) - o1[2]*(C[1]-C1[1]);
    double vp1_y = v1[1] + o1[2]*(C[0]-C1[0]) - o1[0]*(C[2]-C1[2]);
    double vp1_z = v1[2] + o1[0]*(C[1]-C1[1]) - o1[1]*(C[0]-C1[0]);
    double hvp1_x = hv1[0] + ho1[1]*(C[2]-C1[2]) - ho1[2]*(C[1]-C1[1]);
    double hvp1_y = hv1[1] + ho1[2]*(C[0]-C1[0]) - ho1[0]*(C[2]-C1[2]);
    double hvp1_z = hv1[2] + ho1[0]*(C[1]-C1[1]) - ho1[1]*(C[0]-C1[0]);


    double vp2_x = v2[0] + o2[1]*(C[2]-C2[2]) - o2[2]*(C[1]-C2[1]);
    double vp2_y = v2[1] + o2[2]*(C[0]-C2[0]) - o2[0]*(C[2]-C2[2]);
    double vp2_z = v2[2] + o2[0]*(C[1]-C2[1]) - o2[1]*(C[0]-C2[0]);
    double hvp2_x = hv2[0] + ho2[1]*(C[2]-C2[2]) - ho2[2]*(C[1]-C2[1]);
    double hvp2_y = hv2[1] + ho2[2]*(C[0]-C2[0]) - ho2[0]*(C[2]-C2[2]);
    double hvp2_z = hv2[2] + ho2[0]*(C[1]-C2[1]) - ho2[1]*(C[0]-C2[0]);


    double vv   =( vp2_x- vp1_x)*NX + ( vp2_y- vp1_y)*NY + ( vp2_z- vp1_z)*NZ;
    double hvv = (hvp2_x-hvp1_x)*NX + (hvp2_y-hvp1_y)*NY + (hvp2_z-hvp1_z)*NZ;

    #if 1
    int rank ; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (std::fabs(denom) < 1e-21) std::cout << "DENOM=" << denom << " nom="<< nom << std::endl;
    if (rank == 0) std::cout << "            impulse=" << impulse << " " << impulse*NX/m1 << " " << impulse*NY/m1 << " " << impulse*NZ/m1 << std::endl;
    //Total energy
    const double E = m1*  (v1[0]*v1[0] +       v1[1]*v1[1] +       v1[2]*v1[2]) +
                     I1[0]*o1[0]*o1[0] + I1[3]*o1[0]*o1[1] + I1[4]*o1[0]*o1[2]+
                     I1[3]*o1[1]*o1[0] + I1[1]*o1[1]*o1[1] + I1[5]*o1[1]*o1[2]+
                     I1[4]*o1[2]*o1[0] + I1[5]*o1[2]*o1[1] + I1[2]*o1[2]*o1[2]+
                     m2*  (v2[0]*v2[0] +       v2[1]*v2[1] +       v2[2]*v2[2]) +
                     I2[0]*o2[0]*o2[0] + I2[3]*o2[0]*o2[1] + I2[4]*o2[0]*o2[2]+
                     I2[3]*o2[1]*o2[0] + I2[1]*o2[1]*o2[1] + I2[5]*o2[1]*o2[2]+
                     I2[4]*o2[2]*o2[0] + I2[5]*o2[2]*o2[1] + I2[2]*o2[2]*o2[2];
    const double hE= m1*  (hv1[0]*hv1[0] +       hv1[1]*hv1[1] +       hv1[2]*hv1[2]) +
                     I1[0]*ho1[0]*ho1[0] + I1[3]*ho1[0]*ho1[1] + I1[4]*ho1[0]*ho1[2]+
                     I1[3]*ho1[1]*ho1[0] + I1[1]*ho1[1]*ho1[1] + I1[5]*ho1[1]*ho1[2]+
                     I1[4]*ho1[2]*ho1[0] + I1[5]*ho1[2]*ho1[1] + I1[2]*ho1[2]*ho1[2]+
                     m2*  (hv2[0]*hv2[0] +       hv2[1]*hv2[1] +       hv2[2]*hv2[2]) +
                     I2[0]*ho2[0]*ho2[0] + I2[3]*ho2[0]*ho2[1] + I2[4]*ho2[0]*ho2[2]+
                     I2[3]*ho2[1]*ho2[0] + I2[1]*ho2[1]*ho2[1] + I2[5]*ho2[1]*ho2[2]+
                     I2[4]*ho2[2]*ho2[0] + I2[5]*ho2[2]*ho2[1] + I2[2]*ho2[2]*ho2[2];
    const double Lx  = m1 * ( -  v1[1]*C1[2] +  v1[2]*C1[1]) + m2 * ( -  v2[1]*C2[2] +  v2[2]*C2[1]) + I1[0] *o1[0] + I1[3] *o1[1] + I1[4] *o1[2] + I2[0] *o2[0] + I2[3] *o2[1] + I2[4] *o2[2];
    const double Ly  = m1 * ( -  v1[2]*C1[0] +  v1[0]*C1[2]) + m2 * ( -  v2[2]*C2[0] +  v2[0]*C2[2]) + I1[3] *o1[0] + I1[1] *o1[1] + I1[5] *o1[2] + I2[3] *o2[0] + I2[1] *o2[1] + I2[5] *o2[2];
    const double Lz  = m1 * ( -  v1[0]*C1[1] +  v1[1]*C1[0]) + m2 * ( -  v2[0]*C2[1] +  v2[1]*C2[0]) + I1[4] *o1[0] + I1[5] *o1[1] + I1[2] *o1[2] + I2[4] *o2[0] + I2[5] *o2[1] + I2[2] *o2[2];
    const double hLx = m1 * ( - hv1[1]*C1[2] + hv1[2]*C1[1]) + m2 * ( - hv2[1]*C2[2] + hv2[2]*C2[1]) + I1[0]*ho1[0] + I1[3]*ho1[1] + I1[4]*ho1[2] + I2[0]*ho2[0] + I2[3]*ho2[1] + I2[4]*ho2[2];
    const double hLy = m1 * ( - hv1[2]*C1[0] + hv1[0]*C1[2]) + m2 * ( - hv2[2]*C2[0] + hv2[0]*C2[2]) + I1[3]*ho1[0] + I1[1]*ho1[1] + I1[5]*ho1[2] + I2[3]*ho2[0] + I2[1]*ho2[1] + I2[5]*ho2[2];
    const double hLz = m1 * ( - hv1[0]*C1[1] + hv1[1]*C1[0]) + m2 * ( - hv2[0]*C2[1] + hv2[1]*C2[0]) + I1[4]*ho1[0] + I1[5]*ho1[1] + I1[2]*ho1[2] + I2[4]*ho2[0] + I2[5]*ho2[1] + I2[2]*ho2[2];
    const double Mx = m1*v1[0] + m2*v2[0];
    const double My = m1*v1[1] + m2*v2[1];
    const double Mz = m1*v1[2] + m2*v2[2];
    const double hMx = m1*hv1[0] + m2*hv2[0];
    const double hMy = m1*hv1[1] + m2*hv2[1];
    const double hMz = m1*hv1[2] + m2*hv2[2];

    if (rank == 0)
    {
        std::cout << "vv = " << vv << std::endl;
        std::cout << "hvv = " << hvv << std::endl;

        std::cout << "E = " << E  << std::endl;
        std::cout << "hE= " << hE << std::endl;
        std::cout << " L= " << Lx << " " << Ly << " " << Lz << std::endl;
        std::cout << "hL= " << hLx << " " << hLy << " " << hLz << std::endl;
        std::cout << " M= " << Mx << " " << My << " " << Mz << std::endl;
        std::cout << "hM= " << hMx << " " << hMy << " " << hMz << std::endl;
    }
    #endif

}
#endif

void ComputeJ(const double * Rc, const double * R, const double * N, const double * I, double *J)
{
    //Invert I
    const Real m00 = I[0];
    const Real m01 = I[3];
    const Real m02 = I[4];
    const Real m11 = I[1];
    const Real m12 = I[5];
    const Real m22 = I[2];
    Real a00 = m22*m11 - m12*m12;
    Real a01 = m02*m12 - m22*m01;
    Real a02 = m01*m12 - m02*m11;
    Real a11 = m22*m00 - m02*m02;
    Real a12 = m01*m02 - m00*m12;
    Real a22 = m00*m11 - m01*m01;
    const Real determinant =  1.0/((m00 * a00) + (m01 * a01) + (m02 * a02));
    a00 *= determinant;
    a01 *= determinant;
    a02 *= determinant;
    a11 *= determinant;
    a12 *= determinant;
    a22 *= determinant;

    const double aux_0 = ( Rc[1] - R[1] )*N[2] - ( Rc[2] - R[2] )*N[1];
    const double aux_1 = ( Rc[2] - R[2] )*N[0] - ( Rc[0] - R[0] )*N[2];
    const double aux_2 = ( Rc[0] - R[0] )*N[1] - ( Rc[1] - R[1] )*N[0];
    J[0] = a00*aux_0 + a01*aux_1 + a02*aux_2;
    J[1] = a01*aux_0 + a11*aux_1 + a12*aux_2;
    J[2] = a02*aux_0 + a12*aux_1 + a22*aux_2;
}


void ElasticCollision1(const double  m1,const double  m2,
                       const double *I1,const double *I2,
                       const double *v1,const double *v2,
                       const double *o1,const double *o2,
                       double *hv1,double *hv2,
                       double *ho1,double *ho2,
                       const double *C1,const double *C2,
                       const double  NX,const double  NY,const double NZ,
                       const double  CX,const double  CY,const double CZ,
                       double *vc1,double *vc2)
{
    const double e = 1.0; // coefficient of restitution
    const double N[3] ={NX,NY,NZ};
    const double C[3] ={CX,CY,CZ};

    const double k1[3] = { N[0]/m1, N[1]/m1, N[2]/m1};
    const double k2[3] = {-N[0]/m2,-N[1]/m2,-N[2]/m2};
    double J1[3];
    double J2[3]; 
    ComputeJ(C,C1,N,I1,J1);
    ComputeJ(C,C2,N,I2,J2);
    J2[0] = -J2[0];
    J2[1] = -J2[1];
    J2[2] = -J2[2];

    double u1DEF[3];
    u1DEF[0] = vc1[0] - v1[0] - ( o1[1]*(C[2]-C1[2]) - o1[2]*(C[1]-C1[1]) );
    u1DEF[1] = vc1[1] - v1[1] - ( o1[2]*(C[0]-C1[0]) - o1[0]*(C[2]-C1[2]) );
    u1DEF[2] = vc1[2] - v1[2] - ( o1[0]*(C[1]-C1[1]) - o1[1]*(C[0]-C1[0]) );
    double u2DEF[3];
    u2DEF[0] = vc2[0] - v2[0] - ( o2[1]*(C[2]-C2[2]) - o2[2]*(C[1]-C2[1]) );
    u2DEF[1] = vc2[1] - v2[1] - ( o2[2]*(C[0]-C2[0]) - o2[0]*(C[2]-C2[2]) );
    u2DEF[2] = vc2[2] - v2[2] - ( o2[0]*(C[1]-C2[1]) - o2[1]*(C[0]-C2[0]) );

    const double nom = e*( (vc1[0]-vc2[0])*N[0] + 
                           (vc1[1]-vc2[1])*N[1] + 
                           (vc1[2]-vc2[2])*N[2] )
                       + ( (v1[0]-v2[0] + u1DEF[0] - u2DEF[0] )*N[0] + 
                           (v1[1]-v2[1] + u1DEF[1] - u2DEF[1] )*N[1] + 
                           (v1[2]-v2[2] + u1DEF[2] - u2DEF[2] )*N[2] )
                  +( (o1[1]*(C[2]-C1[2]) - o1[2]*(C[1]-C1[1]) )* N[0]+
                     (o1[2]*(C[0]-C1[0]) - o1[0]*(C[2]-C1[2]) )* N[1]+
                     (o1[0]*(C[1]-C1[1]) - o1[1]*(C[0]-C1[0]) )* N[2])
                  -( (o2[1]*(C[2]-C2[2]) - o2[2]*(C[1]-C2[1]) )* N[0]+
                     (o2[2]*(C[0]-C2[0]) - o2[0]*(C[2]-C2[2]) )* N[1]+
                     (o2[0]*(C[1]-C2[1]) - o2[1]*(C[0]-C2[0]) )* N[2]);

    const double denom = -(1.0/m1+1.0/m2) + 
               +( ( J1[1]*(C[2]-C1[2]) - J1[2]*(C[1]-C1[1]) ) *(-N[0])+
                  ( J1[2]*(C[0]-C1[0]) - J1[0]*(C[2]-C1[2]) ) *(-N[1])+
                  ( J1[0]*(C[1]-C1[1]) - J1[1]*(C[0]-C1[0]) ) *(-N[2]))
               -( ( J2[1]*(C[2]-C2[2]) - J2[2]*(C[1]-C2[1]) ) *(-N[0])+
                  ( J2[2]*(C[0]-C2[0]) - J2[0]*(C[2]-C2[2]) ) *(-N[1])+
                  ( J2[0]*(C[1]-C2[1]) - J2[1]*(C[0]-C2[0]) ) *(-N[2]));
    const double impulse = nom/(denom+1e-21);
    hv1[0] = v1[0] + k1[0]*impulse;
    hv1[1] = v1[1] + k1[1]*impulse;
    hv1[2] = v1[2] + k1[2]*impulse;
    hv2[0] = v2[0] + k2[0]*impulse;
    hv2[1] = v2[1] + k2[1]*impulse;
    hv2[2] = v2[2] + k2[2]*impulse;
    ho1[0] = o1[0] + J1[0]*impulse;
    ho1[1] = o1[1] + J1[1]*impulse;
    ho1[2] = o1[2] + J1[2]*impulse;
    ho2[0] = o2[0] + J2[0]*impulse;
    ho2[1] = o2[1] + J2[1]*impulse;
    ho2[2] = o2[2] + J2[2]*impulse;

    #if 0
    double vp1_x = v1[0] + o1[1]*(C[2]-C1[2]) - o1[2]*(C[1]-C1[1]);
    double vp1_y = v1[1] + o1[2]*(C[0]-C1[0]) - o1[0]*(C[2]-C1[2]);
    double vp1_z = v1[2] + o1[0]*(C[1]-C1[1]) - o1[1]*(C[0]-C1[0]);
    double hvp1_x = hv1[0] + ho1[1]*(C[2]-C1[2]) - ho1[2]*(C[1]-C1[1]);
    double hvp1_y = hv1[1] + ho1[2]*(C[0]-C1[0]) - ho1[0]*(C[2]-C1[2]);
    double hvp1_z = hv1[2] + ho1[0]*(C[1]-C1[1]) - ho1[1]*(C[0]-C1[0]);

    double vp2_x = v2[0] + o2[1]*(C[2]-C2[2]) - o2[2]*(C[1]-C2[1]);
    double vp2_y = v2[1] + o2[2]*(C[0]-C2[0]) - o2[0]*(C[2]-C2[2]);
    double vp2_z = v2[2] + o2[0]*(C[1]-C2[1]) - o2[1]*(C[0]-C2[0]);
    double hvp2_x = hv2[0] + ho2[1]*(C[2]-C2[2]) - ho2[2]*(C[1]-C2[1]);
    double hvp2_y = hv2[1] + ho2[2]*(C[0]-C2[0]) - ho2[0]*(C[2]-C2[2]);
    double hvp2_z = hv2[2] + ho2[0]*(C[1]-C2[1]) - ho2[1]*(C[0]-C2[0]);

    double vv   =( vp2_x- vp1_x)*NX + ( vp2_y- vp1_y)*NY + ( vp2_z- vp1_z)*NZ;
    double hvv = (hvp2_x-hvp1_x)*NX + (hvp2_y-hvp1_y)*NY + (hvp2_z-hvp1_z)*NZ;

    int rank ; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (std::fabs(denom) < 1e-21) std::cout << "DENOM=" << denom << " nom="<< nom << std::endl;
    if (rank == 0) std::cout << "            impulse=" << impulse << " " << impulse*NX/m1 << " " << impulse*NY/m1 << " " << impulse*NZ/m1 << std::endl;
    //Total energy
    const double E = m1*  (v1[0]*v1[0] +       v1[1]*v1[1] +       v1[2]*v1[2]) +
                     I1[0]*o1[0]*o1[0] + I1[3]*o1[0]*o1[1] + I1[4]*o1[0]*o1[2]+
                     I1[3]*o1[1]*o1[0] + I1[1]*o1[1]*o1[1] + I1[5]*o1[1]*o1[2]+
                     I1[4]*o1[2]*o1[0] + I1[5]*o1[2]*o1[1] + I1[2]*o1[2]*o1[2]+
                     m2*  (v2[0]*v2[0] +       v2[1]*v2[1] +       v2[2]*v2[2]) +
                     I2[0]*o2[0]*o2[0] + I2[3]*o2[0]*o2[1] + I2[4]*o2[0]*o2[2]+
                     I2[3]*o2[1]*o2[0] + I2[1]*o2[1]*o2[1] + I2[5]*o2[1]*o2[2]+
                     I2[4]*o2[2]*o2[0] + I2[5]*o2[2]*o2[1] + I2[2]*o2[2]*o2[2];
    const double hE= m1*  (hv1[0]*hv1[0] +       hv1[1]*hv1[1] +       hv1[2]*hv1[2]) +
                     I1[0]*ho1[0]*ho1[0] + I1[3]*ho1[0]*ho1[1] + I1[4]*ho1[0]*ho1[2]+
                     I1[3]*ho1[1]*ho1[0] + I1[1]*ho1[1]*ho1[1] + I1[5]*ho1[1]*ho1[2]+
                     I1[4]*ho1[2]*ho1[0] + I1[5]*ho1[2]*ho1[1] + I1[2]*ho1[2]*ho1[2]+
                     m2*  (hv2[0]*hv2[0] +       hv2[1]*hv2[1] +       hv2[2]*hv2[2]) +
                     I2[0]*ho2[0]*ho2[0] + I2[3]*ho2[0]*ho2[1] + I2[4]*ho2[0]*ho2[2]+
                     I2[3]*ho2[1]*ho2[0] + I2[1]*ho2[1]*ho2[1] + I2[5]*ho2[1]*ho2[2]+
                     I2[4]*ho2[2]*ho2[0] + I2[5]*ho2[2]*ho2[1] + I2[2]*ho2[2]*ho2[2];
    const double Lx  = m1 * ( -  v1[1]*C1[2] +  v1[2]*C1[1]) + m2 * ( -  v2[1]*C2[2] +  v2[2]*C2[1]) + I1[0] *o1[0] + I1[3] *o1[1] + I1[4] *o1[2] + I2[0] *o2[0] + I2[3] *o2[1] + I2[4] *o2[2];
    const double Ly  = m1 * ( -  v1[2]*C1[0] +  v1[0]*C1[2]) + m2 * ( -  v2[2]*C2[0] +  v2[0]*C2[2]) + I1[3] *o1[0] + I1[1] *o1[1] + I1[5] *o1[2] + I2[3] *o2[0] + I2[1] *o2[1] + I2[5] *o2[2];
    const double Lz  = m1 * ( -  v1[0]*C1[1] +  v1[1]*C1[0]) + m2 * ( -  v2[0]*C2[1] +  v2[1]*C2[0]) + I1[4] *o1[0] + I1[5] *o1[1] + I1[2] *o1[2] + I2[4] *o2[0] + I2[5] *o2[1] + I2[2] *o2[2];
    const double hLx = m1 * ( - hv1[1]*C1[2] + hv1[2]*C1[1]) + m2 * ( - hv2[1]*C2[2] + hv2[2]*C2[1]) + I1[0]*ho1[0] + I1[3]*ho1[1] + I1[4]*ho1[2] + I2[0]*ho2[0] + I2[3]*ho2[1] + I2[4]*ho2[2];
    const double hLy = m1 * ( - hv1[2]*C1[0] + hv1[0]*C1[2]) + m2 * ( - hv2[2]*C2[0] + hv2[0]*C2[2]) + I1[3]*ho1[0] + I1[1]*ho1[1] + I1[5]*ho1[2] + I2[3]*ho2[0] + I2[1]*ho2[1] + I2[5]*ho2[2];
    const double hLz = m1 * ( - hv1[0]*C1[1] + hv1[1]*C1[0]) + m2 * ( - hv2[0]*C2[1] + hv2[1]*C2[0]) + I1[4]*ho1[0] + I1[5]*ho1[1] + I1[2]*ho1[2] + I2[4]*ho2[0] + I2[5]*ho2[1] + I2[2]*ho2[2];
    const double Mx = m1*v1[0] + m2*v2[0];
    const double My = m1*v1[1] + m2*v2[1];
    const double Mz = m1*v1[2] + m2*v2[2];
    const double hMx = m1*hv1[0] + m2*hv2[0];
    const double hMy = m1*hv1[1] + m2*hv2[1];
    const double hMz = m1*hv1[2] + m2*hv2[2];

    if (rank == 0)
    {
        std::cout << "vv = " << vv << std::endl;
        std::cout << "hvv = " << hvv << std::endl;

        std::cout << "E = " << E  << std::endl;
        std::cout << "hE= " << hE << std::endl;
        std::cout << " L= " << Lx << " " << Ly << " " << Lz << std::endl;
        std::cout << "hL= " << hLx << " " << hLy << " " << hLz << std::endl;
        std::cout << " M= " << Mx << " " << My << " " << Mz << std::endl;
        std::cout << "hM= " << hMx << " " << hMy << " " << hMz << std::endl;
    }
    #endif
}

}

void Penalization::preventCollidingObstacles() const
{
    using CHI_MAT = Real[FluidBlock::sizeZ][FluidBlock::sizeY][FluidBlock::sizeX];
    using UDEFMAT = Real[FluidBlock::sizeZ][FluidBlock::sizeY][FluidBlock::sizeX][3];

    const auto & shapes = sim.obstacle_vector->getObstacleVector();
    const auto & infos  = sim.grid->getBlocksInfo();
    const size_t N = sim.obstacle_vector->nObstacles();

    struct CollisionInfo // hitter and hittee, symmetry but we do things twice
    {
        Real iM = 0;
        Real iPosX = 0;
        Real iPosY = 0;
        Real iPosZ = 0;
        Real iMomX = 0;
        Real iMomY = 0;
        Real iMomZ = 0;
        Real ivecX = 0;
        Real ivecY = 0;
        Real ivecZ = 0;
        Real jM = 0;
        Real jPosX = 0;
        Real jPosY = 0;
        Real jPosZ = 0;
        Real jMomX = 0;
        Real jMomY = 0;
        Real jMomZ = 0;
        Real jvecX = 0;
        Real jvecY = 0;
        Real jvecZ = 0;
    };
    std::vector<CollisionInfo> collisions(N);

    std::vector <double> n_vec(3*N,0.0);

    #pragma omp parallel for schedule(static)
    for (size_t i=0; i<N; ++i)
    for (size_t j=0; j<N; ++j)
    {
        if(i==j) continue;
        auto & coll = collisions[i];

        const auto& iBlocks = shapes[i]->obstacleBlocks;
        const Real iU0      = shapes[i]->transVel[0];
        const Real iU1      = shapes[i]->transVel[1];
        const Real iU2      = shapes[i]->transVel[2];
        const Real iomega0  = shapes[i]->angVel  [0];
        const Real iomega1  = shapes[i]->angVel  [1];
        const Real iomega2  = shapes[i]->angVel  [2];
        const Real iCx      = shapes[i]->centerOfMass[0];
        const Real iCy      = shapes[i]->centerOfMass[1];
        const Real iCz      = shapes[i]->centerOfMass[2];

        const auto& jBlocks = shapes[j]->obstacleBlocks;
        const Real jU0      = shapes[j]->transVel[0];
        const Real jU1      = shapes[j]->transVel[1];
        const Real jU2      = shapes[j]->transVel[2];
        const Real jomega0  = shapes[j]->angVel  [0];
        const Real jomega1  = shapes[j]->angVel  [1];
        const Real jomega2  = shapes[j]->angVel  [2];
        const Real jCx      = shapes[j]->centerOfMass[0];
        const Real jCy      = shapes[j]->centerOfMass[1];
        const Real jCz      = shapes[j]->centerOfMass[2];

        assert(iBlocks.size() == jBlocks.size());

        const size_t nBlocks = iBlocks.size();
        for (size_t k=0; k<nBlocks; ++k)
        {
            if ( iBlocks[k] == nullptr || jBlocks[k] == nullptr ) continue;

            const auto & iSDF  = iBlocks[k]->sdfLab;
            const auto & jSDF  = jBlocks[k]->sdfLab;

            const CHI_MAT & iChi  = iBlocks[k]->chi;
            const CHI_MAT & jChi  = jBlocks[k]->chi;

            const UDEFMAT & iUDEF = iBlocks[k]->udef;
            const UDEFMAT & jUDEF = jBlocks[k]->udef;

            for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
            for(int ix=0; ix<FluidBlock::sizeX; ++ix)
            {
                if(iChi[iz][iy][ix] <= 0.0 || jChi[iz][iy][ix] <= 0.0 ) continue;

                const auto pos = infos[k].pos<Real>(ix, iy, iz);

                const Real iUr0 = iomega1* (pos[2] - iCz) - iomega2*(pos[1]-iCy);
                const Real iUr1 = iomega2* (pos[0] - iCx) - iomega0*(pos[2]-iCz);
                const Real iUr2 = iomega0* (pos[1] - iCy) - iomega1*(pos[0]-iCx);
                coll.iM    += iChi[iz][iy][ix];
                coll.iPosX += iChi[iz][iy][ix] * pos[0];
                coll.iPosY += iChi[iz][iy][ix] * pos[1];
                coll.iPosZ += iChi[iz][iy][ix] * pos[2];
                coll.iMomX += iChi[iz][iy][ix] * (iU0 + iUr0 + iUDEF[iz][iy][ix][0]);
                coll.iMomY += iChi[iz][iy][ix] * (iU1 + iUr1 + iUDEF[iz][iy][ix][1]);
                coll.iMomZ += iChi[iz][iy][ix] * (iU2 + iUr2 + iUDEF[iz][iy][ix][2]);

                const Real jUr0 = jomega1* (pos[2] - jCz) - jomega2*(pos[1]-jCy);
                const Real jUr1 = jomega2* (pos[0] - jCx) - jomega0*(pos[2]-jCz);
                const Real jUr2 = jomega0* (pos[1] - jCy) - jomega1*(pos[0]-jCx);
                coll.jM    += jChi[iz][iy][ix];
                coll.jPosX += jChi[iz][iy][ix] * pos[0];
                coll.jPosY += jChi[iz][iy][ix] * pos[1];
                coll.jPosZ += jChi[iz][iy][ix] * pos[2];
                coll.jMomX += jChi[iz][iy][ix] * (jU0 + jUr0 + jUDEF[iz][iy][ix][0]);
                coll.jMomY += jChi[iz][iy][ix] * (jU1 + jUr1 + jUDEF[iz][iy][ix][1]);
                coll.jMomZ += jChi[iz][iy][ix] * (jU2 + jUr2 + jUDEF[iz][iy][ix][2]);

                coll.ivecX += iChi[iz][iy][ix] * 0.5*(iSDF[iz+1][iy+1][ix+2] - iSDF[iz+1][iy+1][ix  ]);
                coll.ivecY += iChi[iz][iy][ix] * 0.5*(iSDF[iz+1][iy+2][ix+1] - iSDF[iz+1][iy  ][ix+1]);
                coll.ivecZ += iChi[iz][iy][ix] * 0.5*(iSDF[iz+2][iy+1][ix+1] - iSDF[iz  ][iy+1][ix+1]);

                coll.jvecX += jChi[iz][iy][ix] * 0.5*(jSDF[iz+1][iy+1][ix+2] - jSDF[iz+1][iy+1][ix  ]);
                coll.jvecY += jChi[iz][iy][ix] * 0.5*(jSDF[iz+1][iy+2][ix+1] - jSDF[iz+1][iy  ][ix+1]);
                coll.jvecZ += jChi[iz][iy][ix] * 0.5*(jSDF[iz+2][iy+1][ix+1] - jSDF[iz  ][iy+1][ix+1]);
            }
        }
    }

    std::vector<double> buffer(20*N); //CollisionInfo holds 20 doubles
    for (size_t i = 0 ; i < N ; i++)
    {
        auto & coll = collisions[i];
        buffer[20*i     ] = coll.iM   ;
        buffer[20*i + 1 ] = coll.iPosX;
        buffer[20*i + 2 ] = coll.iPosY;
        buffer[20*i + 3 ] = coll.iPosZ;
        buffer[20*i + 4 ] = coll.iMomX;
        buffer[20*i + 5 ] = coll.iMomY;
        buffer[20*i + 6 ] = coll.iMomZ;
        buffer[20*i + 7 ] = coll.ivecX;
        buffer[20*i + 8 ] = coll.ivecY;
        buffer[20*i + 9 ] = coll.ivecZ;
        buffer[20*i + 10] = coll.jM   ;
        buffer[20*i + 11] = coll.jPosX;
        buffer[20*i + 12] = coll.jPosY;
        buffer[20*i + 13] = coll.jPosZ;
        buffer[20*i + 14] = coll.jMomX;
        buffer[20*i + 15] = coll.jMomY;
        buffer[20*i + 16] = coll.jMomZ;
        buffer[20*i + 17] = coll.jvecX;
        buffer[20*i + 18] = coll.jvecY;
        buffer[20*i + 19] = coll.jvecZ;

    }
    MPI_Allreduce(MPI_IN_PLACE, buffer.data(), buffer.size(), MPI_DOUBLE, MPI_SUM, grid->getCartComm());
    for (size_t i = 0 ; i < N ; i++)
    {
        auto & coll = collisions[i];
        coll.iM    = buffer[20*i     ];
        coll.iPosX = buffer[20*i + 1 ];
        coll.iPosY = buffer[20*i + 2 ];
        coll.iPosZ = buffer[20*i + 3 ];
        coll.iMomX = buffer[20*i + 4 ];
        coll.iMomY = buffer[20*i + 5 ];
        coll.iMomZ = buffer[20*i + 6 ];
        coll.ivecX = buffer[20*i + 7 ];
        coll.ivecY = buffer[20*i + 8 ];
        coll.ivecZ = buffer[20*i + 9 ];
        coll.jM    = buffer[20*i + 10];
        coll.jPosX = buffer[20*i + 11];
        coll.jPosY = buffer[20*i + 12];
        coll.jPosZ = buffer[20*i + 13];
        coll.jMomX = buffer[20*i + 14];
        coll.jMomY = buffer[20*i + 15];
        coll.jMomZ = buffer[20*i + 16];
        coll.jvecX = buffer[20*i + 17];
        coll.jvecY = buffer[20*i + 18];
        coll.jvecZ = buffer[20*i + 19];
    }

    #pragma omp parallel for schedule(static)
    for (size_t i=0; i<N; ++i)
    for (size_t j=i+1; j<N; ++j)
    {
        if (i==j) continue;
        const double m1 = shapes[i]->mass;
        const double m2 = shapes[j]->mass;
        const double v1[3]={shapes[i]->transVel[0],shapes[i]->transVel[1],shapes[i]->transVel[2]};
        const double o1[3]={shapes[i]->  angVel[0],shapes[i]->  angVel[1],shapes[i]->  angVel[2]};
        const double v2[3]={shapes[j]->transVel[0],shapes[j]->transVel[1],shapes[j]->transVel[2]};
        const double o2[3]={shapes[j]->  angVel[0],shapes[j]->  angVel[1],shapes[j]->  angVel[2]};
        const double I1[6]={shapes[i]->J[0],shapes[i]->J[1],shapes[i]->J[2],shapes[i]->J[3],shapes[i]->J[4],shapes[i]->J[5]};
        const double I2[6]={shapes[j]->J[0],shapes[j]->J[1],shapes[j]->J[2],shapes[j]->J[3],shapes[j]->J[4],shapes[j]->J[5]};
        const double C1[3]={shapes[i]->centerOfMass[0],shapes[i]->centerOfMass[1],shapes[i]->centerOfMass[2]};
        const double C2[3]={shapes[j]->centerOfMass[0],shapes[j]->centerOfMass[1],shapes[j]->centerOfMass[2]};

        auto & coll       = collisions[i];
        auto & coll_other = collisions[j];
        // less than one fluid element of overlap: wait to get closer. no hit
        if(coll.iM       < 8.0 || coll.jM       < 8.0) continue; //object i did not collide
        if(coll_other.iM < 8.0 || coll_other.jM < 8.0) continue; //object j did not collide

        if (std::fabs(coll.iPosX/coll.iM  - coll_other.iPosX/coll_other.iM ) > 0.2 ||
            std::fabs(coll.iPosY/coll.iM  - coll_other.iPosY/coll_other.iM ) > 0.2 ||
            std::fabs(coll.iPosZ/coll.iM  - coll_other.iPosZ/coll_other.iM ) > 0.2 ) //used 0.2 because fish lenght is 0.2 usually!
        {
            continue; // then both objects i and j collided, but not with each other!
        }

        // A collision happened!
        sim.bCollision = true;

        const bool iForcedX = shapes[i]->bForcedInSimFrame[0];
        const bool iForcedY = shapes[i]->bForcedInSimFrame[1];
        const bool iForcedZ = shapes[i]->bForcedInSimFrame[2];
        const bool jForcedX = shapes[j]->bForcedInSimFrame[0];
        const bool jForcedY = shapes[j]->bForcedInSimFrame[1];
        const bool jForcedZ = shapes[j]->bForcedInSimFrame[2];
        if (iForcedX || iForcedY || iForcedZ || jForcedX || jForcedY || jForcedZ)
        {
            std::cout << "Forced objects not supported for collision." << std::endl;
            MPI_Abort(grid->getCartComm(),1);
        }

        double ho1[3];
        double ho2[3];
        double hv1[3];
        double hv2[3];

        //1. Compute collision normal vector (NX,NY,NZ)
        const Real norm_i = std::sqrt(coll.ivecX*coll.ivecX + coll.ivecY*coll.ivecY + coll.ivecZ*coll.ivecZ);
        const Real norm_j = std::sqrt(coll.jvecX*coll.jvecX + coll.jvecY*coll.jvecY + coll.jvecZ*coll.jvecZ);
        const Real mX = coll.ivecX/norm_i - coll.jvecX/norm_j;
        const Real mY = coll.ivecY/norm_i - coll.jvecY/norm_j;
        const Real mZ = coll.ivecZ/norm_i - coll.jvecZ/norm_j;
        const Real inorm = 1.0/std::sqrt(mX*mX + mY*mY + mZ*mZ);
        const Real NX = mX * inorm;
        const Real NY = mY * inorm;
        const Real NZ = mZ * inorm;

        //If objects are already moving away from each other, don't do anything
        //if( (v2[0]-v1[0])*NX + (v2[1]-v1[1])*NY + (v2[2]-v1[2])*NZ <= 0 ) continue;
        const Real hitVelX = coll.jMomX / coll.jM - coll.iMomX / coll.iM;
        const Real hitVelY = coll.jMomY / coll.jM - coll.iMomY / coll.iM;
        const Real hitVelZ = coll.jMomZ / coll.jM - coll.iMomZ / coll.iM;
        const Real projVel = hitVelX * NX + hitVelY * NY + hitVelZ * NZ;

        /*const*/ double vc1[3] = {coll.iMomX/coll.iM, coll.iMomY/coll.iM, coll.iMomZ/coll.iM};
        /*const*/ double vc2[3] = {coll.jMomX/coll.jM, coll.jMomY/coll.jM, coll.jMomZ/coll.jM};


        if(projVel<=0) continue; // vel goes away from collision: no need to bounce

        //2. Compute collision location
        const Real inv_iM = 1.0/coll.iM;
        const Real inv_jM = 1.0/coll.jM;
        const Real iPX = coll.iPosX * inv_iM; // object i collision location
        const Real iPY = coll.iPosY * inv_iM;
        const Real iPZ = coll.iPosZ * inv_iM;
        const Real jPX = coll.jPosX * inv_jM; // object j collision location
        const Real jPY = coll.jPosY * inv_jM;
        const Real jPZ = coll.jPosZ * inv_jM;
        const Real CX = 0.5*(iPX+jPX);
        const Real CY = 0.5*(iPY+jPY);
        const Real CZ = 0.5*(iPZ+jPZ);

        //3. Take care of the collision. Assume elastic collision (kinetic energy is conserved)
        ElasticCollision1(m1,m2,I1,I2,v1,v2,o1,o2,hv1,hv2,ho1,ho2,C1,C2,NX,NY,NZ,CX,CY,CZ,vc1,vc2);
        shapes[i]->transVel[0] = hv1[0];
        shapes[i]->transVel[1] = hv1[1];
        shapes[i]->transVel[2] = hv1[2];
        shapes[j]->transVel[0] = hv2[0];
        shapes[j]->transVel[1] = hv2[1];
        shapes[j]->transVel[2] = hv2[2];
        shapes[i]->angVel[0] = ho1[0];
        shapes[i]->angVel[1] = ho1[1];
        shapes[i]->angVel[2] = ho1[2];
        shapes[j]->angVel[0] = ho2[0];
        shapes[j]->angVel[1] = ho2[1];
        shapes[j]->angVel[2] = ho2[2];

        //if (sim.verbose)
        {
            #pragma omp critical
            {
                std::cout << "Collision between objects " << i << " and " << j << std::endl;
                std::cout << " iM   (0) = " << collisions[i].iM    << " jM   (1) = " << collisions[j].jM << std::endl;
                std::cout << " jM   (0) = " << collisions[i].jM    << " jM   (1) = " << collisions[j].iM << std::endl;
                std::cout << " Normal vector = (" << NX << "," << NY << "," << NZ << std::endl;
                std::cout << " Location      = (" << CX << "," << CY << "," << CZ << std::endl;
                std::cout << " Shape " << i << " before collision u    =(" <<  v1[0] << "," <<  v1[1] << "," <<  v1[2] << ")" << std::endl;
                std::cout << " Shape " << i << " after  collision u    =(" << hv1[0] << "," << hv1[1] << "," << hv1[2] << ")" << std::endl;
                std::cout << " Shape " << j << " before collision u    =(" <<  v2[0] << "," <<  v2[1] << "," <<  v2[2] << ")" << std::endl;
                std::cout << " Shape " << j << " after  collision u    =(" << hv2[0] << "," << hv2[1] << "," << hv2[2] << ")" << std::endl;
                std::cout << " Shape " << i << " before collision omega=(" <<  o1[0] << "," <<  o1[1] << "," <<  o1[2] << ")" << std::endl;
                std::cout << " Shape " << i << " after  collision omega=(" << ho1[0] << "," << ho1[1] << "," << ho1[2] << ")" << std::endl;
                std::cout << " Shape " << j << " before collision omega=(" <<  o2[0] << "," <<  o2[1] << "," <<  o2[2] << ")" << std::endl;
                std::cout << " Shape " << j << " after  collision omega=(" << ho2[0] << "," << ho2[1] << "," << ho2[2] << ")" << std::endl;
            }
        }
    }
}


Penalization::Penalization(SimulationData & s) : Operator(s) {}

void Penalization::operator()(const double dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  preventCollidingObstacles();

  std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
  #pragma omp parallel
  { // each thread needs to call its own non-const operator() function
    if(sim.bImplicitPenalization)
    {
      KernelPenalization<1> K(sim.dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
    else
    {
      KernelPenalization<0> K(sim.dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
  }

  ObstacleVisitor*K = new KernelFinalizePenalizationForce(sim.grid);
  sim.obstacle_vector->Accept(K); // accept you son of a french cow
  delete K;

  check("Penalization");
}

CubismUP_3D_NAMESPACE_END
