//
//  Cubism3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Operator.h"

#include "../SimulationData.h"
#include "../operators/Operator.h"
#include "../poisson/DiffusionSolverAMRKernels.h"
#include <Cubism/BlockInfo.h>
#include <vector>
#include <cassert>
#include <cstring>

CubismUP_3D_NAMESPACE_BEGIN

class DiffusionSolver
{
 public: 
  int mydirection = 0;
  Real dt;
 protected:
  SimulationData& sim;

  template <typename Lab>
  struct KernelLHSDiffusion
  {
    const SimulationData & sim;
    KernelLHSDiffusion(const SimulationData&s, const Real _dt) : sim(s), dt(_dt) {}
    const std::vector<BlockInfo> & lhsInfo = sim.lhsInfo();
    const Real dt;
    const int Nx = ScalarBlock::sizeX;
    const int Ny = ScalarBlock::sizeY;
    const int Nz = ScalarBlock::sizeZ;
    const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};

    void operator()(const Lab & lab, const BlockInfo& info) const
    {
      ScalarBlock & __restrict__ o = (*sim.lhs)(info.blockID);
      const Real h = info.h; 
      const Real coef = -1.0/(dt*sim.nu)*h*h*h;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        o(x,y,z) = h*( lab(x-1,y,z) + lab(x+1,y,z) + lab(x,y-1,z) + lab(x,y+1,z) +lab(x,y,z-1) + lab(x,y,z+1) - 6.0*lab(x,y,z)) + coef*lab(x,y,z);
      }

      BlockCase<ScalarBlock> * tempCase = (BlockCase<ScalarBlock> *)(lhsInfo[info.blockID].auxiliary);

      if (tempCase == nullptr) return; //no flux corrections needed for this block

      ScalarElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
      ScalarElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
      ScalarElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
      ScalarElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
      ScalarElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
      ScalarElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;

      if (faceXm != nullptr)
      {
        const int x = 0;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
          faceXm[y + Ny*z] = h*(lab(x,y,z) - lab(x-1,y,z));
      }
      if (faceXp != nullptr)
      {
        const int x = Nx-1;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
          faceXp[y + Ny*z] = h*(lab(x,y,z) - lab(x+1,y,z));
      }
      if (faceYm != nullptr)
      {
        const int y = 0;
        for(int z=0; z<Nz; ++z)
        for(int x=0; x<Nx; ++x)
          faceYm[x + Nx*z] = h*(lab(x,y,z) - lab(x,y-1,z));
      }
      if (faceYp != nullptr)
      {
        const int y = Ny-1;
        for(int z=0; z<Nz; ++z)
        for(int x=0; x<Nx; ++x)
          faceYp[x + Nx*z] = h*(lab(x,y,z) - lab(x,y+1,z));
      }
      if (faceZm != nullptr)
      {
        const int z = 0;
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          faceZm[x + Nx*y] = h*(lab(x,y,z) - lab(x,y,z-1));
      }
      if (faceZp != nullptr)
      {
        const int z = Nz-1;
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          faceZp[x + Nx*y] = h*(lab(x,y,z) - lab(x,y,z+1));
      }

    }
  };

  void _preconditioner(const std::vector<Real> & input, std::vector<Real> & output)
  {
    auto &  zInfo         = sim.pres->getBlocksInfo(); //used for preconditioning
    const size_t Nblocks  = zInfo.size();
    const int BSX         = VectorBlock::sizeX;
    const int BSY         = VectorBlock::sizeY;
    const int BSZ         = VectorBlock::sizeZ;

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      ScalarBlock & __restrict__ bb = (*sim.pres)(i);
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        bb(ix,iy,iz).s = input[j];
      }
    }

    #pragma omp parallel
    {
      cubismup3d::diffusion_kernels::getZImplParallel(sim.presInfo(),sim.nu,dt);
    }

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      const ScalarBlock & __restrict__ bb = (*sim.pres)(i);
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        output[j] = bb(ix,iy,iz).s;;
      }
    }

  }

  void _lhs(std::vector<Real> & input, std::vector<Real> & output)
  {
    auto &  zInfo         = sim.pres->getBlocksInfo(); //used for preconditioning
    auto & AxInfo         = sim.lhs ->getBlocksInfo(); //will store the LHS result
    const size_t Nblocks  = zInfo.size();
    const int BSX         = VectorBlock::sizeX;
    const int BSY         = VectorBlock::sizeY;
    const int BSZ         = VectorBlock::sizeZ;

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      ScalarBlock & __restrict__ zz = *(ScalarBlock*) zInfo[i].ptrBlock;
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        zz(ix,iy,iz).s  = input[j];
      }
    }

    using Lab0 = cubism::BlockLabMPI<BlockLabBC<ScalarBlock, aligned_block_allocator,0>, ScalarGrid>;
    using Lab1 = cubism::BlockLabMPI<BlockLabBC<ScalarBlock, aligned_block_allocator,1>, ScalarGrid>;
    using Lab2 = cubism::BlockLabMPI<BlockLabBC<ScalarBlock, aligned_block_allocator,2>, ScalarGrid>;

    if (mydirection == 0) compute<Lab0>(KernelLHSDiffusion<Lab0>(sim,dt),sim.pres,sim.lhs);
    if (mydirection == 1) compute<Lab1>(KernelLHSDiffusion<Lab1>(sim,dt),sim.pres,sim.lhs);
    if (mydirection == 2) compute<Lab2>(KernelLHSDiffusion<Lab2>(sim,dt),sim.pres,sim.lhs);

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      ScalarBlock & __restrict__ Ax = *(ScalarBlock*) AxInfo[i].ptrBlock;
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        output[j]   = Ax(ix,iy,iz).s;
      }
    }
  }

  std::vector<Real> b;
  std::vector<Real> phat;
  std::vector<Real> rhat;
  std::vector<Real> shat;
  std::vector<Real> what;
  std::vector<Real> zhat;
  std::vector<Real> qhat;
  std::vector<Real> s;
  std::vector<Real> w;
  std::vector<Real> z;
  std::vector<Real> t;
  std::vector<Real> v;
  std::vector<Real> q;
  std::vector<Real> r;
  std::vector<Real> y;
  std::vector<Real> x;
  std::vector<Real> r0;
  std::vector<Real> x_opt;

 public:
  DiffusionSolver(SimulationData& ss): sim(ss){}
  DiffusionSolver(const DiffusionSolver& c) = delete; 

  void solve()
  {
    //Algorithm 11 from the paper:
    //"The communication-hiding pipelined BiCGstab method for the parallel solution of large unsymmetric linear systems"
    //by S. Cools, W. Vanroose
    //This is a BiCGstab with less global communication (reductions) that are overlapped with computation.
  
    //Warning: 'input'  initially contains the RHS of the system!
    //Warning: 'output' initially contains the initial solution guess x0!
    const auto & AxInfo      = sim.lhsInfo();//input ->getBlocksInfo(); //will store the LHS result
    const auto &  zInfo      = sim.presInfo();//output->getBlocksInfo(); //will store the input 'x' when LHS is computed
    const size_t Nblocks     = zInfo.size();            //total blocks of this rank
    const int BSX            = VectorBlock::sizeX;      //block size in x direction
    const int BSY            = VectorBlock::sizeY;      //block size in y direction
    const int BSZ            = VectorBlock::sizeZ;      //block size in z direction
    const size_t N           = BSX*BSY*BSZ*Nblocks;     //total number of variables of this rank
    const Real eps           = 1e-100;                  //used in denominators, to not divide by zero
    const Real max_error     = sim.DiffusionErrorTol;          //error tolerance for Linf norm of residual
    const Real max_rel_error = sim.DiffusionErrorTolRel;       //relative error tolerance for Linf(r)/Linf(r0)
    bool serious_breakdown   = false; //shows if the solver will restart in this iteration
    bool useXopt             = false; //(is almost always true) use the solution that had the smallest residual
    Real min_norm            = 1e50;  //residual norm (for best solution, see also 'useXopt')
    Real norm_1              = 0.0;   //used to decide if the solver will restart 
    Real norm_2              = 0.0;   //used to decide if the solver will restart
    const MPI_Comm m_comm    = sim.comm;
    const bool verbose       = sim.rank == 0;
  
    phat.resize(N);
    rhat.resize(N);
    shat.resize(N);
    what.resize(N);
    zhat.resize(N);
    qhat.resize(N);
    s   .resize(N);
    w   .resize(N);
    z   .resize(N);
    t   .resize(N);
    v   .resize(N);
    q   .resize(N);
    r   .resize(N);
    y   .resize(N);
    x   .resize(N);
    r0  .resize(N);
    b   .resize(N); // RHS of the system will be stored here
    x_opt.resize(N);// solution with minimum residual
  
    //initialize b,r,x
    #pragma omp parallel for
    for(size_t i=0; i< Nblocks; i++)
    {    
      ScalarBlock & __restrict__ rhs  = *(ScalarBlock*) AxInfo[i].ptrBlock;
      const ScalarBlock & __restrict__ zz = *(ScalarBlock*)  zInfo[i].ptrBlock;
  
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        b[j] = rhs(ix,iy,iz).s;
        r[j] = rhs(ix,iy,iz).s;
        x[j] = zz (ix,iy,iz).s;
      }
    }
  
  
    //In what follows, we indicate by (*n*) the n-th step of the algorithm
  
    //(*2*) r0 = b - A*x0, r0hat = M^{-1}*r0, w0=A*r0hat, w0hat=M^{-1}w0
    _lhs(x,r0);
    #pragma omp parallel for
    for (size_t i=0; i < N; i++)
    {
      r0[i] = r [i] - r0[i];
      r [i] = r0[i];
    }
    _preconditioner(r0,rhat);
    _lhs(rhat,w);
    _preconditioner(w,what);
  
    //(*3*) t0=A*w0hat, alpha0 = (r0,r0) / (r0,w0), beta=0
    _lhs(what,t);
    Real alpha = 0.0;
    Real norm  = 0.0;
    Real beta  = 0.0;
    Real omega = 0.0;
    Real r0r_prev;
    {
      Real temp0 = 0.0;
      Real temp1 = 0.0;
      #pragma omp parallel for reduction (+:temp0,temp1,norm)
      for (size_t j=0; j < N; j++)
      {
        temp0 += r0[j]*r0[j];
        temp1 += r0[j]*w [j];
          norm += r0[j]*r0[j];//std::max(norm,std::fabs(r[j]));
        }
        Real temporary[2] = {temp0,temp1};
        MPI_Allreduce(MPI_IN_PLACE,temporary,2,MPI_Real,MPI_SUM,m_comm);
        MPI_Allreduce(MPI_IN_PLACE,&norm    ,1,MPI_Real,MPI_SUM,m_comm);
        alpha = temporary[0]/(temporary[1]+eps);
        r0r_prev = temporary[0];
        norm = std::sqrt(norm);
        if (verbose) std::cout << "[Diffusion solver, direction: " << mydirection << " ]: initial error norm:" << norm << "\n";
      }
      const Real init_norm = norm;
    
      //(*4*) for k=0,1,...
      int k;
      for ( k = 0 ; k < 1000; k++)
      {
        Real qy = 0.0;
        Real yy = 0.0;
    
        //(*5*),(*6*),...,(*11*)
        if (k%50 != 0)
        { 
          #pragma omp parallel for reduction (+:qy,yy)
          for (size_t j=0; j < N; j++)
          {
            phat[j] = rhat[j] + beta * (phat[j] - omega * shat[j]);
            s   [j] = w   [j] + beta * (s   [j] - omega * z   [j]);
            shat[j] = what[j] + beta * (shat[j] - omega * zhat[j]);
            z   [j] = t   [j] + beta * (z   [j] - omega * v   [j]);
            q   [j] = r   [j] - alpha*  s   [j];
            qhat[j] = rhat[j] - alpha*  shat[j];
            y   [j] = w   [j] - alpha*  z   [j];
            qy += q[j]*y[j];
            yy += y[j]*y[j];
          }
        }
        else
        {
          //every 50 iterations we use the residual replacement strategy, to prevent loss of accuracy
          //and compute stuff with the exact (not pipelined) versions
          #pragma omp parallel for
          for (size_t j=0; j < N; j++)
          {
            phat[j] = rhat[j] + beta * (phat[j] - omega * shat[j]);
          }
          _lhs(phat,s);
          _preconditioner(s,shat);
          _lhs(shat,z);
          #pragma omp parallel for
          for (size_t j=0; j < N; j++)
          {
            q   [j] = r   [j] - alpha*  s   [j];
            qhat[j] = rhat[j] - alpha*  shat[j];
            y   [j] = w   [j] - alpha*  z   [j];
            qy += q[j]*y[j];
            yy += y[j]*y[j];
          }
        }
    
        //(*12*) begin reduction (q,y),(y,y)
        MPI_Request request;
        Real quantities[6];
        quantities[0] = qy;
        quantities[1] = yy;
        MPI_Iallreduce(MPI_IN_PLACE,&quantities,2,MPI_Real,MPI_SUM,m_comm,&request);
    
        //(*13*) computation zhat = M^{-1}*z
        _preconditioner(z,zhat);
    
        //(*14*) computation v = A*zhat
        _lhs(zhat,v);
    
        //(*15*) end reduction
        MPI_Waitall(1,&request,MPI_STATUSES_IGNORE);
        qy = quantities[0];
        yy = quantities[1];
    
        //(*16*) omega = (q,y)/(y,y)
        omega = qy / (yy+eps);
    
        //(*17*),(*18*),(*19*),(*20*)
        Real r0r = 0.0;
        Real r0w = 0.0;
        Real r0s = 0.0;
        Real r0z = 0.0;
        norm = 0.0;
        norm_1 = 0.0;
        norm_2 = 0.0;
        if (k%50 != 0)
        {
          #pragma omp parallel for reduction (+:r0r,r0w,r0s,r0z,norm_1,norm_2,norm)
          for (size_t j=0; j < N; j++)
          {
            x   [j] = x   [j] + alpha *  phat[j] + omega * qhat[j] ;
            r   [j] = q   [j] - omega *  y   [j]                   ;
            rhat[j] = qhat[j] - omega * (what[j] - alpha * zhat[j]);
            w   [j] = y   [j] - omega * (t   [j] - alpha * v   [j]);
            r0r += r0[j]*r[j];
            r0w += r0[j]*w[j];
            r0s += r0[j]*s[j];
            r0z += r0[j]*z[j];
            norm += r[j]*r[j];//std::max(norm,std::fabs(r[j]));
            norm_1 += r [j] * r [j];
            norm_2 += r0[j] * r0[j];
          }
        }
        else
        {
          //every 50 iterations we use the residual replacement strategy, to prevent loss of accuracy
          //and compute stuff with the exact (not pipelined) versions
          #pragma omp parallel for
          for (size_t j=0; j < N; j++)
          {
            x   [j] = x   [j] + alpha *  phat[j] + omega * qhat[j] ;
          }
          _lhs(x,r);
          #pragma omp parallel for
          for (size_t j=0; j < N; j++)
          {
            r[j] = b[j] - r[j] ;
          }
          _preconditioner(r,rhat);
          _lhs(rhat,w);
          #pragma omp parallel for reduction (+:r0r,r0w,r0s,r0z,norm_1,norm_2,norm)
          for (size_t j=0; j < N; j++)
          {
            r0r += r0[j]*r[j];
            r0w += r0[j]*w[j];
            r0s += r0[j]*s[j];
            r0z += r0[j]*z[j];
            norm += r[j]*r[j];
            //norm = std::max(norm,std::fabs(r[j]));
            norm_1 += r [j] * r [j];
            norm_2 += r0[j] * r0[j];
          }
        }
        quantities[0] = r0r;
        quantities[1] = r0w;
        quantities[2] = r0s;
        quantities[3] = r0z;
        quantities[4] = norm_1;
        quantities[5] = norm_2;
    
        //(*21*) begin reductions
        MPI_Request request2;
        MPI_Iallreduce(MPI_IN_PLACE,&quantities,6,MPI_Real,MPI_SUM,m_comm,&request );
        MPI_Iallreduce(MPI_IN_PLACE,&norm      ,1,MPI_Real,MPI_SUM,m_comm,&request2);
    
        //(*22*) computation what = M^{-1}*w
        _preconditioner(w,what);
    
        //(*23*) computation t = A*what
        _lhs(what,t);
    
        //(*24*) end reductions
        MPI_Waitall(1,&request,MPI_STATUSES_IGNORE);
        r0r = quantities[0];
        r0w = quantities[1];
        r0s = quantities[2];
        r0z = quantities[3];
        norm_1 = quantities[4];
        norm_2 = quantities[5];
    
        //(*25*)
        beta   = alpha / (omega + eps)* r0r / (r0r_prev + eps);
    
        //(*26*)
        alpha  = r0r / (r0w+beta*r0s-beta*omega*r0z);
        Real alphat = 1.0/(omega+eps) + r0w/(r0r+eps) - beta*omega*r0z/(r0r+eps);
        alphat = 1.0/(alphat+eps);
        if (std::fabs(alphat) < 10*std::fabs(alpha)) alpha = alphat;
    
        r0r_prev = r0r;
        MPI_Waitall(1,&request2,MPI_STATUSES_IGNORE);
        norm = std::sqrt(norm);
        //Check if restart should be made. If so, current solution estimate is used as an initial
        //guess and solver starts again.
        serious_breakdown = r0r * r0r < 1e-16 * norm_1 * norm_2;
        if (serious_breakdown)// && restarts < max_restarts)
        {
          //restarts ++;
          if (verbose)
            std::cout << "  [Diffusion solver, direction: " << mydirection << " ]: Restart at iteration: " << k << " norm: " << norm << std::endl;
    
          #pragma omp parallel for
          for (size_t i=0; i < N; i++)
            r0[i] = r[i];
        
          _preconditioner(r0,rhat);
          _lhs(rhat,w);
        
          alpha = 0.0;
          Real temp0 = 0.0;
          Real temp1 = 0.0;
          #pragma omp parallel for reduction (+:temp0,temp1)
          for (size_t j=0; j < N; j++)
          {
            temp0 += r0[j]*r0[j];
            temp1 += r0[j]*w [j];
          }
          Real temporary[2] = {temp0,temp1};
          MPI_Iallreduce(MPI_IN_PLACE, temporary,2,MPI_Real,MPI_SUM,m_comm,&request2);
        
          _preconditioner(w,what);
          _lhs(what,t);
        
          MPI_Waitall(1,&request2,MPI_STATUSES_IGNORE);
        
          alpha = temporary[0]/(temporary[1]+eps);
          r0r_prev = temporary[0];
          beta  = 0.0;
          omega = 0.0;
        }
    
        if (norm < min_norm)
        {
          useXopt = true;
          min_norm = norm;
          #pragma omp parallel for
          for (size_t i=0; i < N; i++)
            x_opt[i] = x[i];
        }
        //if (verbose) std::cout << k <<"  norm: " << norm << "\n";
    
    
        if ( norm < max_error || norm/(init_norm+eps) < max_rel_error )
        {
          if (verbose)
            std::cout << "  [Diffusion solver, direction: " << mydirection << " ]: Converged after " << k << " iterations.\n";
          break;
        }
      }
    
      if (verbose)
      {
        std::cout <<  " Error norm (relative) = " << min_norm << "/" << max_error << std::endl;
      }
    
      Real * xsol = useXopt ? x_opt.data():x.data();
      #pragma omp parallel for
      for (size_t i=0; i < Nblocks; i++)
      {
        ScalarBlock & __restrict__ bb  = (*sim.pres)(i);
        for(int iz=0; iz<BSZ; iz++)
        for(int iy=0; iy<BSY; iy++)
        for(int ix=0; ix<BSX; ix++)
        {
            const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
            bb(ix,iy,iz).s = xsol[j];
        }
      }
    }
};


class AdvectionDiffusionImplicit : public Operator
{
  std::vector<Real> pressure;
  std::vector<Real> velocity;
public:
  AdvectionDiffusionImplicit(SimulationData & s) : Operator(s) { }

  ~AdvectionDiffusionImplicit() { }

  void operator()(const Real dt);
  void euler(const Real dt);

  std::string getName() { return "AdvectionDiffusionImplicit"; }
};

CubismUP_3D_NAMESPACE_END
