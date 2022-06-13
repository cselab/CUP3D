//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PoissonSolverAMR.h"
#include "PoissonSolverAMRKernels.h"

namespace cubismup3d {

PoissonSolverAMR::PoissonSolverAMR(SimulationData& s): sim(s),findLHS(s){}

static void getZImplParallel(
    SimulationData &sim,
    const std::vector<cubism::BlockInfo>& vInfo)
{
  using namespace cubismup3d::poisson_kernels;
  const size_t Nblocks = vInfo.size();

  // We could enable this, we don't really care about denormals.
  // _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  // A struct to enforce relative alignment between matrices. The relative
  // alignment of Ax and r MUST NOT be a multiple of 4KB due to cache bank
  // conflicts. See "Haswell and Broadwell pipeline", section
  // "Cache and memory access" here:
  // https://www.agner.org/optimize/microarchitecture.pdf
  struct Tmp {
    // It seems like some offsets with respect to the page boundary of 4KB are
    // faster than the others. (This is accomplished by adding an offset here
    // and using alignas(4096) below). However, this is likely CPU-dependent,
    // so we don't hardcode such fine-tunings here.
    // char offset[0xec0];
    Block r;
    // Ensure p[0+1][0+1][0+xPad] is 64B-aligned for AVX-512 to work.
    char padding1[64 - xPad * sizeof(Real)];
    PaddedBlock p;
    char padding2[xPad * sizeof(Real)];
    Block Ax;
  };
  alignas(64) Tmp tmp{};  // See the kernels cpp file for required alignments.
  Block &r = tmp.r;
  Block &Ax = tmp.Ax;
  PaddedBlock &p = tmp.p;

  #pragma omp for
  for (size_t i = 0; i < Nblocks; ++i) {
    static_assert(sizeof(ScalarBlock) == sizeof(Block));
    assert((uintptr_t)vInfo[i].ptrBlock % kBlockAlignment == 0);
    Block &block = *(Block *)__builtin_assume_aligned(
        vInfo[i].ptrBlock, kBlockAlignment);

    const Real invh = 1 / vInfo[i].h;
    Real rrPartial[NX] = {};
    for (int iz = 0; iz < NZ; ++iz)
    for (int iy = 0; iy < NY; ++iy)
    for (int ix = 0; ix < NX; ++ix) {
      r[iz][iy][ix] = invh * block[iz][iy][ix];
      rrPartial[ix] += r[iz][iy][ix] * r[iz][iy][ix];
      p[iz + 1][iy + 1][ix + xPad] = r[iz][iy][ix];
      block[iz][iy][ix] = 0;
    }
    Real rr = sum(rrPartial);

    const Real sqrNorm0 = (Real)1 / (N * N) * rr;

    if (sqrNorm0 < 1e-32)
      continue;

    const Real *pW = &p[0][0][0] - 1;
    const Real *pE = &p[0][0][0] + 1;

    for (int k = 0; k < 100; ++k) {
      // rr = kernelPoissonGetZInnerReference(p,Ax, r, block, sqrNorm0, rr);
      rr = kernelPoissonGetZInner(p, pW, pE, Ax, r, block, sqrNorm0, rr);
      if (rr <= 0)
        break;
    }
  }
}


void PoissonSolverAMR::getZ()
{
  #pragma omp parallel
  {
    getZImplParallel(sim, sim.presInfo());
  }
}

void PoissonSolverAMR::solve()
{
    const std::vector<BlockInfo>& vInfo_lhs = sim.lhsInfo();
    const std::vector<BlockInfo>& vInfo_z   = sim.presInfo();
    const Real eps = 1e-21;
    const size_t Nblocks = vInfo_z.size();
    const size_t N = BlockType::sizeX*BlockType::sizeY*BlockType::sizeZ*Nblocks;
    const int Nx = ScalarBlock::sizeX;
    const int Ny = ScalarBlock::sizeY;
    const int Nz = ScalarBlock::sizeZ;
    std::vector<Real> x   (N,0.0);
    std::vector<Real> r   (N,0.0);
    std::vector<Real> p   (N,0.0);
    std::vector<Real> v   (N,0.0);
    std::vector<Real> s   (N,0.0);
    std::vector<Real> rhat(N,0.0);

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
        const ScalarBlock & __restrict__ Z    = *(ScalarBlock*) vInfo_z  [i].ptrBlock;
        ScalarBlock & __restrict__ LHS  = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
	if (sim.bMeanConstraint == 1 || sim.bMeanConstraint > 2)
          if (vInfo_z[i].index[0] == 0 && vInfo_z[i].index[1] == 0 && vInfo_z[i].index[2] == 0) LHS(0,0,0).s = 0.0;
        for(int iz=0; iz<Nz; iz++)
        for(int iy=0; iy<Ny; iy++)
        for(int ix=0; ix<Nx; ix++)
        {
            const size_t src_index = _dest(vInfo_lhs[i], iz, iy, ix);
            x[src_index] = Z(ix,iy,iz).s;
            r[src_index] = LHS(ix,iy,iz).s;//this is the rhs now
        }
    }


    //1. rVector <-- b - AxVector
    //2. rhatVector = rVector
    //3. rho = a = omega = 1.0
    //4. vVector = pVector = 0

    Real norm = 0.0; 
    findLHS(0); // AxVector <-- A*x_{0}, x_0 = pressure
    #pragma omp parallel for reduction (+:norm)
    for(size_t i=0; i< Nblocks; i++)
    {
        const ScalarBlock & __restrict__ LHS  = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
        for(int iz=0; iz<Nz; iz++)
        for(int iy=0; iy<Ny; iy++)
        for(int ix=0; ix<Nx; ix++)
        {
            const size_t src_index = _dest(vInfo_z[i], iz, iy, ix);
            r[src_index] -= LHS(ix,iy,iz).s;
            rhat[src_index] = r[src_index];
            norm+= r[src_index]*r[src_index];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_Real,MPI_SUM,sim.comm);
    norm = std::sqrt(norm);

    Real rho = 1.0;
    Real alpha = 1.0;
    Real omega = 1.0;
    Real rho_m1;

    std::vector <Real> x_opt (N);
    bool useXopt = false;
    Real min_norm = 1e50;
    Real init_norm=norm;
    const Real max_error     = sim.PoissonErrorTol;
    const Real max_rel_error = sim.PoissonErrorTolRel;
    bool serious_breakdown = false;
    int iter_opt = 0;

    //5. for k = 1,2,...
    for (size_t k = 1; k < 1000; k++)
    {
        //1. rho_i = (rhat_0,r_{k-1})
        //2. beta = rho_{i}/rho_{i-1} * alpha/omega_{i-1}
        rho_m1 = rho;
        rho = 0.0;
        Real norm_1 = 0.0;
        Real norm_2 = 0.0;
        #pragma omp parallel for reduction(+:rho,norm_1,norm_2)
        for(size_t i=0; i< N; i++)
        {
            rho += r[i]*rhat[i];
            norm_1 += r[i]*r[i];
            norm_2 += rhat[i]*rhat[i];
        }
        Real aux_norm [3] = {rho,norm_1,norm_2};
        MPI_Allreduce(MPI_IN_PLACE,&aux_norm,3,MPI_Real,MPI_SUM,sim.comm);
        rho = aux_norm[0];
        norm_1 = aux_norm[1];
        norm_2 = aux_norm[2];
        Real beta = rho / (rho_m1+eps) * alpha / (omega+eps) ;

        norm_1 = sqrt(norm_1);
        norm_2 = sqrt(norm_2);
        Real cosTheta = rho/norm_1/norm_2; 
        serious_breakdown = std::fabs(cosTheta) < 1e-8;
        if (serious_breakdown)
        {
            beta = 0.0;
            rho = 0.0;
            #pragma omp parallel for reduction(+:rho)
            for(size_t i=0; i< N; i++)
            {
                rhat[i] = r[i];
                rho += r[i]*rhat[i];
                p[i] = 0.0;
                v[i] = 0.0;
            }
            if (sim.rank==0) 
                std::cout << "  [Poisson solver]: restart at iteration:" << k << 
                             "  norm:"<< norm <<" init_norm:" << init_norm << std::endl;
            MPI_Allreduce(MPI_IN_PLACE,&rho,1,MPI_Real,MPI_SUM,sim.comm);
            alpha = 1.;
            omega = 1.;
            rho_m1 = 1.;
            beta = rho / (rho_m1+eps) * alpha / (omega+eps) ;
        }

        //3. p_i = r_{i-1} + beta*(p_{i-1}-omega *v_{i-1})
        //4. z = K_2^{-1} p
        #pragma omp parallel for
        for (size_t i=0; i < Nblocks; i++)
        {
            ScalarBlock & __restrict__ Z    = *(ScalarBlock*) vInfo_z  [i].ptrBlock;
            for(int iz=0; iz<Nz; iz++)
            for(int iy=0; iy<Ny; iy++)
            for(int ix=0; ix<Nx; ix++)
            {
                const size_t src_index = _dest(vInfo_z[i], iz, iy, ix);
                p[src_index] = r[src_index] + beta*(p[src_index]-omega*v[src_index]);
                Z(ix,iy,iz).s = p[src_index];
            }
        }
        getZ();

        //5. v = A z
        //6. alpha = rho_i / (rhat_0,v_i)
        alpha = 0.0;
        findLHS(0); //v stored in AxVector
        #pragma omp parallel for reduction(+:alpha)
        for (size_t i=0; i < Nblocks; i++)
        {
            const ScalarBlock & __restrict__ LHS  = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
            for(int iz=0; iz<Nz; iz++)
            for(int iy=0; iy<Ny; iy++)
            for(int ix=0; ix<Nx; ix++)
            {
                const size_t src_index = _dest(vInfo_lhs[i], iz, iy, ix);
                v[src_index] = LHS(ix,iy,iz).s;
                alpha += rhat[src_index] * v[src_index];
            }
        }  
        MPI_Allreduce(MPI_IN_PLACE,&alpha,1,MPI_Real,MPI_SUM,sim.comm);
        alpha = rho / (alpha + eps);
        //7. x += a z
        //8. 
        //9. s = r_{i-1}-alpha * v_i
        //10. z = K_2^{-1} s
        #pragma omp parallel for
        for (size_t i=0; i < Nblocks; i++)
        {
            ScalarBlock & __restrict__ Z    = *(ScalarBlock*) vInfo_z  [i].ptrBlock;
            for(int iz=0; iz<Nz; iz++)
            for(int iy=0; iy<Ny; iy++)
            for(int ix=0; ix<Nx; ix++)
            {
                const size_t src_index = _dest(vInfo_z[i], iz, iy, ix);
                x[src_index] += alpha * Z(ix,iy,iz).s;
                s[src_index] = r[src_index] - alpha * v[src_index];
                Z(ix,iy,iz).s = s[src_index];
            }
        }
        getZ();

        //11. t = Az
        findLHS(0); // t stored in AxVector

        //12. omega = ...
        Real aux1 = 0;
        Real aux2 = 0;
        #pragma omp parallel for reduction (+:aux1,aux2)
        for (size_t i=0; i < Nblocks; i++)
        {
            const ScalarBlock & __restrict__ LHS  = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
            for(int iz=0; iz<Nz; iz++)
            for(int iy=0; iy<Ny; iy++)
            for(int ix=0; ix<Nx; ix++)
            {
                const size_t src_index = _dest(vInfo_lhs[i], iz, iy, ix);
                aux1 += LHS(ix,iy,iz).s * s[src_index];
                aux2 += LHS(ix,iy,iz).s * LHS(ix,iy,iz).s;
            }
        }
        Real aux_12[2] = {aux1,aux2};
        MPI_Allreduce(MPI_IN_PLACE,&aux_12,2,MPI_Real,MPI_SUM,sim.comm);
        aux1 = aux_12[0];
        aux2 = aux_12[1];
        omega = aux1 / (aux2+eps); 

        //13. x += omega * z
        //14.
        //15. r = s - omega * t
        norm = 0.0; 
        #pragma omp parallel for reduction(+:norm)
        for (size_t i=0; i < Nblocks; i++)
        {
            const ScalarBlock & __restrict__ LHS  = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
            const ScalarBlock & __restrict__ Z    = *(ScalarBlock*) vInfo_z  [i].ptrBlock;
            for(int iz=0; iz<Nz; iz++)
            for(int iy=0; iy<Ny; iy++)
            for(int ix=0; ix<Nx; ix++)
            {
                const size_t src_index = _dest(vInfo_z[i], iz, iy, ix);
                x[src_index] += omega * Z(ix,iy,iz).s;
                r[src_index] = s[src_index]-omega*LHS(ix,iy,iz).s;
                norm+= r[src_index]*r[src_index];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_Real,MPI_SUM,sim.comm);
        norm = std::sqrt(norm);

        if (norm < min_norm)
        {
            useXopt = true;
            iter_opt = k;
            min_norm = norm;
            #pragma omp parallel for
            for (size_t i=0; i < N; i++)
            {
                x_opt[i] = x[i];
            }
        }
  
        if (k==1) init_norm = norm;

        if (norm / (init_norm+eps) > 100000.0 && k > 10)
        {
            useXopt = true;
            if (sim.rank==0)
                std::cout <<  "XOPT Poisson solver converged after " <<  k << " iterations. Error norm = " << norm << "  iter_opt="<< iter_opt << std::endl;
            break;
        }

        if ( (norm < max_error || norm/init_norm < max_rel_error ) )
        {
            if (sim.rank==0)
                std::cout <<  "  [Poisson solver]: Converged after " <<  k << " iterations. Error norm = " << norm << std::endl;
            break;
        }

    }//k-loop

    Real * xsol = useXopt ? x_opt.data():x.data();
    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
        const int m = vInfo_z[i].level;
        const long long n = vInfo_z[i].Z;
        const BlockInfo & info = sim.pres->getBlockInfoAll(m,n);
        ScalarBlock & __restrict__ b  = *(ScalarBlock*) info.ptrBlock;
        for(int iz=0; iz<Nz; iz++)
        for(int iy=0; iy<Ny; iy++)
        for(int ix=0; ix<Nx; ix++)
        {
            const size_t src_index = _dest(vInfo_z[i], iz, iy, ix);
            b(ix,iy,iz).s = xsol[src_index];
        }
    }
  
}
}//namespace cubismup3d
