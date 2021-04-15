//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PoissonSolverAMR.h"

namespace cubismup3d {

Real PoissonSolverAMR::computeAverage() const
{
  Real avgP = 0;
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static) reduction(+ : avgP)
  for(size_t i=0; i<vInfo.size(); ++i)
  {
    const double h3 = vInfo[i].h_gridpoint * vInfo[i].h_gridpoint * vInfo[i].h_gridpoint;
    BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++)
      avgP += h3 * b.data[iz][iy][ix].p;
  }
  MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPIREAL, MPI_SUM, m_comm);
  avgP /= sim.extent[0] * sim.extent[1] * sim.extent[2];
  return avgP;
}

#ifdef PRECOND
double PoissonSolverAMR::getA_local(int I1,int I2)
{
  static constexpr int BSX = BlockType::sizeX;
  static constexpr int BSY = BlockType::sizeY;
  
  int k1 =  I1 / (BSX*BSY);
  int j1 = (I1 - k1*(BSX*BSY) ) / BSX;
  int i1 =  I1 - k1*(BSX*BSY) - j1*BSX;

  int k2 =  I2 / (BSX*BSY);
  int j2 = (I2 - k2*(BSX*BSY) ) / BSX;
  int i2 =  I2 - k2*(BSX*BSY) - j2*BSX;

  if (i1==i2 && j1==j2 && k1==k2)
  {
    return 6.0;
  }
  else if (abs(i1-i2) + abs(j1-j2) + abs(k1-k2) == 1)
  {
    return -1.0;
  }
  else
  {
    return 0.0;
  }
}

PoissonSolverAMR::PoissonSolverAMR(SimulationData& s): sim(s),findLHS(s)
{
  iter_min = 20;
  if (StreamerDiv::channels != 1) {
    fprintf(stderr, "PoissonSolverScalar_MPI(): Error: StreamerDiv::channels is %d (should be 1)\n",
            StreamerDiv::channels);
    fflush(0); exit(1);
  }

  std::vector<std::vector<double>> L;

  int BSX = BlockType::sizeX;
  int BSY = BlockType::sizeY;
  int BSZ = BlockType::sizeZ;
  int N = BSX*BSY*BSZ;

  compute_ix.resize(omp_get_max_threads());
  compute_iy.resize(omp_get_max_threads());
  compute_iz.resize(omp_get_max_threads());
  #pragma omp parallel
  {
      int tid = omp_get_thread_num();
      compute_ix[tid].resize(N,0);
      compute_iy[tid].resize(N,0);
      compute_iz[tid].resize(N,0);
      for (int I = 0; I < N ; I++)
      {
        int iz =  I / (BSX*BSY);
        int iy = (I - iz*(BSX*BSY) )/ BSX;
        int ix =  I - iz*(BSX*BSY) - iy * BSX;
        compute_ix[tid][I] = ix;
        compute_iy[tid][I] = iy;
        compute_iz[tid][I] = iz;
      }
  }

  L.resize(N);

  for (int i = 0 ; i<N ; i++)
  {
    L[i].resize(i+1);
  }
  for (int i = 0 ; i<N ; i++)
  {
    double s1=0;
    for (int k=0; k<=i-1; k++)
      s1 += L[i][k]*L[i][k];
    L[i][i] = sqrt(getA_local(i,i) - s1);
    for (int j=i+1; j<N; j++)
    {
      double s2 = 0;
      for (int k=0; k<=i-1; k++)
        s2 += L[i][k]*L[j][k];
      L[j][i] = (getA_local(j,i)-s2) / L[i][i];
    }
  }
  L_row.resize(omp_get_max_threads());
  L_col.resize(omp_get_max_threads());
  Ld.resize(omp_get_max_threads());
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    L_row[tid].resize(N);
    L_col[tid].resize(N);

    for (int i = 0 ; i<N ; i++)
    {
      Ld[tid].push_back(1.0/L[i][i]);
      for (int j = 0 ; j < i ; j++)
      {
        if ( abs(L[i][j]) > 1e-10 ) L_row[tid][i].push_back({j,L[i][j]});
      }
    }    
    for (int j = 0 ; j<N ; j++)
    {
      for (int i = j ; i < N ; i++)
      {
        if ( abs(L[i][j]) > 1e-10 ) L_col[tid][j].push_back({i,L[i][j]});
      }
    }
  }
}
#else
PoissonSolverAMR::PoissonSolverAMR(SimulationData& s): sim(s),findLHS(s)
{
  iter_min = 20;
  if (StreamerDiv::channels != 1) {
    fprintf(stderr, "PoissonSolverScalar_MPI(): Error: StreamerDiv::channels is %d (should be 1)\n",
            StreamerDiv::channels);
    fflush(0); exit(1);
  }
}
#endif


void PoissonSolverAMR::getZ()
{
  static constexpr int N = BlockType::sizeX * BlockType::sizeY * BlockType::sizeZ;
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  const size_t Nblocks = vInfo.size();
  #pragma omp parallel
  { 
    const int tid = omp_get_thread_num();
    #pragma omp for schedule(runtime)
    for (size_t i=0; i < Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      const double invh = 1.0/vInfo[i].h_gridpoint;
  
      //1. z = L^{-1}r
      for (int I = 0; I < N ; I++)
      {
        double rhs = 0.0;
        for (size_t jj = 0 ; jj < L_row[tid][I].size(); jj++)
        {
          const int J = L_row[tid][I][jj].first;
          const double LIJ = L_row[tid][I][jj].second;
          const int iz = compute_iz[tid][J];
          const int iy = compute_iy[tid][J];
          const int ix = compute_ix[tid][J];
          rhs += LIJ*b(ix,iy,iz).zVector;
        }    
        const int iz = compute_iz[tid][I];
        const int iy = compute_iy[tid][I];
        const int ix = compute_ix[tid][I];
        b(ix,iy,iz).zVector = (b(ix,iy,iz).zVector - rhs) * Ld[tid][I];
      }
  
      //2. z = -(1/h)*L^T{-1}r
      for (int I = N-1; I >= 0 ; I--)
      {
        double rhs = 0.0;
        for (size_t jj = 0 ; jj < L_col[tid][I].size(); jj++)
        {
          const int J = L_col[tid][I][jj].first;
          const double LJI = L_col[tid][I][jj].second;
          const int iz = compute_iz[tid][J];
          const int iy = compute_iy[tid][J];
          const int ix = compute_ix[tid][J];
          rhs += LJI*b(ix,iy,iz).zVector;
        }
        const int iz = compute_iz[tid][I];
        const int iy = compute_iy[tid][I];
        const int ix = compute_ix[tid][I];     
        b(ix,iy,iz).zVector = -(invh*b(ix,iy,iz).zVector + rhs) * Ld[tid][I];
      }
    }
  }
}


void PoissonSolverAMR::solve()
{
    sim.startProfiler("Poisson solve");

    if (iter_min >= 1) iter_min --;

    std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
    const double eps = 1e-21;
    const size_t Nblocks = vInfo.size();
    const size_t N = BlockType::sizeX*BlockType::sizeY*BlockType::sizeZ*Nblocks;
    std::vector <Real> storeGridElements (8*N);
    #pragma omp parallel for schedule(runtime)
    for (size_t i=0; i < Nblocks; i++)
    {
        BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
        const size_t offset = _offset( vInfo[i] );
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            const size_t src_index = _dest(offset, iz, iy, ix);
            storeGridElements[8*src_index  ] = b(ix,iy,iz).chi;
            storeGridElements[8*src_index+1] = b(ix,iy,iz).u;
            storeGridElements[8*src_index+2] = b(ix,iy,iz).v;
            storeGridElements[8*src_index+3] = b(ix,iy,iz).w;
            storeGridElements[8*src_index+4] = b(ix,iy,iz).p;
            storeGridElements[8*src_index+5] = b(ix,iy,iz).tmpU;
            storeGridElements[8*src_index+6] = b(ix,iy,iz).tmpV;
            storeGridElements[8*src_index+7] = b(ix,iy,iz).tmpW;
            b(ix,iy,iz).zVector = b(ix,iy,iz).xVector;//this is done because Get_LHS works with zVector
        }
    }

    //1. rVector <-- b - AxVector
    //2. rhatVector = rVector
    //3. rho = a = omega = 1.0
    //4. vVector = pVector = 0

    //skip LHS for second order time, as x0 = 0 is the initial guess (x := P^{n+1}-P^{n} for 2nd order time)
    double norm = 0.0; 
    if (sim.TimeOrder == 1 || sim.step < sim.step_2nd_start)
    {
        findLHS(0); // AxVector <-- A*x_{0}, x_0 = pressure
        #pragma omp parallel for reduction (+:norm)
        for(size_t i=0; i< Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                b(ix,iy,iz).rVector = b.tmp[iz][iy][ix] - b(ix,iy,iz).AxVector;
                b(ix,iy,iz).rhatVector = b(ix,iy,iz).rVector;
                b(ix,iy,iz).pVector = 0.0;
                b(ix,iy,iz).vVector = 0.0;
                norm+= b(ix,iy,iz).rVector*b(ix,iy,iz).rVector;
            }
        }
    }
    else
    {
        #pragma omp parallel for reduction (+:norm)
        for(size_t i=0; i< Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                b(ix,iy,iz).rVector = b.tmp[iz][iy][ix];
                b(ix,iy,iz).rhatVector = b(ix,iy,iz).rVector;
                b(ix,iy,iz).pVector = 0.0;
                b(ix,iy,iz).vVector = 0.0;
                norm+= b(ix,iy,iz).rVector*b(ix,iy,iz).rVector;
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_DOUBLE,MPI_SUM,m_comm);
    norm = std::sqrt(norm);

    double rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    double rho_m1;

    std::vector <Real> x_opt (N);
    bool useXopt = false;
    double min_norm = 1e50;
    double init_norm=norm;
    const double max_error     = sim.step < 100 ? 0.0 : sim.PoissonErrorTol;
    const double max_rel_error = sim.step < 100 ? 0.0 : sim.PoissonErrorTolRel;
    bool serious_breakdown = false;
    int iter_opt = 0;

    //5. for k = 1,2,...
    //bool skip_loop = (norm < max_error) && (0 == iter_min); 
    //if (skip_loop)
    //{
    //  if (m_rank==0)
    //    std::cout <<  "Poisson solver converged after " <<  0 << " iterations. Error norm = " << norm << std::endl;
    //}
    //else
    for (size_t k = 1; k < 300; k++)
    {
        //1. rho_i = (rhat_0,r_{k-1})
        //2. beta = rho_{i}/rho_{i-1} * alpha/omega_{i-1}
        rho_m1 = rho;
        rho = 0.0;
        double norm_1 = 0.0;
        double norm_2 = 0.0;
        #pragma omp parallel for reduction(+:rho,norm_1,norm_2)
        for(size_t i=0; i< Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                rho += b(ix,iy,iz).rVector * b(ix,iy,iz).rhatVector;
                norm_1 += b(ix,iy,iz).rVector*b(ix,iy,iz).rVector;
                norm_2 += b(ix,iy,iz).rhatVector*b(ix,iy,iz).rhatVector;
            }
        }
        double aux_norm [3] = {rho,norm_1,norm_2};
        MPI_Allreduce(MPI_IN_PLACE,&aux_norm,3,MPI_DOUBLE,MPI_SUM,m_comm);
        rho = aux_norm[0];
        norm_1 = aux_norm[1];
        norm_2 = aux_norm[2];
        double beta = rho / (rho_m1+eps) * alpha / (omega+eps) ;

        norm_1 = sqrt(norm_1);
        norm_2 = sqrt(norm_2);
        double cosTheta = rho/norm_1/norm_2; 
        serious_breakdown = std::fabs(cosTheta) < 1e-10;
        if (serious_breakdown)
        {
            beta = 0.0;
            rho = 0.0;
            #pragma omp parallel for reduction(+:rho)
            for(size_t i=0; i< Nblocks; i++)
            {
                BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
                for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    b(ix,iy,iz).rhatVector = b(ix,iy,iz).rVector;
                    rho += b(ix,iy,iz).rVector * b(ix,iy,iz).rhatVector;
                }
            }
            if (m_rank == 0) 
                std::cout << "  [Poisson solver]: restart at iteration:" << k << 
                             "  norm:"<< norm <<" init_norm:" << init_norm << std::endl;
            MPI_Allreduce(MPI_IN_PLACE,&rho,1,MPI_DOUBLE,MPI_SUM,m_comm);
        }

        //3. p_i = r_{i-1} + beta*(p_{i-1}-omega *v_{i-1})
        //4. z = K_2^{-1} p
        #pragma omp parallel for schedule(runtime)
        for (size_t i=0; i < Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                b(ix,iy,iz).pVector = b(ix,iy,iz).rVector + beta* ( b(ix,iy,iz).pVector - omega * b(ix,iy,iz).vVector);
                b(ix,iy,iz).zVector = b(ix,iy,iz).pVector;
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
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                b(ix,iy,iz).vVector = b(ix,iy,iz).AxVector;
                alpha += b(ix,iy,iz).rhatVector * b(ix,iy,iz).vVector;
            }
        }  
        MPI_Allreduce(MPI_IN_PLACE,&alpha,1,MPI_DOUBLE,MPI_SUM,m_comm);
        alpha = rho / (alpha + eps);

        //7. x += a z
        //8. 
        //9. s = r_{i-1}-alpha * v_i
        //10. z = K_2^{-1} s
        #pragma omp parallel for schedule(runtime)
        for (size_t i=0; i < Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                b(ix,iy,iz).xVector += alpha * b(ix,iy,iz).zVector;
                b(ix,iy,iz).sVector = b(ix,iy,iz).rVector - alpha * b(ix,iy,iz).vVector;
                b(ix,iy,iz).zVector = b(ix,iy,iz).sVector;
            }
        }
        getZ();

        //11. t = Az
        findLHS(0); // t stored in AxVector

        //12. omega = ...
        double aux1 = 0;
        double aux2 = 0;
        #pragma omp parallel for reduction (+:aux1,aux2)
        for (size_t i=0; i < Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                aux1 += b(ix,iy,iz).AxVector * b(ix,iy,iz).sVector;
                aux2 += b(ix,iy,iz).AxVector * b(ix,iy,iz).AxVector;
            }
        }
        double aux_12[2] = {aux1,aux2};
        MPI_Allreduce(MPI_IN_PLACE,&aux_12,2,MPI_DOUBLE,MPI_SUM,m_comm);
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
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                b(ix,iy,iz).xVector += omega * b(ix,iy,iz).zVector;
                b(ix,iy,iz).rVector  = b(ix,iy,iz).sVector- omega * b(ix,iy,iz).AxVector;
                norm+= b(ix,iy,iz).rVector*b(ix,iy,iz).rVector; 
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_DOUBLE,MPI_SUM,m_comm);
        norm = std::sqrt(norm);

        if (norm < min_norm)
        {
            useXopt = true;
            iter_opt = k;
            min_norm = norm;
            #pragma omp parallel for schedule(runtime)
            for (size_t i=0; i < Nblocks; i++)
            {
                BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
                const size_t offset = _offset( vInfo[i] );
                for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = _dest(offset, iz, iy, ix);
                    x_opt[src_index] = b(ix,iy,iz).xVector;
                }
            }
        }
    
        if (k==1) init_norm = norm;

        if (norm / (init_norm+eps) > 2.0 && k > 10)
        {
            useXopt = true;
            if (m_rank==0)
                std::cout <<  "XOPT Poisson solver converged after " <<  k << " iterations. Error norm = " << norm << "  iter_opt="<< iter_opt << std::endl;
            break;
        }

        if ( (norm < max_error || norm/init_norm < max_rel_error ) && k > iter_min )
        {
            if (m_rank==0)
                std::cout <<  "Poisson solver converged after " <<  k << " iterations. Error norm = " << norm << std::endl;
            break;
        }

    }//k-loop

    if (useXopt)
    {
        #pragma omp parallel for schedule(runtime)
        for (size_t i=0; i < Nblocks; i++)
        {
            BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
            const size_t offset = _offset( vInfo[i] );
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const size_t src_index = _dest(offset, iz, iy, ix);
                b(ix,iy,iz).xVector = x_opt[src_index];
            }
        }
    }
  
    // Subtract average pressure from all gridpoints and copy back stuff
    const Real avgP = computeAverage();
    #pragma omp parallel for schedule(runtime)
    for (size_t i=0; i < Nblocks; i++)
    {
        BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
        const size_t offset = _offset( vInfo[i] );
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            b(ix,iy,iz).p -= avgP;//p is already the xVector!!
            const size_t src_index = _dest(offset, iz, iy, ix);     
            b(ix,iy,iz).chi  = storeGridElements[8*src_index  ];
            b(ix,iy,iz).u    = storeGridElements[8*src_index+1];
            b(ix,iy,iz).v    = storeGridElements[8*src_index+2];
            b(ix,iy,iz).w    = storeGridElements[8*src_index+3];
            //b(ix,iy,iz).p    = storeGridElements[8*src_index+4];
            b(ix,iy,iz).tmpU = storeGridElements[8*src_index+5];
            b(ix,iy,iz).tmpV = storeGridElements[8*src_index+6];
            b(ix,iy,iz).tmpW = storeGridElements[8*src_index+7];
        }
    } 
    sim.stopProfiler();
}

}//namespace cubismup3d
