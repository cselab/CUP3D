//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PoissonSolverAMR.h"

namespace cubismup3d {

void PoissonSolverAMR::_cub2fftw() const
{
  assert(data.size()>0);
  #ifndef NDEBUG
  {
    computeRelativeCorrection();
  }
  #endif
}

Real PoissonSolverAMR::computeRelativeCorrection() const
{
  Real sumRHS = 0, sumABS = 0;
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static) reduction(+ : sumRHS, sumABS)
  for(size_t i=0; i<vInfo.size(); ++i)
  {
    const double h3 = vInfo[i].h_gridpoint * vInfo[i].h_gridpoint * vInfo[i].h_gridpoint;
    const size_t offset = _offset( vInfo[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      sumABS += h3 * std::fabs( data[src_index] );
      sumRHS += h3 * data[src_index];
    }
  }
  double sums[2] = {sumRHS, sumABS};
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE,MPI_SUM, m_comm);
  sums[1] = std::max(std::numeric_limits<double>::epsilon(), sums[1]);
  const Real correction = sums[0] / sums[1];
  if(sim.verbose)
    printf("Relative RHS correction:%e / %e\n", sums[0], sums[1]);
  return correction;
}

Real PoissonSolverAMR::computeAverage() const
{
  Real avgP = 0;
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static) reduction(+ : avgP)
  for(size_t i=0; i<vInfo.size(); ++i)
  {
    const double h3 = vInfo[i].h_gridpoint * vInfo[i].h_gridpoint * vInfo[i].h_gridpoint;
    const size_t offset = _offset( vInfo[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      avgP += h3 * data[src_index];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPIREAL, MPI_SUM, m_comm);
  avgP /= sim.extent[0] * sim.extent[1] * sim.extent[2];
  return avgP;
}

void PoissonSolverAMR::_fftw2cub() const
{
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<vInfo.size(); ++i) {
    BlockType& b = *(BlockType*) vInfo[i].ptrBlock;
    const size_t offset = _offset( vInfo[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      b(ix,iy,iz).p = data[src_index];
    }
  }
}



struct KernelLHS
{
  const SimulationData & sim;
  bool useX;
  KernelLHS(const SimulationData&s, bool _useX) : sim(s), useX(_useX) {}

  const StencilInfo stencil{-1,-1,-1,2,2,2,
                            false, {FE_CHI, FE_U, FE_V, FE_W, FE_P, FE_TMPU, FE_TMPV, FE_TMPW} };
  void operator()(LabMPI & lab, const BlockInfo& info, FluidBlock& o) const
  {
    const double h = info.h_gridpoint; 
    if (useX)
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        o(ix,iy,iz).AxVector = h*( lab(ix-1,iy,iz).xVector + 
                                   lab(ix+1,iy,iz).xVector + 
                                   lab(ix,iy-1,iz).xVector + 
                                   lab(ix,iy+1,iz).xVector +
                                   lab(ix,iy,iz-1).xVector +
                                   lab(ix,iy,iz+1).xVector - 6.0*lab(ix,iy,iz).xVector);
      }
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        o(ix,iy,iz).AxVector = h*( lab(ix-1,iy,iz).pVector + 
                                   lab(ix+1,iy,iz).pVector + 
                                   lab(ix,iy-1,iz).pVector + 
                                   lab(ix,iy+1,iz).pVector +
                                   lab(ix,iy,iz-1).pVector +
                                   lab(ix,iy,iz+1).pVector - 6.0*lab(ix,iy,iz).pVector);
      }
    }
  }
};

void PoissonSolverAMR::Get_LHS (bool useX)
{
    //compute A*x and store it into LHS
    //here A corresponds to the discrete Laplacian operator for an AMR mesh

    //copy - pasted from Operator.h
    // (good enough for now)
    const KernelLHS kernel(sim,useX);

    cubism::SynchronizerMPI_AMR<Real,FluidGridMPI>& Synch = *sim.grid->sync(kernel);

    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i) {
      labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
      labs[i].prepare(* sim.grid, Synch);
    }

    int rank;
    MPI_Comm_rank(sim.grid->getCartComm(), &rank);
    MPI_Barrier(sim.grid->getCartComm());
    std::vector<cubism::BlockInfo*> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(static)
      for(int i=0; i<Ninner; i++) {
        const cubism::BlockInfo I = *avail0[i];
        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
      }
    }

    if(sim.nprocs>1)
    {
      std::vector<cubism::BlockInfo*> avail1 = Synch.avail_halo();
      const int Nhalo = avail1.size();

      #pragma omp parallel
      {
        int tid = omp_get_thread_num();
        LabMPI& lab = labs[tid];

        #pragma omp for schedule(static)
        for(int i=0; i<Nhalo; i++) {
          const cubism::BlockInfo I = *avail1[i];
          FluidBlock& b = *(FluidBlock*)I.ptrBlock;
          lab.load(I, 0);
          kernel(lab, I, b);
        }
      }
    }

    if(labs != nullptr) {
      delete [] labs;
      labs = nullptr;
    }

    MPI_Barrier(sim.grid->getCartComm());
}

void PoissonSolverAMR::solve()
{
  sim.startProfiler("PoissonSolverAMR cub2rhs");
  _cub2fftw();
  sim.stopProfiler();

  static constexpr int BSX = BlockType::sizeX;
  static constexpr int BSY = BlockType::sizeY;
  static constexpr int BSZ = BlockType::sizeZ;

  std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();

  const size_t Nblocks = vInfo.size();
  const size_t N = BSX*BSY*BSZ*Nblocks;
  size_t Nsystem;
  MPI_Allreduce(&N,&Nsystem,1,MPI_LONG,MPI_SUM,m_comm);

  assert (N == data.size());

  /*
  The following vectors are needed:

  p: conjugate vector
  r: residual vector
  x: current solution estimate
  s: best solution estimate
  data: rhs of system, or the current Ax (LHS * x)

  They are stored in tmpU,tmpV,tmpW, u and v. 

  */

  std::vector <Real> storeGridElements (6*N);

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
      storeGridElements[6*src_index  ] = b(ix,iy,iz).pVector ;
      storeGridElements[6*src_index+1] = b(ix,iy,iz).rVector ;
      storeGridElements[6*src_index+2] = b(ix,iy,iz).xVector ;
      storeGridElements[6*src_index+3] = b(ix,iy,iz).sVector ;
      storeGridElements[6*src_index+4] = b(ix,iy,iz).AxVector;
      storeGridElements[6*src_index+5] = b(ix,iy,iz).zVector;
    }
  }

  const double max_error = 1e-6;
  double err_min         = 1e50;
  double err             = 1;
  double alpha           = 1.0;
  double beta            = 0.0;
  double rk_rk,rkp1_rkp1,pAp;
  #ifdef PRECOND
  double rk_zk;
  double rkp1_zkp1;
  #endif

  //Copy system RHS (stored in data) to rVector
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
      b(ix,iy,iz).rVector = data[src_index];
    }
  }

  Get_LHS(true); // AxVector <-- A*x_{0}, x_0 = pressure

  //rVector <-- rVector - alpha * AxVector
  #pragma omp parallel for schedule(runtime)
  for(size_t i=0; i< Nblocks; i++)
  {
    BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++)
      b(ix,iy,iz).rVector -= alpha * b(ix,iy,iz).AxVector;
  }

  //rk_rk = rVector /cdot rVector
  rk_rk = 0.0;
  #pragma omp parallel for reduction(+:rk_rk)
  for(size_t i=0; i< Nblocks; i++)
  {
    BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++)
      rk_rk += b(ix,iy,iz).rVector * b(ix,iy,iz).rVector;
  }
  MPI_Allreduce(MPI_IN_PLACE,&rk_rk,1,MPI_DOUBLE,MPI_SUM,m_comm);
  
  #ifdef PRECOND
    FindZ();
    //rk_zk = rVector /cdot rVector
    rk_zk = 0.0;
    #pragma omp parallel for reduction(+:rk_zk)
    for(size_t i=0; i< Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        rk_zk += b(ix,iy,iz).rVector * b(ix,iy,iz).zVector;
    }
    MPI_Allreduce(MPI_IN_PLACE,&rk_zk,1,MPI_DOUBLE,MPI_SUM,m_comm);

    //pVector = zVector
    #pragma omp parallel for schedule(runtime)
    for (size_t i=0; i < Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        b(ix,iy,iz).pVector = b(ix,iy,iz).zVector;
    }
  #else
    //pVector = rVector
    #pragma omp parallel for schedule(runtime)
    for (size_t i=0; i < Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        b(ix,iy,iz).pVector = b(ix,iy,iz).rVector;
    }
  #endif
 
  bool flag = false;
  int count = 0;


  for (size_t k = 1; k < 10000; k++)
  {    
    count++;
    err = std::sqrt(rk_rk)/Nsystem;

    if (m_rank == 0 && (k % 500 == 0|| k<5) ) std::cout << "k=" << k<<" err="<<err << std::endl;

    if (err < err_min && k >=2)
    {
      flag = true;
      err_min = err;
      #pragma omp parallel for schedule(runtime)
      for (size_t i=0; i < Nblocks; i++)
      {
        BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
          b(ix,iy,iz).sVector = b(ix,iy,iz).xVector;
      }
    }

    if (err < max_error && k > 10) break;
    if (  err/(err_min+1e-21) > 10.0 && k > 20) break; //error grows, stop iterations!

    Get_LHS(false);// AxVector <-- A*pVector

    //pAp = pVector /cdot AxVector
    pAp = 0.0;
    #pragma omp parallel for reduction(+:pAp)
    for(size_t i=0; i< Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        pAp += b(ix,iy,iz).pVector * b(ix,iy,iz).AxVector;
    }
    MPI_Allreduce(MPI_IN_PLACE,&pAp,1,MPI_DOUBLE,MPI_SUM,m_comm);

    if ( std::fabs(pAp) < 1e-21) {std::cout << "CG:pAp is small"<< std::endl; break;}

    #ifdef PRECOND
       alpha = rk_zk/pAp;
    #else
       alpha = rk_rk/pAp;
    #endif

    //x_{k+1} <-- x_{k} + alpha * p_{k}
    #pragma omp parallel for schedule(runtime)
    for(size_t i=0; i< Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        b(ix,iy,iz).xVector += alpha * b(ix,iy,iz).pVector;
    }

    //r_{k+1} <-- r_{k} - alpha * Ax
    #pragma omp parallel for schedule(runtime)
    for(size_t i=0; i< Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        b(ix,iy,iz).rVector -= alpha * b(ix,iy,iz).AxVector;
    }

    //rkp1_rkp1 = rVector /cdot rVector
    rkp1_rkp1 = 0.0;
    #pragma omp parallel for reduction(+:rkp1_rkp1)
    for(size_t i=0; i< Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        rkp1_rkp1 += b(ix,iy,iz).rVector * b(ix,iy,iz).rVector;
    }
    MPI_Allreduce(MPI_IN_PLACE,&rkp1_rkp1,1,MPI_DOUBLE,MPI_SUM,m_comm);


    #ifdef PRECOND
      FindZ();
      rkp1_zkp1 = 0.0;
      #pragma omp parallel for reduction(+:rkp1_zkp1)
      for(size_t i=0; i< Nblocks; i++)
      {
        BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
          rkp1_zkp1 += b(ix,iy,iz).rVector * b(ix,iy,iz).zVector;
      }
      MPI_Allreduce(MPI_IN_PLACE,&rkp1_zkp1,1,MPI_DOUBLE,MPI_SUM,m_comm);
      beta = rkp1_zkp1 / rk_zk;
      rk_zk = rkp1_zkp1;  
      
      //p_{k+1} <-- z_{k+1} + beta * p_{k}   
      #pragma omp parallel for schedule(runtime)
      for(size_t i=0; i< Nblocks; i++)
      {
        BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
          b(ix,iy,iz).pVector = beta * b(ix,iy,iz).pVector + b(ix,iy,iz).zVector;
      }

    #else
      beta = rkp1_rkp1 / rk_rk;
      //p_{k+1} <-- r_{k+1} + beta * p_{k}    
      #pragma omp parallel for schedule(runtime)
      for(size_t i=0; i< Nblocks; i++)
      {
        BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
          b(ix,iy,iz).pVector = beta * b(ix,iy,iz).pVector + b(ix,iy,iz).rVector;
      }
    #endif

    rk_rk = rkp1_rkp1;
  }

  if (flag)
  {
    #pragma omp parallel for schedule(runtime)
    for(size_t i=0; i< Nblocks; i++)
    {
      BlockType & __restrict__ b  = *(BlockType*) vInfo[i].ptrBlock;
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
        b(ix,iy,iz).xVector = b(ix,iy,iz).sVector;
    }
  }
  else
  {
    err_min = err;
  }

  std::cout << "CG Poisson solver took "<<count << " iterations. Final residual norm = "<< err_min << std::endl;


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
      data[src_index] = b(ix,iy,iz).xVector;
    }
  }
  const Real avgP = computeAverage();
  // Subtract average pressure from all gridpoints
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < data.size(); i++) data[i] -= avgP;
  std::cout <<"average = " << avgP << std::endl;
  //NOT DONE YET!!!!
  //// // Set this new mean-0 pressure as next guess
  //// // Save pressure of a corner of the grid so that it can be imposed next time
  //// pLast = data[fixed_idx];
  //// if(sim.verbose) printf("Avg Pressure:%f\n", avgP);


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
      b(ix,iy,iz).pVector  = storeGridElements[6*src_index  ];
      b(ix,iy,iz).rVector  = storeGridElements[6*src_index+1];
      b(ix,iy,iz).xVector  = storeGridElements[6*src_index+2];
      b(ix,iy,iz).sVector  = storeGridElements[6*src_index+3];
      b(ix,iy,iz).AxVector = storeGridElements[6*src_index+4];
      b(ix,iy,iz).zVector  = storeGridElements[6*src_index+5];
    }
  }
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

PoissonSolverAMR::PoissonSolverAMR(SimulationData& s): sim(s)
{
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

void PoissonSolverAMR::FindZ()
{
  static constexpr int BSX = BlockType::sizeX;
  static constexpr int BSY = BlockType::sizeY;
  static constexpr int BSZ = BlockType::sizeZ;
  static constexpr int N   = BSX*BSY*BSZ;

  std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  const size_t Nblocks = vInfo.size();

  #pragma omp parallel
  { 
    int tid = omp_get_thread_num();
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
          int J = L_row[tid][I][jj].first;
          double LIJ = L_row[tid][I][jj].second;
          int iz =  J / (BSX*BSY);
          int iy = (J - iz*(BSX*BSY) )/ BSX;
          int ix =  J - iz*(BSX*BSY) - iy * BSX;
          rhs += LIJ*b(ix,iy,iz).zVector;
        }    
        int iz =  I / (BSX*BSY);
        int iy = (I - iz*(BSX*BSY) )/ BSX;
        int ix =  I - iz*(BSX*BSY) - iy * BSX;
        b(ix,iy,iz).zVector = (b(ix,iy,iz).rVector - rhs) *Ld[tid][I];
      }
  

      //2. z = L^T{-1}r
      for (int I = N-1; I >= 0 ; I--)
      {
        double rhs = 0.0;
        for (size_t jj = 0 ; jj < L_col[tid][I].size(); jj++)
        {
          int J = L_col[tid][I][jj].first;
          double LJI = L_col[tid][I][jj].second;

          int iz =  J / (BSX*BSY);
          int iy = (J - iz*(BSX*BSY) )/ BSX;
          int ix =  J - iz*(BSX*BSY) - iy * BSX;
          rhs += LJI*b(ix,iy,iz).zVector;
        }
        int iz =  I / (BSX*BSY);
        int iy = (I - iz*(BSX*BSY) )/ BSX;
        int ix =  I - iz*(BSX*BSY) - iy * BSX;
        b(ix,iy,iz).zVector = (b(ix,iy,iz).zVector - rhs) *Ld[tid][I];
      }

      for (int iz=0;iz<BSZ;iz++)
      for (int iy=0;iy<BSY;iy++)
      for (int ix=0;ix<BSX;ix++)
        b(ix,iy,iz).zVector = -invh*b(ix,iy,iz).zVector;
    }
  }
}
PoissonSolverAMR::~PoissonSolverAMR()
{
  std::cout << "~PoissonSolverAMR" << std::endl;
  data.clear();
}
#else
PoissonSolverAMR::PoissonSolverAMR(SimulationData& s): sim(s)
{
  if (StreamerDiv::channels != 1) {
    fprintf(stderr, "PoissonSolverScalar_MPI(): Error: StreamerDiv::channels is %d (should be 1)\n",
            StreamerDiv::channels);
    fflush(0); exit(1);
  }
}

PoissonSolverAMR::~PoissonSolverAMR()
{
  std::cout << "~PoissonSolverAMR" << std::endl;
  data.clear();
}
#endif

}//namespace cubismup3d