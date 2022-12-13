//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PoissonSolverAMR.h"

namespace cubismup3d {

#if 0 //non-pipelined version

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
    x   .resize(N,0.0);
    r   .resize(N,0.0);
    p   .resize(N,0.0);
    v   .resize(N,0.0);
    s   .resize(N,0.0);
    rhat.resize(N,0.0);

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
        const ScalarBlock & __restrict__ Z = (*sim.pres)(i);
        ScalarBlock & __restrict__ LHS  = (*sim.lhs)(i);
	if (sim.bMeanConstraint == 1 || sim.bMeanConstraint > 2)
          if (vInfo_z[i].index[0] == 0 && vInfo_z[i].index[1] == 0 && vInfo_z[i].index[2] == 0) LHS(0,0,0).s = 0.0;
        for(int iz=0; iz<Nz; iz++)
        for(int iy=0; iy<Ny; iy++)
        for(int ix=0; ix<Nx; ix++)
        {
            const size_t src_index = _dest(vInfo_lhs[i], iz, iy, ix);
            p[src_index] = 0.0;
            v[src_index] = 0.0;
            s[src_index] = 0.0;
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
        const ScalarBlock & __restrict__ LHS  = (*sim.lhs)(i);
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
            ScalarBlock & __restrict__ Z = (*sim.pres)(i);
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
            const ScalarBlock & __restrict__ LHS  = (*sim.lhs)(i);
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
            ScalarBlock & __restrict__ Z = (*sim.pres)(i);
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
            const ScalarBlock & __restrict__ LHS  = (*sim.lhs)(i);
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
            const ScalarBlock & __restrict__ LHS = (*sim.lhs)(i);
            const ScalarBlock & __restrict__ Z   = (*sim.pres)(i);
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
        ScalarBlock & __restrict__ b  = (*sim.pres)(i);
        for(int iz=0; iz<Nz; iz++)
        for(int iy=0; iy<Ny; iy++)
        for(int ix=0; ix<Nx; ix++)
        {
            const size_t src_index = _dest(vInfo_z[i], iz, iy, ix);
            b(ix,iy,iz).s = xsol[src_index];
        }
    }
  
}
#else

void PoissonSolverAMR::solve()
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
  const Real max_error     = sim.PoissonErrorTol;          //error tolerance for Linf norm of residual
  const Real max_rel_error = sim.PoissonErrorTolRel;       //relative error tolerance for Linf(r)/Linf(r0)
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

    if (sim.bMeanConstraint == 1 || sim.bMeanConstraint > 2)
        if (AxInfo[i].index[0] == 0 && AxInfo[i].index[1] == 0 && AxInfo[i].index[2] == 0)
            rhs(0,0,0).s = 0.0;

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
    if (verbose) std::cout << "[Poisson solver]: initial error norm:" << norm << "\n";
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
    norm = std::sqrt(norm);

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
    //Check if restart should be made. If so, current solution estimate is used as an initial
    //guess and solver starts again.
    serious_breakdown = r0r * r0r < 1e-16 * norm_1 * norm_2;
    if (serious_breakdown)// && restarts < max_restarts)
    {
      //restarts ++;
      if (verbose)
        std::cout << "  [Poisson solver]: Restart at iteration: " << k << " norm: " << norm << std::endl;

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
    if ( norm < max_error || norm/(init_norm+eps) < max_rel_error )
    {
      if (verbose)
        std::cout << "  [Poisson solver]: Converged after " << k << " iterations.\n";
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
#endif

}//namespace cubismup3d
