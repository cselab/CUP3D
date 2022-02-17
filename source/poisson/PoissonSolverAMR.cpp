//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PoissonSolverAMR.h"

namespace cubismup3d {

PoissonSolverAMR::PoissonSolverAMR(SimulationData& s): sim(s),findLHS(s)
{
    iter_min = 20;
}

void PoissonSolverAMR::getZ()
{
  const std::vector<cubism::BlockInfo>& vInfo = gridPoisson.getBlocksInfo();
  const size_t Nblocks = vInfo.size();

  #pragma omp parallel
  {
    const int nx = BlockType::sizeX;
    const int ny = BlockType::sizeY;
    const int nz = BlockType::sizeZ;
    const int nx2 = nx + 2;
    const int ny2 = ny + 2;
    const int nz2 = nz + 2;
    const int N = nx*ny*nz;
    const int N2 = nx2*ny2*nz2;
    std::vector<double> p (N2,0.0);
    std::vector<double> r (N ,0.0);
    std::vector<double> x (N ,0.0);
    std::vector<double> Ax(N ,0.0);

    #pragma omp for
    for (size_t i=0; i < Nblocks; i++)
    {
        BlockTypePoisson & __restrict__ b  = *(BlockTypePoisson*) vInfo[i].ptrBlock;
        const double invh = 1.0/vInfo[i].h;
        double norm0 = 0;
        double rr = 0;
        double a2 = 0;
        double beta = 0;
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            x[ix+iy*nx+iz*nx*ny] = 0.0;
            r[ix+iy*nx+iz*nx*ny] = invh*b(ix,iy,iz).s;
            norm0 += r[ix+iy*nx+iz*nx*ny]*r[ix+iy*nx+iz*nx*ny];
            p[(ix+1)+(iy+1)*(nx+2)+(iz+1)*(nx+2)*(ny+2)] = r[ix+iy*nx+iz*nx*ny];
            rr += r[ix+iy*nx+iz*nx*ny]*r[ix+iy*nx+iz*nx*ny];
        }
        norm0 = sqrt(norm0)/N;
        double norm = 0;
        if (norm0 > 1e-16)
        for (int k = 0 ; k < 100 ; k ++)
        {
            a2 = 0;
            for(int iz=0; iz<nz; iz++)
            for(int iy=0; iy<ny; iy++)
            for(int ix=0; ix<nx; ix++)
            {
                const int index1 = ix+iy*nx+iz*nx*ny;
                const int index2 = (ix+1)+(iy+1)*(nx+2)+(iz+1)*(nx+2)*(ny+2);
                Ax[index1] = p[index2 + 1] + p[index2 - 1] + p[index2 + nx2] + p[index2 - nx2] + p[index2 + nx2*ny2] + p[index2 - nx2*ny2] - 6.0*p[index2];
                a2 += p[index2]*Ax[index1];
            }
            const double a = rr/(a2+1e-55);
            double norm_new = 0;
            beta = 0;
            for(int iz=0; iz<nz; iz++)
            for(int iy=0; iy<ny; iy++)
            for(int ix=0; ix<nx; ix++)
            {
                const int index1 = ix+iy*nx+iz*nx*ny;
                const int index2 = (ix+1)+(iy+1)*(nx+2)+(iz+1)*(nx+2)*(ny+2);
                x[index1] += a*p [index2];
                r[index1] -= a*Ax[index1];
                norm_new += r[index1]*r[index1];
            }
            beta = norm_new;
            norm_new = sqrt(norm_new)/N;
            //if (k>max(max(nx,ny),nz) && std::fabs(norm/norm_new - 1.0) < 1e-7) break;
            norm = norm_new;
            if (norm/norm0< 1e-7) break;
            double temp = rr;
            rr = beta;
            beta /= (temp+1e-55);
            for(int iz=0; iz<nz; iz++)
            for(int iy=0; iy<ny; iy++)
            for(int ix=0; ix<nx; ix++)
            {
                const int index1 = ix+iy*nx+iz*nx*ny;
                const int index2 = (ix+1)+(iy+1)*(nx+2)+(iz+1)*(nx+2)*(ny+2);
                p[index2] =r[index1] + beta*p[index2];
            }
        }
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            const int index1 = ix+iy*nx+iz*nx*ny;
            b(ix,iy,iz).s = x[index1];
        }
    }    
  }
}

void PoissonSolverAMR::solve()
{

    if (iter_min >= 1) iter_min --;

    std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
    const double eps = 1e-21;
    const size_t Nblocks = vInfo.size();
    const size_t N = BlockType::sizeX*BlockType::sizeY*BlockType::sizeZ*Nblocks;
    std::vector<double> x   (N,0.0);
    std::vector<double> r   (N,0.0);
    std::vector<double> p   (N,0.0);
    std::vector<double> v   (N,0.0);
    std::vector<double> s   (N,0.0);
    std::vector<double> rhat(N,0.0);
    std::vector<cubism::BlockInfo>& vInfoPoisson = gridPoisson.getBlocksInfo();

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
        const int m = vInfoPoisson[i].level;
        const long long n = vInfoPoisson[i].Z;
        const BlockInfo & info = grid.getBlockInfoAll(m,n);
        BlockType & __restrict__ b  = *(BlockType*) info.ptrBlock;
        BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
        const size_t offset = _offset( vInfoPoisson[i] );
        if (vInfoPoisson[i].index[0] == 0 && 
            vInfoPoisson[i].index[1] == 0 && 
            vInfoPoisson[i].index[2] == 0)
            b.tmp[4][4][4] = 0;

        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            const size_t src_index = _dest(offset, iz, iy, ix);
            x[src_index]         = b(ix,iy,iz).p;
            bPoisson(ix,iy,iz).s = b(ix,iy,iz).p; //this is done because Get_LHS works with zVector
        }
    }


    //1. rVector <-- b - AxVector
    //2. rhatVector = rVector
    //3. rho = a = omega = 1.0
    //4. vVector = pVector = 0

    double norm = 0.0; 
    findLHS(0); // AxVector <-- A*x_{0}, x_0 = pressure
    #pragma omp parallel for reduction (+:norm)
    for(size_t i=0; i< Nblocks; i++)
    {
        const int m = vInfoPoisson[i].level;
        const long long n = vInfoPoisson[i].Z;
        const BlockInfo & info = grid.getBlockInfoAll(m,n);
        BlockType & __restrict__ b  = *(BlockType*) info.ptrBlock;
        BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
        const size_t offset = _offset( vInfoPoisson[i] );
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            const size_t src_index = _dest(offset, iz, iy, ix);
            r[src_index] = b.tmp[iz][iy][ix] - bPoisson(ix,iy,iz).lhs;
            rhat[src_index] = r[src_index];
            norm+= r[src_index]*r[src_index];
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
    const double max_error     = sim.PoissonErrorTol;
    const double max_rel_error = sim.PoissonErrorTolRel;
    bool serious_breakdown = false;
    int iter_opt = 0;

    //5. for k = 1,2,...
    for (size_t k = 1; k < 1000; k++)
    {
        //1. rho_i = (rhat_0,r_{k-1})
        //2. beta = rho_{i}/rho_{i-1} * alpha/omega_{i-1}
        rho_m1 = rho;
        rho = 0.0;
        double norm_1 = 0.0;
        double norm_2 = 0.0;
        #pragma omp parallel for reduction(+:rho,norm_1,norm_2)
        for(size_t i=0; i< N; i++)
        {
            rho += r[i]*rhat[i];
            norm_1 += r[i]*r[i];
            norm_2 += rhat[i]*rhat[i];
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
            if (m_rank == 0) 
                std::cout << "  [Poisson solver]: restart at iteration:" << k << 
                             "  norm:"<< norm <<" init_norm:" << init_norm << std::endl;
            MPI_Allreduce(MPI_IN_PLACE,&rho,1,MPI_DOUBLE,MPI_SUM,m_comm);
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
            BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
            const size_t offset = _offset( vInfoPoisson[i] );
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const size_t src_index = _dest(offset, iz, iy, ix);
                p[src_index] = r[src_index] + beta*(p[src_index]-omega*v[src_index]);
                bPoisson(ix,iy,iz).s = p[src_index];
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
            BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
            const size_t offset = _offset( vInfoPoisson[i] );
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const size_t src_index = _dest(offset, iz, iy, ix);
                v[src_index] = bPoisson(ix,iy,iz).lhs;//bPoisson.LHS[iz][iy][ix];
                alpha += rhat[src_index] * v[src_index];
            }
        }  
        MPI_Allreduce(MPI_IN_PLACE,&alpha,1,MPI_DOUBLE,MPI_SUM,m_comm);
        alpha = rho / (alpha + eps);
        //7. x += a z
        //8. 
        //9. s = r_{i-1}-alpha * v_i
        //10. z = K_2^{-1} s
        #pragma omp parallel for
        for (size_t i=0; i < Nblocks; i++)
        {
            BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
            const size_t offset = _offset( vInfoPoisson[i] );
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const size_t src_index = _dest(offset, iz, iy, ix);
                x[src_index] += alpha * bPoisson(ix,iy,iz).s;
                s[src_index] = r[src_index] - alpha * v[src_index];
                bPoisson(ix,iy,iz).s = s[src_index];
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
            BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
            const size_t offset = _offset( vInfoPoisson[i] );
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const size_t src_index = _dest(offset, iz, iy, ix);
                aux1 += bPoisson(ix,iy,iz).lhs * s[src_index];
                aux2 += bPoisson(ix,iy,iz).lhs * bPoisson(ix,iy,iz).lhs;
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
            BlockTypePoisson & __restrict__ bPoisson  = *(BlockTypePoisson*) vInfoPoisson[i].ptrBlock;
            const size_t offset = _offset( vInfoPoisson[i] );
            for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const size_t src_index = _dest(offset, iz, iy, ix);
                x[src_index] += omega * bPoisson(ix,iy,iz).s;
                r[src_index] = s[src_index]-omega*bPoisson(ix,iy,iz).lhs;
                norm+= r[src_index]*r[src_index];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_DOUBLE,MPI_SUM,m_comm);
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
            if (m_rank==0)
                std::cout <<  "XOPT Poisson solver converged after " <<  k << " iterations. Error norm = " << norm << "  iter_opt="<< iter_opt << std::endl;
            break;
        }

        if ( (norm < max_error || norm/init_norm < max_rel_error ) && k > iter_min )
        {
            if (m_rank==0)
                std::cout <<  "  [Poisson solver]: Converged after " <<  k << " iterations. Error norm = " << norm << std::endl;
            break;
        }

    }//k-loop

    double * xsol = useXopt ? x_opt.data():x.data();
    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
        const int m = vInfoPoisson[i].level;
        const long long n = vInfoPoisson[i].Z;
        const BlockInfo & info = grid.getBlockInfoAll(m,n);
        BlockType & __restrict__ b  = *(BlockType*) info.ptrBlock;
        const size_t offset = _offset( vInfoPoisson[i] );
        for(int iz=0; iz<BlockType::sizeZ; iz++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int ix=0; ix<BlockType::sizeX; ix++)
        {
            const size_t src_index = _dest(offset, iz, iy, ix);
            b(ix,iy,iz).p = xsol[src_index];
        }
    }
  
}
}//namespace cubismup3d
