//
//  MultigridHypre.h
//  CubismUP_3D
//
//	Operates on
//		tmp, divU
//
//  Created by Christian Conti on 1/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorMultigrid_h
#define CubismUP_3D_OperatorMultigrid_h

#include <HYPRE_struct_ls.h>
#include <vector>
#include <omp.h>
#include "Layer.h"

class MultigridHypre
{
private:
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   A;
    HYPRE_StructVector   b;
    HYPRE_StructVector   x;
    HYPRE_StructSolver   solver;
    
    bool bSetup;
    bool bConstantCoefficients;
    
    int sizeX, sizeY, size2;
    int sX, sY;
    int rankX, rankY, nranks;
    const int nprocsX, nprocsY;
    
    FluidGridMPI * fluidgrid;
    
    inline void _neighbor(Real c, Real e, Real w, Real n, Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS)
    {
        avgE = c;
        avgW = w;
        avgN = n;
        avgS = c;
    }
    
    inline void _mean(Real c, Real e, Real w, Real n, Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS)
    {
        avgE = .5 * (c + e);
        avgW = .5 * (c + w);
        avgN = .5 * (c + n);
        avgS = .5 * (c + s);
    }
    
    inline void _harmonicAvg(Real c, Real e, Real w, Real n, Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS)
    {
        avgE = 2. * c * e / (c + e);
        avgW = 2. * c * w / (c + w);
        avgN = 2. * c * n / (c + n);
        avgS = 2. * c * s / (c + s);
    }
    
    // fills with 1/rho
    template<typename Lab>
    void _fillMatrixVariableOMP(vector<double>& values)
    {
        vector<BlockInfo> myInfo = fluidgrid->getBlocksInfo();
        BlockInfo * ary = &myInfo.front();
        
        int Nx = fluidgrid->getBlocksPerDimension(0);
        int Ny = fluidgrid->getBlocksPerDimension(1);
        
        if (nranks == nprocsX*nprocsY)
        {
            Nx = fluidgrid->getBlocksPerDimension(0)/nprocsX;
            Ny = fluidgrid->getBlocksPerDimension(1)/nprocsY;
            assert(Nx*FluidBlock::sizeX == sX);
            assert(Ny*FluidBlock::sizeY == sY);
        }
        
#pragma omp parallel
        {
            Lab lab;
            const int stencil_start[3] = {-1,-1,0};
            const int stencil_end[3] = {2,2,1};
            lab.prepare(*fluidgrid, stencil_start, stencil_end, false);
            
#pragma omp for schedule(static) collapse(2)
            for(int by=0; by<Ny; by++)
                for(int bx=0; bx<Nx; bx++)
                {
                    const int i = (bx+rankX*Nx) + (by+rankY*Ny)*fluidgrid->getBlocksPerDimension(0);
                    if (i>=myInfo.size())
                        cout << bx << " " << rankX << " " << Nx << " " << by << " " << rankY << " " << Ny << " " << fluidgrid->getBlocksPerDimension(0) << " " << i << "\t" << myInfo.size() << endl;
                    assert(i<myInfo.size());
                    lab.load(ary[i], 0);
                    
                    BlockInfo info = myInfo[i];
                    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
                    
                    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                        {
                            const int gix = ix + info.index[0]*FluidBlock::sizeX;
                            const int giy = iy + info.index[1]*FluidBlock::sizeY;
                            const int lix = gix - sX*rankX;
                            const int liy = giy - sY*rankY;
                            const int lidx = (lix + liy*sX)*5;
                            
                            assert(lidx>=0);
                            assert(lidx<size2*5);
                            assert(gix>=0);
                            assert(gix<sizeX);
                            assert(giy>=0);
                            assert(giy<sizeY);
                            
                            FluidElement& phi  = lab(ix  ,iy  );
                            FluidElement& phiN = lab(ix  ,iy+1);
                            FluidElement& phiS = lab(ix  ,iy-1);
                            FluidElement& phiE = lab(ix+1,iy  );
                            FluidElement& phiW = lab(ix-1,iy  );
                            
                            Real rhoE, rhoW, rhoN, rhoS;
                            //_neighbor(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
                            _mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
                            //_harmonicAvg(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
                            
#ifdef _MIXED_
                            if (giy>0 && giy<sizeY-1)
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS);
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;
                                values[lidx+2] = 1./rhoE;
                                values[lidx+3] = 1./rhoS;
                                values[lidx+4] = 1./rhoN;
                            }
                            else if (giy==0) // neumann bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN);
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;
                                values[lidx+2] = 1./rhoE;
                                values[lidx+3] = 0;
                                values[lidx+4] = 1./rhoN;
                            }
                            else if (giy==sizeY-1) // dirichlet bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS);
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;
                                values[lidx+2] = 1./rhoE;
                                values[lidx+3] = 1./rhoS;
                                values[lidx+4] = 0;
                            }
#endif // _MIXED_
                            
#ifdef _PIPE_
                            if (giy>0 && giy<sizeY-1)
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS);
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;
                                values[lidx+2] = 1./rhoE;
                                values[lidx+3] = 1./rhoS;
                                values[lidx+4] = 1./rhoN;
                            }
                            else if (giy==0) // neumann bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN);
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;
                                values[lidx+2] = 1./rhoE;
                                values[lidx+3] = 0;
                                values[lidx+4] = 1./rhoN;
                            }
                            else if (giy==sizeY-1) // neumann bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoS);
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;
                                values[lidx+2] = 1./rhoE;
                                values[lidx+3] = 1./rhoS;
                                values[lidx+4] = 0;
                            }
#endif // _PIPE_
                            
#if 0
                            // needed for walls
                            if (giy>0 && giy<sizeY-1 && gix>0 && gix<sizeX-1)
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS);//-(rhoE + rhoW + rhoN + rhoS);//
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;//rhoW;//
                                values[lidx+2] = 1./rhoE;//rhoE;//
                                values[lidx+3] = 1./rhoS;//rhoS;//
                                values[lidx+4] = 1./rhoN;//rhoN;//
                            }
                            
                            if (giy==0) // neumann bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN);//-(rhoE + rhoW + rhoN);//
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;//rhoW;//
                                values[lidx+2] = 1./rhoE;//rhoE;//
                                values[lidx+3] = 0;
                                values[lidx+4] = 1./rhoN;//rhoN;//
                            }
                            else if (giy==sizeY-1) // dirichlet bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS);//-(rhoE + rhoW + rhoN + rhoS);//
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;//rhoW;//
                                values[lidx+2] = 1./rhoE;//rhoE;//
                                values[lidx+3] = 1./rhoS;//rhoS;//
                                values[lidx+4] = 0;
                            }
                            
                            if (gix==0) // neumann bc
                            {
                                values[lidx  ] = -(1./rhoE + 1./rhoN + 1./rhoS);//-(rhoE + rhoN + rhoS);//
                                
                                // W,E,S,N
                                values[lidx+1] = 0;
                                values[lidx+2] = 1./rhoE;//rhoE;//
                                values[lidx+3] = 1./rhoS;//rhoS;//
                                values[lidx+4] = 1./rhoN;//rhoN;//
                            }
                            else if (gix==sizeX-1) // dirichlet bc
                            {
                                values[lidx  ] = -(1./rhoW + 1./rhoN + 1./rhoS);//-(rhoE + rhoW + rhoN + rhoS);//
                                
                                // W,E,S,N
                                values[lidx+1] = 1./rhoW;//rhoW;//
                                values[lidx+2] = 0;
                                values[lidx+3] = 1./rhoS;//rhoS;//
                                values[lidx+4] = 1./rhoN;//rhoN;//
                            }
#endif
                            
#ifdef _PERIODIC_
                            values[lidx  ] = -(1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS);
                            
                            // W,E,S,N
                            values[lidx+1] = 1./rhoW;
                            values[lidx+2] = 1./rhoE;
                            values[lidx+3] = 1./rhoS;
                            values[lidx+4] = 1./rhoN;
#endif // _PERIODIC_
                        }
                }
        }
    }
    
    template<typename Lab>
    void _fillMatrixConstOMP(vector<double>& values)
    {
        /* We have size2 grid points, each with 5 stencil entries */
#pragma omp parallel for
        for (int i=0; i<values.size(); i+=5)
        {
            values[i  ] = -4;
            values[i+1] = 1;
            values[i+2] = 1;
            values[i+3] = 1;
            values[i+4] = 1;
        }
        
#ifndef _PERIODIC_
        // need to set up boundary conditions
        if (rankY==0)
#pragma omp parallel for
            for (int ix=0; ix<sX; ix++)
            {
                // bottom BC - 0-Neumann for pressure
                values[ix*5  ] = -3;
                values[ix*5+1] = 1;
                values[ix*5+2] = 1;
                values[ix*5+3] = 0;
                values[ix*5+4] = 1;
            }
        else if (rankY==nprocsY-1)
#pragma omp parallel for
            for (int ix=0; ix<sX; ix++)
            {
                // top BC - 0-Dirichlet for pressure
                values[(sY-1)*sX*5+ix*5  ] = -4;
                values[(sY-1)*sX*5+ix*5+1] = 1;
                values[(sY-1)*sX*5+ix*5+2] = 1;
                values[(sY-1)*sX*5+ix*5+3] = 1;
                values[(sY-1)*sX*5+ix*5+4] = 0;
            }
#endif // _PERIODIC_
    }
    
    void _getResultsOMP()
    {
        vector<BlockInfo> myInfo = fluidgrid->getBlocksInfo();
        BlockInfo * ary = &myInfo.front();
        const int N = myInfo.size();
        
        int ilower[2] = { 0, 0 };
        int iupper[2]= { sizeX-1, sizeY-1 };
        
        if (nranks>1)
        {
            ilower[0] = rankX * sX;
            ilower[1] = rankY * sY;
            iupper[0] = (rankX+1) * sX - 1;
            iupper[1] = (rankY+1) * sY - 1;
        }
        
        vector<double> values(size2);
        HYPRE_StructVectorGetBoxValues(x, ilower, iupper, &values[0]);
        
#ifdef _PERIODIC_
        double avg = 0;
#pragma omp parallel for schedule(static) reduction(+:avg)
        for(int i=0; i<size2; i++)
            avg += values[i];
        avg /= size2;
#endif // _PERIODIC_
        
        int Nx = fluidgrid->getBlocksPerDimension(0);
        int Ny = fluidgrid->getBlocksPerDimension(1);
        
        if (nranks == nprocsX*nprocsY)
        {
            Nx = fluidgrid->getBlocksPerDimension(0)/nprocsX;
            Ny = fluidgrid->getBlocksPerDimension(1)/nprocsY;
        }
        
#pragma omp parallel for schedule(static) collapse(2)
        for(int by=0; by<Ny; by++)
            for(int bx=0; bx<Nx; bx++)
            {
                const int i = bx+rankX*Nx + (by+rankY*Ny)*fluidgrid->getBlocksPerDimension(0);
                assert(i<N);
                
                BlockInfo info = myInfo[i];
                FluidBlock& b = *(FluidBlock*)info.ptrBlock;
                
                for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                    {
                        const int gix = ix + info.index[0]*FluidBlock::sizeX - sX*rankX;
                        const int giy = iy + info.index[1]*FluidBlock::sizeY - sY*rankY;
                        const int idx = gix + giy*sX;
                        assert(idx>=0);
                        assert(idx<size2);
                        
                        assert(!std::isnan(values[idx]));
                        
#ifndef _PERIODIC_
                        b(ix,iy).tmp  = values[idx]; // this is used for debugging
                        b(ix,iy).divU = values[idx];
#else // _PERIODIC_
                        b(ix,iy).tmp  = values[idx]-avg; // this is used for debugging
                        b(ix,iy).divU = values[idx]-avg;
#endif // _PERIODIC_
                    }
            }
        
        MPI_Request request[N];
        MPI_Status status[N];
        if (rankX!=0 || rankY!=0)
        {
            // rank != 0 - send
            for(int by=0; by<Ny; by++)
                for(int bx=0; bx<Nx; bx++)
                {
                    const int blockID = bx+rankX*Nx + (by+rankY*Ny)*fluidgrid->getBlocksPerDimension(0);
                    assert(blockID<N);
                    FluidBlock& b = (*fluidgrid)(bx+rankX*Nx,by+rankY*Ny);
                    int ierror = MPI_Isend(&b, sizeof(b), MPI_BYTE, 0, blockID, MPI_COMM_WORLD, &request[blockID]);
                }
        }
        else
        {
            // rank == 0 - receive
            for(int by=0; by<fluidgrid->getBlocksPerDimension(1); by++)
                for(int bx=0; bx<fluidgrid->getBlocksPerDimension(0); bx++)
                {
                    // skip self
                    if (bx<Nx && by<Ny)
                        continue;
                    
                    const int blockID = bx + by*fluidgrid->getBlocksPerDimension(0);
                    assert(blockID<N);
                    FluidBlock& b = (*fluidgrid)(bx,by);
                    const int senderID = bx/Nx + (by/Ny)*nprocsX;
                    assert(senderID<nprocsX*nprocsY);
                    int ierror = MPI_Recv(&b, sizeof(b), MPI_BYTE, senderID, blockID, MPI_COMM_WORLD, &status[blockID]);
                }
        }
    }
    
    inline void _setupMatrix()
    {
        /* Create an empty matrix object */
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
        
        /* Indicate that the matrix coefficients are ready to be set */
        HYPRE_StructMatrixInitialize(A);
        
        /* Set the matrix coefficients.  Each processor assigns coefficients
         for the boxes in the grid that it owns. Note that the coefficients
         associated with each stencil entry may vary from grid point to grid
         point if desired.  Here, we first set the same stencil entries for
         each grid point.  Then we make modifications to grid points near
         the boundary. */
        int ilower[2] = { 0, 0 };
        int iupper[2]= { sizeX-1, sizeY-1 };
        
        if (nranks>1)
        {
            ilower[0] = rankX * sX;
            ilower[1] = rankY * sY;
            iupper[0] = (rankX+1) * sX - 1;
            iupper[1] = (rankY+1) * sY - 1;
        }
        
        int stencil_indices[5] = {0,1,2,3,4}; /* labels for the stencil entries -
                                               these correspond to the offsets
                                               defined above */
        int nentries = 5;
        int nvalues  = size2 * nentries;
        vector<double> values(nvalues);
        
        /* We have size2 grid points, each with 5 stencil entries */
        if (bConstantCoefficients)
            _fillMatrixConstOMP<Lab>(values);
        else
            _fillMatrixVariableOMP<Lab>(values);
        
        HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                       stencil_indices, &values[0]);
        
        /* This is a collective call finalizing the matrix assembly.
         The matrix is now ``ready to be used'' */
        HYPRE_StructMatrixAssemble(A);
    }
    
    inline void _setupVectors()
    {
        /* Create an empty vector object */
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
        
        /* Indicate that the vector coefficients are ready to be set */
        HYPRE_StructVectorInitialize(b);
        HYPRE_StructVectorInitialize(x);
        
        /* Set the vector coefficients */
        int ilower[2] = { 0, 0 };
        int iupper[2]= { sizeX-1, sizeY-1 };
        
        int Nx = fluidgrid->getBlocksPerDimension(0);
        int Ny = fluidgrid->getBlocksPerDimension(1);
        
        if (nranks>1)
        {
            ilower[0] = rankX * sX;
            ilower[1] = rankY * sY;
            iupper[0] = (rankX+1) * sX - 1;
            iupper[1] = (rankY+1) * sY - 1;
            Nx = fluidgrid->getBlocksPerDimension(0)/nprocsX;
            Ny = fluidgrid->getBlocksPerDimension(1)/nprocsY;
        }
        
        vector<double> values(size2);
        
        vector<BlockInfo> myInfo = fluidgrid->getBlocksInfo();
        BlockInfo * ary = &myInfo.front();
        
        // right hand side - divU
#pragma omp parallel for schedule(static) collapse(2)
        for(int by=0; by<Ny; by++)
            for(int bx=0; bx<Nx; bx++)
            {
                const int i = bx+rankX*Nx + (by+rankY*Ny)*fluidgrid->getBlocksPerDimension(0);
                
                BlockInfo info = myInfo[i];
                FluidBlock& b = *(FluidBlock*)info.ptrBlock;
                
                const double h2 = info.h_gridpoint * info.h_gridpoint;
                
                for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                    {
                        const int gix = ix + info.index[0]*FluidBlock::sizeX - sX*rankX;
                        const int giy = iy + info.index[1]*FluidBlock::sizeY - sY*rankY;
                        const int idx = gix + giy*sX;
                        values[idx] = h2*b(ix,iy).divU;
                    }
            }
        
        HYPRE_StructVectorSetBoxValues(b, ilower, iupper, &values[0]);
        
        // initial guess
        /*
         #pragma omp parallel for schedule(static)
         for (int i = 0; i < size2; i ++)
         values[i] = 0.0;
         */
        
        HYPRE_StructVectorSetBoxValues(x, ilower, iupper, &values[0]);
        
        /* This is a collective call finalizing the vector assembly.
         The vectors are now ``ready to be used'' */
        HYPRE_StructVectorAssemble(b);
        HYPRE_StructVectorAssemble(x);
    }
    
public:
    MultigridHypre() : nprocsX(8), nprocsY(4), bSetup(false), bConstantCoefficients(false)
	{
		cout << "MG not implemented for 3D\n";
		abort();
	}
    
    ~MultigridHypre() {}
    
    void setup(FluidGridMPI * inputGrid, bool bConstCoeffs, const int rank, const int nprocs)
    {
        nranks = nprocs;
        fluidgrid = inputGrid;
        bConstantCoefficients = bConstCoeffs;
        
        vector<BlockInfo> myInfo = fluidgrid->getBlocksInfo();
        const int N = myInfo.size();
        
        sizeX = fluidgrid->getBlocksPerDimension(0) * FluidBlock::sizeX;
        sizeY = fluidgrid->getBlocksPerDimension(1) * FluidBlock::sizeY;
        size2 = sizeX*sizeY;
        
        // require ghosts! - bcast grid
        // broadcast from rank 0
        for (int j=0; j<fluidgrid->getBlocksPerDimension(1); j++)
            for (int i=0; i<fluidgrid->getBlocksPerDimension(0); i++)
            {
                FluidBlock &b = (*fluidgrid)(i,j);
                MPI_Bcast(&b, sizeof(b), MPI_BYTE, 0, MPI_COMM_WORLD);
            }
        
        /* 1. Set up a grid. Each processor describes the piece
         of the grid that it owns. */
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
        
        int ilower[2] = { 0, 0 };
        int iupper[2] = { sizeX-1, sizeY-1 };
        
        sX = sizeX;
        sY = sizeY;
        rankX = 0;
        rankY = 0;
        
        if (nranks>1)
        {
            if (nranks != nprocsX*nprocsY)
            {
                cout << "Currently only supporting 1 or 32 MPI processes\n";
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            
            if (fluidgrid->getBlocksPerDimension(0)<nprocsX)
            {
                cout << "Grid size too small\n";
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            rankX = rank % nprocsX; // 0-7
            rankY = rank / nprocsX; // 0-3
            sX = sizeX / nprocsX;
            sY = sizeY / nprocsY;
            size2 = sX*sY;
            ilower[0] = rankX * sX;
            ilower[1] = rankY * sY;
            iupper[0] = (rankX+1) * sX - 1;
            iupper[1] = (rankY+1) * sY - 1;
        }
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
        
#ifndef _PERIODIC_
        int periodicity[2] = { sizeX,0 };
        //int periodicity[2] = { 0,0 }; // for walled sides
#else // _PERIODIC_
        int periodicity[2] = { sizeX,sizeY };
#endif // _PERIODIC_
        HYPRE_StructGridSetPeriodic(grid, periodicity);
        
        HYPRE_StructGridAssemble(grid);
        
        
        /* 2. Define the discretization stencil */
        HYPRE_StructStencilCreate(2, 5, &stencil);
        
        int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
        
        for (int entry = 0; entry < 5; entry++)
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
        
        
        if (nranks>1)
            omp_set_num_threads(1);
        
        bSetup = true;
    }
    
    void operator()()
    {
        assert(bSetup);
        
        /* 3. Set up a Struct Matrix */
        _setupMatrix();
        
        /* 4. Set up Struct Vectors for b and x. */
        _setupVectors();
        
        /* 5. Set up and use a solver (See the Reference Manual for descriptions
         of all of the options.) */
        /* Create an empty PCG Struct solver */
        int n_pre = 1;
        int n_post = 1;
        int num_iterations;
        double final_res_norm;
        
        HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructSMGSetMemoryUse(solver, 0);
        HYPRE_StructSMGSetMaxIter(solver, 100);
        HYPRE_StructSMGSetTol(solver, 1.0e-08);
        HYPRE_StructSMGSetRelChange(solver, 0);
        HYPRE_StructSMGSetNonZeroGuess(solver);
        HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
        HYPRE_StructSMGSetNumPostRelax(solver, n_post);
        /* Logging must be on to get iterations and residual norm info below */
        HYPRE_StructSMGSetLogging(solver, 1);
        HYPRE_StructSMGSetPrintLevel(solver, 1);
        
        /* Setup and solve */
        HYPRE_StructSMGSetup(solver, A, b, x);
        HYPRE_StructSMGSolve(solver, A, b, x);
        
        /* Get some info on the run */
        //HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
        //HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        
        //printf("\n");
        //printf("Iterations = %d\n", num_iterations);
        //printf("Final Relative Residual Norm = %g\n", final_res_norm);
        //printf("\n");
        
        
        /* 6. Copy back the results */
        _getResultsOMP();
        
        HYPRE_StructGridDestroy(grid);
        HYPRE_StructStencilDestroy(stencil);
        HYPRE_StructMatrixDestroy(A);
        HYPRE_StructVectorDestroy(b);
        HYPRE_StructVectorDestroy(x);
        HYPRE_StructSMGDestroy(solver);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (nranks>1 && rankX==0 && rankY==0)
            omp_set_num_threads(NTHREADS);
    }
};

#endif
