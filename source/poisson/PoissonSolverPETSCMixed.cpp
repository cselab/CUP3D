//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#include "PoissonSolverPETSCMixed.h"
#ifdef CUP_PETSC

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
/*
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
*/
#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

extern PetscErrorCode ComputeRHS(KSP solver, Vec RHS, void * Sptr);
extern PetscErrorCode ComputeJacobian(KSP solver, Mat J, Mat JAC, void *Sptr);
std::vector<char*> readRunArgLst();

struct PoissonSolverMixed_PETSC::PetscData
{
  Real * const cupRHS;
  const size_t cupRHS_size;
  const int myNx, myNy, myNz;
  const int gNx, gNy, gNz;
  const Real h;
  const MPI_Comm m_comm;
  const BCflag BCx, BCy, BCz;
  KSP solver;
  DM grid;
  Vec SOL;

  PetscData(MPI_Comm c, size_t gx, size_t gy, size_t gz, size_t nx, size_t ny,
  size_t nz, size_t nsize, BCflag BCX, BCflag BCY, BCflag BCZ, Real H, Real*ptr)
  : cupRHS(ptr), cupRHS_size(nsize), myNx(nx), myNy(ny), myNz(nz), gNx(gx),
  gNy(gy), gNz(gz), h(H), m_comm(c), BCx(BCX), BCy(BCY), BCz(BCZ) { }
  ~PetscData() {
    VecDestroy(& SOL);
    DMDestroy(& grid);
    KSPDestroy(& solver);
  }
};

PoissonSolverMixed_PETSC::PoissonSolverMixed_PETSC(SimulationData&s) : PoissonSolver(s)
{
  PETSC_COMM_WORLD = m_comm;
  data = new Real[myN[0] * myN[1] * myN[2]];
  data_size = (size_t) myN[0] * (size_t) myN[1] * (size_t) myN[2];
  stridez = myN[1] * myN[0]; // slow
  stridey = myN[0];
  stridex = 1; // fast

  S= new PetscData(m_comm, gsize[0], gsize[1], gsize[2], myN[0], myN[1], myN[2],
    data_size, sim.BCx_flag, sim.BCy_flag, sim.BCz_flag, h, data);

  // (int argc,char **argv)
  std::vector<char*> args = readRunArgLst();
  int argc = args.size()-1; char ** argv = args.data();
  //static char help[] = "Why on earth would it ask me for a string now?\n";
  PetscInitialize(&argc, &argv, (char*)0, (char*)0);
  PetscErrorCode ierr = DMDACreate3d(m_comm,
    sim.BCx_flag == periodic ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
    sim.BCy_flag == periodic ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
    sim.BCz_flag == periodic ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, gsize[0], gsize[1], gsize[2],
    sim.nprocsx, sim.nprocsy, sim.nprocsz, 1, 1, NULL, NULL, NULL, & S->grid);
  //* */ CHKERRQ(ierr);
  // ierr = DMSetMatType(da_Stokes, MATAIJ); /* What is this for? */
  //* */ CHKERRQ(ierr);
  ierr = DMSetFromOptions(S->grid); /* I'm sure this is doing something.. */
  //* */ CHKERRQ(ierr);
  ierr = DMSetUp(S->grid); /* I just genuinely do not know what. */
  //* */ CHKERRQ(ierr);
  ierr = KSPSetType(S->solver, KSPCR);
  //* */ CHKERRQ(ierr);
  ierr = KSPCreate(m_comm, & S->solver);
  //* */ CHKERRQ(ierr);
  ierr = KSPSetFromOptions(S->solver);
  //* */ CHKERRQ(ierr);
  ierr = KSPSetDM(S->solver, S->grid);
  //* */ CHKERRQ(ierr);
  DMSetApplicationContext(S->grid, S); /* THIS IS GET! THIS I UNDESTAND! */

  KSPSetComputeRHS(S->solver, ComputeRHS, S);
  KSPSetComputeOperators(S->solver, ComputeJacobian, S);
  KSPSetInitialGuessNonzero(S->solver, PETSC_TRUE );
  KSPSetFromOptions(S->solver);
}

PoissonSolverMixed_PETSC::~PoissonSolverMixed_PETSC()
{
  delete S;
  PetscFinalize();
}

void PoissonSolverMixed_PETSC::solve()
{
  sim.startProfiler("PETSC cub2rhs");
  _cub2fftw();
  sim.stopProfiler();

  if(bRankHoldsFixedDOF)
  {
    // set last corner such that last point has pressure pLast
    data[fixed_idx]  = 6*h*pLast; // for conditioning = 1
    // neighbours read value of corner from the RHS:
    data[fixed_m1z] -= h*pLast; //fixed -1dz reads fixed dof (+1dz) from RHS
    data[fixed_p1z] -= h*pLast; //fixed +1dz reads fixed dof (-1dz) from RHS
    data[fixed_m1y] -= h*pLast; //fixed -1dy reads fixed dof (+1dy) from RHS
    data[fixed_p1y] -= h*pLast; //fixed +1dy reads fixed dof (-1dy) from RHS
    data[fixed_m1x] -= h*pLast; //fixed -1dx reads fixed dof (+1dx) from RHS
    data[fixed_p1x] -= h*pLast; //fixed +1dx reads fixed dof (-1dx) from RHS
  }

  KSPSolve(S->solver, NULL, NULL);

  {
    KSPGetSolution(S->solver, & S->SOL);
    PetscScalar ***array;
    DMDAVecGetArray(S->grid, S->SOL, & array);
    PetscInt xStart, yStart, zStart, xSpan, ySpan, zSpan;
    DMDAGetCorners(S->grid, &xStart,&yStart,&zStart, &xSpan,&ySpan,&zSpan);
    assert(myN[0]==(size_t)xSpan&&myN[1]==(size_t)ySpan&&myN[2]==(size_t)zSpan);
    sim.startProfiler("PETSC mean0");
    {
      Real avgP = 0; const Real fac = 1.0 / (gsize[0] * gsize[1] * gsize[2]);
      // Compute average pressure across all ranks:
      #pragma omp parallel for schedule(static) reduction(+ : avgP)
      for(int k=0;k<zSpan;k++) for(int j=0;j<ySpan;j++) for(int i=0;i<xSpan;i++)
        avgP += fac * array[k+zStart][j+yStart][i+xStart];
      MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPIREAL, MPI_SUM, m_comm);
      // Subtract average pressure from all gridpoints
      #pragma omp parallel for schedule(static) reduction(+ : avgP)
      for(int k=0;k<zSpan;k++) for(int j=0;j<ySpan;j++) for(int i=0;i<xSpan;i++)
      {
        array[k+zStart][j+yStart][i+xStart] -= avgP;
        data[i +myN[0]*j +myN[0]*myN[1]*k]= array[k+zStart][j+yStart][i+xStart];
      }
      // Save pressure of a corner of the grid so that it can be imposed next time
      pLast = array[zSpan-1+zStart][ySpan-1+yStart][xSpan-1+xStart];
      DMDAVecRestoreArray(S->grid, S->SOL, & array);
      printf("Avg Pressure:%f\n", avgP);
    }
    sim.stopProfiler();
  }

  sim.startProfiler("PETSC rhs2cub");
  _fftw2cub();
  sim.stopProfiler();
}

PetscErrorCode ComputeRHS(KSP solver, Vec RHS, void * Sptr)
{
  const auto& S = *( PoissonSolverMixed_PETSC::PetscData *) Sptr;
  const Real* const CubRHS = S.cupRHS;
  const size_t cSx = S.myNx, cSy = S.myNy;
  PetscScalar    ***array;

  PetscInt xStart, yStart, zStart, xSpan, ySpan, zSpan;
  DMDAGetCorners(S.grid, &xStart,&yStart,&zStart, &xSpan,&ySpan,&zSpan);
  DMDAVecGetArray(S.grid, RHS, & array);
  #pragma omp parallel for schedule(static)
  for (PetscInt k=0; k<zSpan; k++)
  for (PetscInt j=0; j<ySpan; j++)
  for (PetscInt i=0; i<xSpan; i++)
    array[k+zStart][j+yStart][i+xStart] = CubRHS[i + cSx*j + cSx*cSy*k];

  DMDAVecRestoreArray(S.grid, RHS, & array);
  VecAssemblyBegin(RHS);
  VecAssemblyEnd(RHS);
  // MatNullSpace   nullspace;
  // force right hand side to be consistent for singular matrix
  // it is just a hack that we avoid by fixing one DOF
  // MatNullSpaceCreate(S.m_comm, PETSC_TRUE, 0, 0, &nullspace);
  // MatNullSpaceRemove(nullspace,b);
  // MatNullSpaceDestroy(&nullspace);
  return 0;
}

PetscErrorCode ComputeJacobian(KSP solver, Mat J, Mat JAC, void *Sptr)
{
  const auto& S = *( PoissonSolverMixed_PETSC::PetscData *) Sptr;
  PetscInt xSt, ySt, zSt, xSpan, ySpan, zSpan;
  const Real h = S.h;
  DMDAGetCorners(S.grid, &xSt,&ySt,&zSt, &xSpan,&ySpan,&zSpan);
  #pragma omp parallel for schedule(static)
  for (int k=0; k<zSpan; k++)
  for (int j=0; j<ySpan; j++)
  for (int i=0; i<xSpan; i++)
  {
    MatStencil R, C[7]; R.k = k+zSt; R.j = j+ySt; R.i = i+xSt;
    PetscScalar V[7];
    V[0] = 6*h; C[0].i = xSt+i;   C[0].j = ySt+j;   C[0].k = zSt+k;
    V[1] = h;   C[1].i = xSt+i-1; C[1].j = ySt+j;   C[1].k = zSt+k;
    V[2] = h;   C[2].i = xSt+i+1; C[2].j = ySt+j;   C[2].k = zSt+k;
    V[3] = h;   C[3].i = xSt+i;   C[3].j = ySt+j-1; C[3].k = zSt+k;
    V[4] = h;   C[4].i = xSt+i;   C[4].j = ySt+j+1; C[4].k = zSt+k;
    V[5] = h;   C[5].i = xSt+i;   C[5].j = ySt+j;   C[5].k = zSt+k-1;
    V[6] = h;   C[6].i = xSt+i;   C[6].j = ySt+j;   C[6].k = zSt+k+1;

    if( S.BCz != periodic && zSt+k == S.gNz-1 ) { V[0] += h; V[6] = 0; }
    if( S.BCz != periodic && zSt+k ==       0 ) { V[0] += h; V[5] = 0; }
    if( S.BCy != periodic && ySt+j == S.gNy-1 ) { V[0] += h; V[4] = 0; }
    if( S.BCy != periodic && ySt+j ==       0 ) { V[0] += h; V[3] = 0; }
    if( S.BCx != periodic && xSt+i == S.gNx-1 ) { V[0] += h; V[2] = 0; }
    if( S.BCx != periodic && xSt+i ==       0 ) { V[0] += h; V[1] = 0; }

    if( zSt+k == S.gNz-2 && ySt+j == S.gNy-2 && xSt+i == S.gNx-2 ) {
      // set last corner such that last point has pressure pLast
      V[1]=0; V[2]=0; V[3]=0; V[4]=0; V[5]=0; V[6]=0; V[0]=6*h; // condition = 1
    } // neighbours read value of corner from the RHS:
    if( zSt+k == S.gNz-2 && ySt+j == S.gNy-2 && xSt+i == S.gNx-1 ) V[1]=0;
    if( zSt+k == S.gNz-2 && ySt+j == S.gNy-2 && xSt+i == S.gNx-3 ) V[2]=0;
    if( zSt+k == S.gNz-2 && ySt+j == S.gNy-1 && xSt+i == S.gNx-2 ) V[3]=0;
    if( zSt+k == S.gNz-2 && ySt+j == S.gNy-3 && xSt+i == S.gNx-2 ) V[4]=0;
    if( zSt+k == S.gNz-1 && ySt+j == S.gNy-2 && xSt+i == S.gNx-2 ) V[5]=0;
    if( zSt+k == S.gNz-3 && ySt+j == S.gNy-2 && xSt+i == S.gNx-2 ) V[6]=0;

    int nFilled = 7;
    if( std::fabs(V[6]) < 1e-16 ) { nFilled--;  // row 6 is unneeded
    }
    if( std::fabs(V[5]) < 1e-16 ) { nFilled--;  // row 5 is unneeded
      V[5] = V[6]; C[5].i = C[6].i; C[5].j = C[6].j; C[5].k = C[6].k; // 6 to 5
    }
    if( std::fabs(V[4]) < 1e-16 ) { nFilled--;  // row 4 is unneeded
      V[4] = V[5]; C[4].i = C[5].i; C[4].j = C[5].j; C[4].k = C[5].k; // 5 to 4
      V[5] = V[6]; C[5].i = C[6].i; C[5].j = C[6].j; C[5].k = C[6].k; // 6 to 5
    }
    if( std::fabs(V[3]) < 1e-16 ) { nFilled--;  // row 3 is unneeded
      V[3] = V[4]; C[3].i = C[4].i; C[3].j = C[4].j; C[3].k = C[4].k; // 4 to 3
      V[4] = V[5]; C[4].i = C[5].i; C[4].j = C[5].j; C[4].k = C[5].k; // 5 to 4
      V[5] = V[6]; C[5].i = C[6].i; C[5].j = C[6].j; C[5].k = C[6].k; // 6 to 5
    }
    if( std::fabs(V[2]) < 1e-16 ) { nFilled--;  // row 2 is unneeded
      V[2] = V[3]; C[2].i = C[3].i; C[2].j = C[3].j; C[2].k = C[3].k; // 3 to 2
      V[3] = V[4]; C[3].i = C[4].i; C[3].j = C[4].j; C[3].k = C[4].k; // 4 to 3
      V[4] = V[5]; C[4].i = C[5].i; C[4].j = C[5].j; C[4].k = C[5].k; // 5 to 4
      V[5] = V[6]; C[5].i = C[6].i; C[5].j = C[6].j; C[5].k = C[6].k; // 6 to 5
    }
    if( std::fabs(V[1]) < 1e-16 ) { nFilled--;  // row 1 is unneeded
      V[1] = V[2]; C[1].i = C[2].i; C[1].j = C[2].j; C[1].k = C[2].k; // 2 to 1
      V[2] = V[3]; C[2].i = C[3].i; C[2].j = C[3].j; C[2].k = C[3].k; // 3 to 2
      V[3] = V[4]; C[3].i = C[4].i; C[3].j = C[4].j; C[3].k = C[4].k; // 4 to 3
      V[4] = V[5]; C[4].i = C[5].i; C[4].j = C[5].j; C[4].k = C[5].k; // 5 to 4
      V[5] = V[6]; C[5].i = C[6].i; C[5].j = C[6].j; C[5].k = C[6].k; // 6 to 5
    }
    MatSetValuesStencil(JAC, 1, &R, nFilled, C, V, INSERT_VALUES);
  }
  MatAssemblyBegin(JAC, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(JAC, MAT_FINAL_ASSEMBLY);
  // MatNullSpace   nullspace;
  // MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
  // MatSetNullSpace(J,nullspace);
  // MatNullSpaceDestroy(&nullspace);
  return 0;
}

std::vector<char*> readRunArgLst()
{
  std::vector<char*> args;
  std::vector<std::string> params;
  params.push_back("-ksp_monitor_short");
  #if 1
    params.push_back("-pc_type"); params.push_back("mg");
    params.push_back("-pc_mg_type"); params.push_back("full");
    params.push_back("-ksp_type"); params.push_back("fgmres");
    params.push_back("-pc_mg_levels"); params.push_back("3");
    params.push_back("-mg_coarse_pc_factor_shift_type");
    params.push_back("nonzero");
  #else
    params.push_back("-pc_type"); params.push_back("ksp");
    // -da_grid_x 50 -da_grid_y 50 ???
    params.push_back("-ksp_ksp_type"); params.push_back("cg");
    params.push_back("-kspksp_pc_type_ksp_type"); params.push_back("bjacobi");
    params.push_back("-ksp_ksp_rtol"); params.push_back("1e-1");
    params.push_back("-ksp_ksp_monitor");
    params.push_back("-ksp_type"); params.push_back("pipefgmres");
  #endif
  for(size_t i; i<params.size(); i++) {
    char *arg = new char[params[i].size() + 1];
    copy(params[i].begin(), params[i].end(), arg);  // write into char array
    arg[params[i].size()] = '\0';
    args.push_back(arg);
  }
  args.push_back(0); // push back nullptr as last entry
  return args; // remember to deallocate it!
}

#endif
