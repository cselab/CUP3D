#pragma once
//#include <mpi.h>
//#include <vector>
//#include <cassert>
//#include <cmath>
//#include <iostream>

#include "Definitions.h"
#include "accfft_utils.h"

#ifdef _CUDA_COMP_
#include <cuda_runtime_api.h>
	#ifndef _SP_COMP_
		#include "accfft_gpu.h"
		typedef accfft_plan_gpu myplan;
		typedef Complex myComplex;
	#else
		#include "accfft_gpuf.h"
		typedef accfft_plan_gpuf myplan;
		typedef Complexf myComplex;
	#endif
#else
	#ifndef _SP_COMP_
		#include "accfft.h"
		typedef accfft_plan myplan;
		typedef Complex myComplex;
	#else
		#include "accfftf.h"
		typedef accfft_planf myplan;
		typedef Complexf myComplex;
	#endif
#endif

using namespace std;

#ifdef _CUDA_COMP_
void _fourier_filter_gpu(myComplex *data_hat, const int N[3], const int isize[3], const int istart[3], const Real h);
#endif

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_ACC
{
   typedef typename TGrid::BlockType BlockType;
	const int bs[3], mybpd[3], totbpd[3];
	int nprocs, procid, isize[3],osize[3],istart[3],ostart[3], alloc_max;
	MPI_Comm c_comm;
	Real * rho;
#ifdef _CUDA_COMP_
	Real * rho_gpu;
#endif
	myComplex * phi_hat;
	myplan * plan;

	void _cub2fft(TGrid& grid, Real * out) const
	{
		vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
		const size_t N = local_infos.size();
		const size_t myN[3] = {
				mybpd[0]*bs[0],
				mybpd[1]*bs[1],
				mybpd[2]*bs[2]
		};

#pragma omp parallel for
		for(int i=0; i<N; ++i) {
			const BlockInfo info = local_infos[i];
			BlockType& b = *(BlockType*)info.ptrBlock;
			const int myIstart[3] = {
					bs[0]*info.index[0],
					bs[1]*info.index[1],
					bs[2]*info.index[2]
			};
			const size_t offset = myIstart[2]+
						   myN[2]*myIstart[1]+
					myN[2]*myN[1]*myIstart[0];

			for(int iz=0; iz<bs[0]; iz++)
			for(int iy=0; iy<bs[1]; iy++)
			for(int ix=0; ix<bs[2]; ix++) {
				const size_t dest_index = offset + iz + myN[2] * (iy + myN[1] * ix);
				assert(dest_index >= 0 && dest_index < myN[0] * myN[1] * myN[2]);
				TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
			}
		}
	}

	void _fft2cub(Real * out, TGrid& grid) const
	{
		vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
		const size_t N = local_infos.size();
		const size_t myN[3] = {
				mybpd[0]*bs[0],
				mybpd[1]*bs[1],
				mybpd[2]*bs[2]
		};

#pragma omp parallel for
		for(int i=0; i<N; ++i) {
			const BlockInfo info = local_infos[i];
			BlockType& b = *(BlockType*)info.ptrBlock;
			const int myIstart[3] = {
					bs[0]*info.index[0],
					bs[1]*info.index[1],
					bs[2]*info.index[2]
			};
			const size_t offset = myIstart[2]+
						   myN[2]*myIstart[1]+
					myN[2]*myN[1]*myIstart[0];

			for(int iz=0; iz<bs[0]; iz++)
			for(int iy=0; iy<bs[1]; iy++)
			for(int ix=0; ix<bs[2]; ix++) {
				const size_t src_index = offset + iz + myN[2] * (iy + myN[1] * ix);
				assert(src_index >= 0 && src_index < myN[0] * myN[1] * myN[2]);
				TStreamer::operate(&out[src_index], b.data[iz][iy][ix]);
			}
		}
	}

#ifndef _CUDA_COMP_
	void _fourier_filter(myComplex *data_hat, const Real h)
	{
		const int NX = totbpd[0]*bs[0];
		const int NY = totbpd[1]*bs[1];
		const int NZ = totbpd[2]*bs[2];
        const Real waveFactX = 2.0*M_PI/(h*NX);
        const Real waveFactY = 2.0*M_PI/(h*NY);
        const Real waveFactZ = 2.0*M_PI/(h*NZ);
		const Real norm_factor = 1./Real(NX*NY*NZ);

#pragma omp parallel for collapse(3)
		for (int i=0; i < osize[0]; i++)
		for (int j=0; j < osize[1]; j++)
		for (int k=0; k < osize[2]; k++) {
			const int kx = ostart[0]+i;
			const int ky = ostart[1]+j;
			const int kz = ostart[2]+k;
			const int kkx = (kx>NX/2) ? kx-NX : kx;
			const int kky = (ky>NY/2) ? ky-NY : ky;
			const int kkz = (kz>NZ/2) ? kz-NZ : kz;
			const Real rkx = kkx*waveFactX;
			const Real rky = kky*waveFactY;
			const Real rkz = kkz*waveFactZ;
			const Real kinv = (kkx==0 && kky==0 && kkz==0) ? 0.0
							: -norm_factor/(rkx*rkx+rky*rky+rkz*rkz);
			const int index = (i*osize[1]+j)*osize[2]+k;
			data_hat[index][0] *= kinv;
			data_hat[index][1] *= kinv;
		}
	}
#endif


public:
	PoissonSolverScalarFFTW_ACC(const int desired_threads, TGrid& grid)
	: bs{BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ},
	  mybpd{grid.getResidentBlocksPerDimension(0), grid.getResidentBlocksPerDimension(1), grid.getResidentBlocksPerDimension(2)},
	  totbpd{grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1), grid.getBlocksPerDimension(2)}
	{
		if (totbpd[2]!=mybpd[2]) {
			printf("Poisson solver assumes grid is distrubuted in x and y directions.\n");
			abort();
		}
#ifdef _CUDA_COMP_
		//printf("NO CUDA FOR YOU!\n" ); abort();
#endif
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		int n[3] = {totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2]};
		int c_dims[2] = { totbpd[0]/mybpd[0], totbpd[1]/mybpd[1] };
		assert(totbpd[0]%mybpd[0]==0 && totbpd[1]%mybpd[1]==0);
		accfft_create_comm(MPI_COMM_WORLD,c_dims,&c_comm);

		// Get the local pencil size and the allocation size
#ifdef _CUDA_COMP_
	#ifndef _SP_COMP_
		alloc_max = accfft_local_size_dft_r2c_gpu(n,isize,istart,osize,ostart,c_comm);
	#else
		alloc_max = accfft_local_size_dft_r2c_gpuf(n,isize,istart, osize, ostart,c_comm);
	#endif
#else
	#ifndef _SP_COMP_
		alloc_max = accfft_local_size_dft_r2c(n,isize,istart,osize,ostart,c_comm);
	#else
		alloc_max = accfft_local_size_dft_r2cf(n,isize,istart,osize,ostart,c_comm);
	#endif
#endif
      printf("[mpi rank %d] isize  %3d %3d %3d\n",procid,mybpd[0],mybpd[1],mybpd[2]);
	printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d\n", procid,
			isize[0],isize[1],isize[2],
			osize[0],osize[1],osize[2]
	);
	printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", procid,
			istart[0],istart[1],istart[2],
			ostart[0],ostart[1],ostart[2]
	);
      assert(isize[0]==n[0] && isize[1]==n[1] && isize[2]==n[2]);

#ifdef _CUDA_COMP_
		rho=(Real*)malloc(isize[0]*isize[1]*isize[2]*sizeof(Real));
		cudaMalloc((void**) &rho_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real));
		cudaMalloc((void**) &phi_hat, alloc_max);
#else
		rho=(Real*)accfft_alloc(isize[0]*isize[1]*isize[2]*sizeof(Real));
		phi_hat=(myComplex*)accfft_alloc(alloc_max);
#endif

#ifdef _CUDA_COMP_
	#ifndef _SP_COMP_
			plan = accfft_plan_dft_3d_r2c_gpu(n, rho_gpu, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
	#else
			plan = accfft_plan_dft_3d_r2c_gpuf(n, rho_gpu, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
	#endif
#else
			accfft_init(desired_threads);
	#ifndef _SP_COMP_
			plan = accfft_plan_dft_3d_r2c(n, rho, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
	#else
			plan = accfft_plan_dft_3d_r2cf(n, rho, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
	#endif
#endif

		if (TStreamer::channels != 1) {
			cout << "PoissonSolverScalar_MPI(): Error: TStreamer::channels is " << TStreamer::channels << " (should be 1).\n";
			abort();
		}
	}

	void solve(TGrid& grid)
	{
		_cub2fft(grid, rho);
#ifdef _CUDA_COMP_
		cudaMemcpy(rho_gpu, rho, isize[0]*isize[1]*isize[2]*sizeof(Real), cudaMemcpyHostToDevice);
#endif

		// Perform forward FFT
		MPI_Barrier(c_comm);

#ifdef _CUDA_COMP_
	#ifndef _SP_COMP_
		accfft_execute_r2c_gpu(plan,rho_gpu,phi_hat);
	#else
		accfft_execute_r2c_gpuf(plan,rho_gpu,phi_hat);
	#endif
#else
	#ifndef _SP_COMP_
		accfft_execute_r2c(plan,rho,phi_hat);
	#else
		accfft_execute_r2cf(plan,rho,phi_hat);
	#endif
#endif

		// Spectral solve
		MPI_Barrier(c_comm);
		const Real h = grid.getBlocksInfo().front().h_gridpoint;
#ifdef _CUDA_COMP_
		const int NN[3] = {totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2]}; 
		_fourier_filter_gpu(phi_hat, NN, osize, ostart, h);
#else
		_fourier_filter(phi_hat, h);
#endif

		// Perform backward FFT
#ifdef _CUDA_COMP_
	#ifndef _SP_COMP_
		accfft_execute_c2r_gpu(plan,phi_hat,rho_gpu);
	#else
		accfft_execute_c2r_gpuf(plan,phi_hat,rho_gpu);
	#endif
#else
	#ifndef _SP_COMP_
		accfft_execute_c2r(plan,phi_hat,rho);
	#else
		accfft_execute_c2rf(plan,phi_hat,rho);
	#endif
#endif

#ifdef _CUDA_COMP_
		cudaMemcpy(rho, rho_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real), cudaMemcpyDeviceToHost);
#endif
		_fft2cub(rho, grid);
	}

	~PoissonSolverScalarFFTW_ACC()
	{
#ifndef _CUDA_COMP_
		accfft_free(rho);
		accfft_free(phi_hat);
		accfft_destroy_plan(plan);
		accfft_cleanup();
#else
		free(rho);
		cudaFree(rho_gpu);
		cudaFree(phi_hat);
		accfft_destroy_plan_gpu(plan);
		accfft_cleanup_gpu();
#endif
		MPI_Comm_free(&c_comm);
	}
};
