//#include "PoissonSolverScalarFFTW_ACC.h"
//#include <cuda_runtime_api.h>
#ifndef _SP_COMP_
        typedef double Real;
#else
        typedef float Real;
#endif
typedef Real Complex[2];

__global__
void _fourier_filter_kernel(Complex *data_hat, const dim3 N, const dim3 isize,
                            const dim3 istart, const Real h)
{
    const Real waveFactX = 2.0*M_PI/(h*N.x);
    const Real waveFactY = 2.0*M_PI/(h*N.y);
    const Real waveFactZ = 2.0*M_PI/(h*N.z);
	const Real norm_factor = 1./Real(N.x*N.y*N.z);

	const int halfX = N.x/2;
	const int halfY = N.y/2;
	const int halfZ = N.z/2;

	// take care (direction reverse order for cuda)
	// cuda X dir maps k
	// cuda Y dir maps j
	const unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int j = blockDim.y * blockIdx.y + threadIdx.y;

	if (k < isize.z and j < isize.y)
	{
		const int ky = istart.y+j;
		const int kz = istart.z+k;
		const Real kky = (ky>halfY) ? Real(ky-N.y) : (Real)ky;
		const Real kkz = (kz>halfZ) ? Real(kz-N.z) : (Real)kz;
		const Real rky = kky*waveFactY;
		const Real rkz = kkz*waveFactZ;

		for (int i=0, index=j*isize.z+k; i<isize.x; i++, index+=isize.y*isize.z)
		{
			const int kx = istart.x+i;
			const Real kkx = (kx>halfX) ? Real(kx-N.x) : (Real)kx;
			const Real rkx = kkx*waveFactX;

			const Real kinv = (kkx==0 && kky==0 && kkz==0) ? 0.0
					: -norm_factor/(rkx*rkx+rky*rky+rkz*rkz);
			//const int index = (i*osize[1]+j)*osize[2]+k;
			data_hat[index][0] *= kinv;
			data_hat[index][1] *= kinv;
		}
	}
}

//extern "C" {
#define POISSON_FILTER_DIMX 16
#define POISSON_FILTER_DIMY 16
void _fourier_filter_gpu(Complex *data_hat, const int N[3], const int isize[3], const int istart[3], const Real h)
{
	const dim3 NN     (N[0]     , N[1]     , N[2]);
	const dim3 iisize (isize[0] , isize[1] , isize[2]);
	const dim3 iistart(istart[0], istart[1], istart[2]);

	// CUDA X dir maps isize[2]
	// CUDA Y dir maps isize[1]
	// isize[0] is sweeped inside kernel
	const int blocksInX = (isize[2]+POISSON_FILTER_DIMX-1)/POISSON_FILTER_DIMX;
	const int blocksInY = (isize[1]+POISSON_FILTER_DIMY-1)/POISSON_FILTER_DIMY;

	const dim3 DimGrid(blocksInX, blocksInY, 1);
	const dim3 DimBlock(POISSON_FILTER_DIMX, POISSON_FILTER_DIMY, 1);
	_fourier_filter_kernel<<<DimGrid, DimBlock>>>(data_hat, NN, iisize, iistart, h);
}
//}
