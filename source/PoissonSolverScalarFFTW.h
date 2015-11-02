#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

#include <fftw3.h>

#include "common.h"

#define _DCT_

#ifndef _SP_COMP_
// double
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#else // _SP_COMP_
// float
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#endif // _SP_COMP_

using namespace std;

#include <BlockInfo.h>
#include <Profiler.h>

class FFTWBase
{
	static int registered_objects;
	static bool initialized; //a la singleton
	
	static void _setup(const int desired_threads)
	{
		if (!initialized)
		{
			initialized = true;
			
#ifndef _SP_COMP_
			const int retval = fftw_init_threads();
#else // _SP_COMP_
			const int retval = fftwf_init_threads();
#endif // _SP_COMP_
			if(retval==0)
			{
				cout << "FFTWBase::setup(): Oops the call to fftw_init_threads() returned zero. Aborting\n";
				abort();
			}
			else
#ifndef _SP_COMP_
				fftw_plan_with_nthreads(desired_threads);
#else // _SP_COMP_
                fftwf_plan_with_nthreads(desired_threads);
#endif // _SP_COMP_
		}
		
		registered_objects++;
	}
	
public:
	
	FFTWBase(const int desired_threads) { _setup(desired_threads); }
	
	static void dispose()
	{
		registered_objects--;
		
		if (registered_objects == 0)
		{
#ifndef _SP_COMP_
			fftw_cleanup_threads();
#else // _SP_COMP_
			fftwf_cleanup_threads();
#endif // _SP_COMP_
		}
	}
};

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW : FFTWBase
{
	Profiler profiler;
	
protected:
	typedef typename TGrid::BlockType BlockType;
	
	bool initialized;
	
	myplan fwd, bwd;
	Real * data; // rhs in _setup, out in cub2fftw and fftw2cub
	
protected:
	
	virtual void _setup(Real *& rhs, const size_t nx, const size_t ny)
	{
		if (!initialized)
		{
			initialized = true;
#ifndef _SP_COMP_
			rhs = fftw_alloc_real(2*nx*(ny/2+1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
			
			fwd = fftw_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
			bwd = fftw_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#else // _SP_COMP_
			rhs = fftwf_alloc_real(2*nx*(ny/2+1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
			
			fwd = fftwf_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
			bwd = fftwf_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#endif // _SP_COMP_
		}
	}
	
	void _cub2fftw(TGrid& grid, Real * out, const size_t nx, const size_t ny, const size_t ny_hat) const
	{
		vector<BlockInfo> infos = grid.getBlocksInfo();
		
		const size_t N = infos.size();
		
#pragma omp parallel for
		for(int i=0; i<N; ++i)
		{
			const BlockInfo info = infos[i];
			BlockType& b = *(BlockType*)infos[i].ptrBlock;
			
			// indexing for fftw - [x][y] (y is the fastest running index), block part of the index
			const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
			
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					// indexing for fftw - [x][y] (y is the fastest running index), local part of the index
					const size_t dest_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
					
					assert(dest_index>=0 && dest_index<nx*ny_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
					TStreamer::operate(b(ix,iy), &out[dest_index]);
				}
		}
	}
	
	virtual void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
	{
		if (TStreamer::channels != 1)
		{
			cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
			abort();
		}
		
		const Real h2 = h*h;
		const Real factor = h2*norm_factor;
		
#pragma omp parallel for
		for(int i=0; i<nx; ++i)
			for(int j = 0; j<ny_hat; ++j)
			{
				const int linidx = i*ny_hat+j;
				assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
				
				// based on the 5 point stencil in 1D (h^4 error)
				const Real denom = 32.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 2.*(cos(4.*M_PI*i/nx) + cos(4.*M_PI*j/ny)) - 60.;
				const Real inv_denom = (denom==0)? 0.:1./denom;
				const Real fatfactor = 12. * inv_denom * factor;
				
				// based on the 3 point stencil in 1D (h^2 error)
				//const Real denom = 2.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 4.;
				//const Real inv_denom = (denom==0)? 0.:1./denom;
				//const Real fatfactor = inv_denom * factor;
				
				// this is to check the transform only
				//const Real fatfactor = norm_factor;
				
				in_out[linidx][0] *= fatfactor;
				in_out[linidx][1] *= fatfactor;
			}
		
		//this is sparta!
		in_out[0][0] = in_out[0][1] = 0;
	}
	
	virtual void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
	{
		if (TStreamer::channels != 1)
		{
			cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
			abort();
		}
		
		const Real waveFactX = 2.0*M_PI/(nx*h);
		const Real waveFactY = 2.0*M_PI/(ny*h);
#pragma omp parallel for
		for(int i=0; i<nx; ++i)
			for(int j = 0; j<ny_hat; ++j)
			{
				const int linidx = i*ny_hat+j;
				assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
				
				// wave number
				const int kx = (i <= nx/2) ? i : -(nx-i);
				const int ky = (j <= ny/2) ? j : -(ny-j);
				const Real rkx = kx*waveFactX;
				const Real rky = ky*waveFactY;
				const Real kinv = (kx==0 && ky==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky);
				in_out[linidx][0] *= kinv*norm_factor;
				in_out[linidx][1] *= kinv*norm_factor;
			}
		
		//this is sparta!
		in_out[0][0] = in_out[0][1] = 0;
	}
	
	void _fftw2cub(Real * out, TGrid& grid, const size_t nx, const size_t ny, const size_t ny_hat) const
	{
		vector<BlockInfo> infos = grid.getBlocksInfo();
		
		const size_t N = infos.size();
		
		const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
		
#pragma omp parallel for
		for(int i=0; i<N; ++i)
		{
			const BlockInfo info = infos[i];
			BlockType& b = *(BlockType*)infos[i].ptrBlock;
			
			const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
			
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					const size_t src_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
					
					assert(src_index>=0 && src_index<nx*ny_hat*2); // linking error with openmp
					TStreamer::operate(&out[src_index], b(ix,iy));
				}
		}
	}
	
public:
	
	PoissonSolverScalarFFTW(const int desired_threads, TGrid& grid): FFTWBase(desired_threads), initialized(false)
    {
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1]);
        profiler.pop_stop();
    }
	
	void solve(TGrid& grid, const bool spectral=false)
	{
		const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
		const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
		
		profiler.push_start("CUB2FFTW");
		_cub2fftw(grid, data, gsize[0], gsize[1], gsize[1]/2+1);
		profiler.pop_stop();
		
		profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
		fftw_execute(fwd);
#else // _SP_COMP_
		fftwf_execute(fwd);
#endif // _SP_COMP_
		profiler.pop_stop();
		
		const Real norm_factor = 1./(gsize[0]*gsize[1]);
		const Real h = grid.getBlocksInfo().front().h_gridpoint;
		assert(1./gsize[0]==h);
		
		profiler.push_start("SOLVE");
		if(spectral)
			_solveSpectral((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
		else
			_solve((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
		
		profiler.pop_stop();
		
		profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
		fftw_execute(bwd);
#else // _SP_COMP_
		fftwf_execute(bwd);
#endif // _SP_COMP_
		profiler.pop_stop();
		
		profiler.push_start("FFTW2CUB");
		_fftw2cub(data, grid, gsize[0], gsize[1], gsize[1]/2+1);
		profiler.pop_stop();
		
		//profiler.printSummary();
	}
	
	void dispose()
	{
		if (initialized)
		{
			initialized = false;
			
#ifndef _SP_COMP_
			fftw_destroy_plan(fwd);
			fftw_destroy_plan(bwd);
			
			fftw_free(data);
#else // _SP_COMP_
			fftwf_destroy_plan(fwd);
			fftwf_destroy_plan(bwd);
			
			fftwf_free(data);
#endif // _SP_COMP_
			FFTWBase::dispose();
		}
	}
};

/*
template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_DCT : FFTWBase
{
    Profiler profiler;
    
protected:
    typedef typename TGrid::BlockType BlockType;
    
    bool initialized;
    
    myplan fwd, bwd;
    Real * data; // rhs in _setup, out in cub2fftw and fftw2cub
    
protected:
    
    virtual void _setup(Real *& rhs, const size_t nx, const size_t ny)
    {
        if (!initialized)
        {
            initialized = true;
#ifndef _SP_COMP_
            abort();
            rhs = fftw_alloc_real(4*(2*nx*(ny/2+1)-1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
            
            fwd = fftw_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
            bwd = fftw_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#else // _SP_COMP_
            rhs = fftwf_alloc_real(4*(2*nx*(ny/2+1)-1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
            
            fwd = fftwf_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
            bwd = fftwf_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#endif // _SP_COMP_
        }
    }
    
    void _cub2fftw(TGrid& grid, Real * out, const size_t nx, const size_t ny, const size_t ny_hat) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
            const size_t offset = BlockType::sizeY*info.index[1] + 4*((ny/2+1)*2-1)*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index1 = TStreamer::channels*(offset + iy + 4*(2*(ny/2+1)-1)*ix);
                    
                    assert(dest_index1>=0 && dest_index1<nx*ny_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(ix,iy), &out[dest_index1]);
                }
        }
        
        for (int i=0; i<nx; i++)
            for (int j=0; j<(ny/2+1)*2; j++)
            {
                const int baseIdx = i*4*(2*(ny/2+1)-1)                    + j;
                const int idx1    = i*4*(2*(ny/2+1)-1) + 2*(2*(ny/2+1)-1) - j;
                const int idx2    = i*4*(2*(ny/2+1)-1) + 2*(2*(ny/2+1)-1) + j;
                const int idx3    = i*4*(2*(ny/2+1)-1) + 4*(2*(ny/2+1)-1) - j;
                
                out[idx1] = out[baseIdx];
                out[idx2] = out[baseIdx];
                if (j!=0)
                    out[idx3] = out[baseIdx];
            }
    }
    
    virtual void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<2*(2*ny_hat-1); ++j)
            {
                const int linidx = i*2*(2*ny_hat-1)+j;
                assert(linidx >=0 && linidx<nx*2*(2*ny_hat-1)); // linking error with openmp
                
                // based on the 5 point stencil in 1D (h^4 error)
                const Real denom = 32.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 2.*(cos(4.*M_PI*i/nx) + cos(4.*M_PI*j/ny)) - 60.;
                const Real inv_denom = (denom==0)? 0.:1./denom;
                //const Real fatfactor = 12. * inv_denom * factor;
                
                // based on the 3 point stencil in 1D (h^2 error)
                //const Real denom = 2.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 4.;
                //const Real inv_denom = (denom==0)? 0.:1./denom;
                //const Real fatfactor = inv_denom * factor;
                
                // this is to check the transform only
                const Real fatfactor = norm_factor;
                
                in_out[linidx][0] *= fatfactor;
                in_out[linidx][1] *= fatfactor;
            }
        
        //this is sparta!
        //in_out[0][0] = in_out[0][1] = 0;
    }
    
    virtual void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
    {
        abort();
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real waveFactX = 2.0*M_PI/(nx*h);
        const Real waveFactY = 2.0*M_PI/(ny*h);
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<2*(2*ny_hat-1); ++j)
            {
                const int linidx = i*2*(2*ny_hat-1)+j;
                assert(linidx >=0 && linidx<nx*2*(2*ny_hat-1)); // linking error with openmp
                
                // wave number
                const int kx = (i <= nx/2) ? i : -(nx-i);
                const int ky = (j <= ny/2) ? j : -(ny-j);
                const Real rkx = kx*waveFactX;
                const Real rky = ky*waveFactY;
                const Real kinv = (kx==0 && ky==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky);
                in_out[linidx][0] *= kinv*norm_factor;
                in_out[linidx][1] *= kinv*norm_factor;
            }
        
        //this is sparta!
        in_out[0][0] = in_out[0][1] = 0;
    }
    
    void _fftw2cub(Real * out, TGrid& grid, const size_t nx, const size_t ny, const size_t ny_hat) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
        const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            const size_t offset = BlockType::sizeY*info.index[1] + 4*((ny/2+1)*2-1)*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = TStreamer::channels*(offset + iy + 4*(2*(ny/2+1)*ix-1));
                    
                    assert(src_index>=0 && src_index<nx*4*(ny_hat*2-1)); // linking error with openmp
                    TStreamer::operate(&out[src_index], b(ix,iy));
                }
        }
    }
    
public:
    
    PoissonSolverScalarFFTW_DCT(const int desired_threads, TGrid& grid): FFTWBase(desired_threads), initialized(false)
    {
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1]);
        profiler.pop_stop();
    }
    
    void solve(TGrid& grid, const bool spectral=false)
    {
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("CUB2FFTW");
        _cub2fftw(grid, data, gsize[0], gsize[1], gsize[1]/2+1);
        profiler.pop_stop();
        
        profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
        fftw_execute(fwd);
#else // _SP_COMP_
        fftwf_execute(fwd);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        const Real norm_factor = 1./(gsize[0]*gsize[1]);
        const Real h = grid.getBlocksInfo().front().h_gridpoint;
        assert(1./gsize[0]==h);
        
        profiler.push_start("SOLVE");
        if(spectral)
            _solveSpectral((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
        else
            _solve((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
        
        profiler.pop_stop();
        
        profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
        fftw_execute(bwd);
#else // _SP_COMP_
        fftwf_execute(bwd);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        profiler.push_start("FFTW2CUB");
        _fftw2cub(data, grid, gsize[0], gsize[1], gsize[1]/2+1);
        profiler.pop_stop();
        
        //profiler.printSummary();
    }
    
    void dispose()
    {
        if (initialized)
        {
            initialized = false;
            
#ifndef _SP_COMP_
            fftw_destroy_plan(fwd);
            fftw_destroy_plan(bwd);
            
            fftw_free(data);
#else // _SP_COMP_
            fftwf_destroy_plan(fwd);
            fftwf_destroy_plan(bwd);
            
            fftwf_free(data);
#endif // _SP_COMP_
            FFTWBase::dispose();
        }
    }
};
*/

// 2D, FFT+DCT in 1 step
//*
template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_DCT : FFTWBase
{
    Profiler profiler;
    
protected:
    typedef typename TGrid::BlockType BlockType;
    
    bool initialized;
    
    myplan fwd, bwd;
    Real * data; // rhs in _setup, out in cub2fftw and fftw2cub
    
protected:
    
    virtual void _setup(Real *& rhs, const size_t nx, const size_t ny)
    {
        if (!initialized)
        {
            initialized = true;
            
            //  Forward X:  FFT
            //  Backward X: FFT
            //  Forward Y:  DST-III
            //  Backward Y: DST-II
#ifndef _SP_COMP_
            rhs = fftw_alloc_real(nx*ny);
            
            fwd = fftw_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_R2HC, FFTW_REDFT11, FFTW_MEASURE);
            bwd = fftw_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_HC2R, FFTW_REDFT11, FFTW_MEASURE);
            //fwd = fftw_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_R2HC, FFTW_REDFT01, FFTW_MEASURE);
            //bwd = fftw_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_HC2R, FFTW_REDFT10, FFTW_MEASURE); // is the inverse correct? taken from http://www.manpagez.com/info/fftw3/fftw3-3.1.2/fftw3_35.php#SEC42
#else // _SP_COMP_
            rhs = fftwf_alloc_real(nx*ny);
            
            //fwd = fftwf_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_R2HC, FFTW_REDFT01, FFTW_MEASURE);
            //bwd = fftwf_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_HC2R, FFTW_REDFT10, FFTW_MEASURE);
            fwd = fftwf_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_R2HC, FFTW_REDFT11, FFTW_MEASURE);
            bwd = fftwf_plan_r2r_2d(nx, ny, rhs, rhs, FFTW_HC2R, FFTW_REDFT11, FFTW_MEASURE);
#endif // _SP_COMP_
        }
    }
    
    void _cub2fftw(TGrid& grid, Real * out, const size_t nx, const size_t ny) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
            const size_t offset = BlockType::sizeY*info.index[1] + ny*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index = TStreamer::channels*(offset + iy + ny*ix);
                    
                    assert(dest_index>=0 && dest_index<nx*ny); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(ix,iy), &out[dest_index]);
                }
        }
    }
    
    virtual void _solve(Real * in_out, const size_t nx, const size_t ny, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<ny; ++j)
            {
                const int linidx = i*ny+j;
                assert(linidx >=0 && linidx<nx*ny); // linking error with openmp
                
                // based on the 5 point stencil in 1D (h^4 error)
                //const Real denomY = 32.*cos(1.5*M_PI*(j+.5)/ny) - 2.*cos(2.5*M_PI*(j+.5)/ny) - 30.;
                const Real denomY = 32.*cos(M_PI*(j+.5)/ny) - 2.*cos(2*M_PI*(j+.5)/ny) - 30.;
                const Real denomX = 32.*cos(2.*M_PI*i/nx)   - 2.*cos(4.*M_PI*i/nx)     - 30.;
                const Real denom = denomX + denomY;
                const Real inv_denom = (denom==0)? 0.:1./denom;
                const Real fatfactor = 12. * inv_denom * factor;
                
                // based on the 3 point stencil in 1D (h^2 error)
                //const Real denom = 2.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 4.;
                //const Real inv_denom = (denom==0)? 0.:1./denom;
                //const Real fatfactor = inv_denom * factor;
                
                // this is to check the transform only
                //const Real fatfactor = norm_factor;
                in_out[linidx] *= fatfactor;
            }
        
        //this is sparta!
        //in_out[0] = 0; // WTF? this prevents correct fw/bw if not transformations happen, but has no influence on the 2nd derivative
        // this element should be removed for the FFT but not for the DCT! how to solve this issue?
    }
    
    virtual void _solveSpectral(Real * in_out, const size_t nx, const size_t ny, const Real norm_factor, const Real h)
    {
        abort();
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real waveFactX = 2.0*M_PI/(nx*h);
        const Real waveFactY = 2.0*M_PI/(ny*h);
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<ny; ++j)
            {
                const int linidx = i*ny+j;
                assert(linidx >=0 && linidx<nx*ny); // linking error with openmp
                
                // wave number
                const int kx = (i <= nx/2) ? i : -(nx-i);
                const int ky = (j <= ny/2) ? j : -(ny-j);
                const Real rkx = kx*waveFactX;
                const Real rky = ky*waveFactY;
                const Real kinv = (kx==0 && ky==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky);
                in_out[linidx] *= kinv*norm_factor;
            }
        
        //this is sparta!
        in_out[0] = 0;
    }
    
    void _fftw2cub(Real * out, TGrid& grid, const size_t nx, const size_t ny) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
        const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            const size_t offset = BlockType::sizeY*info.index[1] + ny*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = TStreamer::channels*(offset + iy + ny*ix);
                    
                    assert(src_index>=0 && src_index<nx*ny); // linking error with openmp
                    TStreamer::operate(&out[src_index], b(ix,iy));
                }
        }
    }
    
public:
    
    PoissonSolverScalarFFTW_DCT(const int desired_threads, TGrid& grid): FFTWBase(desired_threads), initialized(false)
    {
        // dimensions of blocks and whole grid
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1]);
        profiler.pop_stop();
    }
    
    void solve(TGrid& grid, const bool spectral=false)
    {
        // dimensions of blocks and whole grid
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("CUB2FFTW");
        _cub2fftw(grid, data, gsize[0], gsize[1]);
        profiler.pop_stop();
        
        profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
        fftw_execute(fwd);
#else // _SP_COMP_
        fftwf_execute(fwd);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        const Real norm_factor = .5/(gsize[0]*gsize[1]);
        const Real h = grid.getBlocksInfo().front().h_gridpoint;
        assert(1./gsize[0]==h);
        
        profiler.push_start("SOLVE");
        if(spectral)
            _solveSpectral(data, gsize[0], gsize[1], norm_factor, h);
        else
            _solve(data, gsize[0], gsize[1], norm_factor, h);
        
        profiler.pop_stop();
        
        profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
        fftw_execute(bwd);
#else // _SP_COMP_
        fftwf_execute(bwd);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        profiler.push_start("FFTW2CUB");
        _fftw2cub(data, grid, gsize[0], gsize[1]);
        profiler.pop_stop();
        
        //profiler.printSummary();
    }
    
    void dispose()
    {
        if (initialized)
        {
            initialized = false;
            
#ifndef _SP_COMP_
            fftw_destroy_plan(fwd);
            fftw_destroy_plan(bwd);
            
            fftw_free(data);
#else // _SP_COMP_
            fftwf_destroy_plan(fwd);
            fftwf_destroy_plan(bwd);
            
            fftwf_free(data);
#endif // _SP_COMP_
            FFTWBase::dispose();
        }
    }
};
//*/
/*
// 2D - separated FFT, DCT
template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_DCT : FFTWBase
{
    Profiler profiler;
    
protected:
    typedef typename TGrid::BlockType BlockType;
    
    bool initialized;
    
    myplan fwdx, fwdy, bwdx, bwdy;
    Real * datax, * datay; // rhs in _setup, out in cub2fftw and fftw2cub
    
protected:
    
    virtual void _setup(Real *& rhsx, Real *& rhsy, const size_t nx, const size_t ny)
    {
        if (!initialized)
        {
            initialized = true;
            
            //  Forward X:  FFT
            //  Backward X: FFT
            //  Forward Y:  DST-III
            //  Backward Y: DST-II
#ifndef _SP_COMP_
            abort();
            rhsx = fftw_alloc_real(2*nx*(ny/2+1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
            
            fwdx = fftw_plan_dft_r2c_2d(nx, ny, rhsx, (mycomplex *)rhsx, FFTW_MEASURE);
            bwdx = fftw_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhsx, rhsx, FFTW_MEASURE);
            
            rhsy = fftw_alloc_real(nx*ny);
            
            fwdy = fftw_plan_r2r_2d(nx, ny, rhsy, rhsy, FFTW_R2HC, FFTW_REDFT01, FFTW_MEASURE);
            bwdy = fftw_plan_r2r_2d(nx, ny, rhsy, rhsy, FFTW_HC2R, FFTW_REDFT10, FFTW_MEASURE); // is the inverse correct? taken from http://www.manpagez.com/info/fftw3/fftw3-3.1.2/fftw3_35.php#SEC42
#else // _SP_COMP_
            rhsx = fftwf_alloc_real(2*nx*(ny/2+1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
            
            fwdx = fftwf_plan_dft_r2c_2d(nx, ny, rhsx, (mycomplex *)rhsx, FFTW_MEASURE);
            bwdx = fftwf_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhsx, rhsx, FFTW_MEASURE);
            
            rhsy = fftwf_alloc_real(nx*ny);
            
            fwdy = fftwf_plan_r2r_2d(nx, ny, rhsy, rhsy, FFTW_R2HC, FFTW_REDFT01, FFTW_MEASURE);
            bwdy = fftwf_plan_r2r_2d(nx, ny, rhsy, rhsy, FFTW_HC2R, FFTW_REDFT10, FFTW_MEASURE);
#endif // _SP_COMP_
        }
    }
    
    void _cub2fftwX(TGrid& grid, Real * out, const size_t nx, const size_t ny, const size_t ny_hat) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
            const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
                    
                    assert(dest_index>=0 && dest_index<nx*ny_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(ix,iy), &out[dest_index]);
                }
        }
    }
    
    void _cub2fftwY(TGrid& grid, Real * out, const size_t nx, const size_t ny) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
            const size_t offset = BlockType::sizeY*info.index[1] + ny*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index = TStreamer::channels*(offset + iy + ny*ix);
                    
                    assert(dest_index>=0 && dest_index<nx*ny); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(ix,iy), &out[dest_index]);
                }
        }
    }
    
    virtual void _solveX(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<ny_hat; ++j)
            {
                const int linidx = i*ny_hat+j;
                assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
                
                // based on the 5 point stencil in 1D (h^4 error)
                const Real denom = 32.*cos(2.*M_PI*i/nx) - 2.*cos(4.*M_PI*i/nx) - 30.;
                const Real inv_denom = (denom==0)? 0.:1./denom;
                const Real fatfactor = 12. * inv_denom * factor;
                
                // this is to check the transform only
                //const Real fatfactor = norm_factor;
                
                in_out[linidx][0] *= fatfactor;
                in_out[linidx][1] *= fatfactor;
            }
        
        //this is sparta!
        in_out[0][0] = in_out[0][1] = 0;
    }
    
    virtual void _solveY(Real * in_out, const size_t nx, const size_t ny, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<ny; ++j)
            {
                const int linidx = i*ny+j;
                assert(linidx >=0 && linidx<nx*ny); // linking error with openmp
                
                // based on the 5 point stencil in 1D (h^4 error)
                const Real denom = 32.*cos(M_PI*(j+.5)/ny) - 2.*cos(2*M_PI*(j+.5)/ny) - 30.;
                const Real inv_denom = (denom==0)? 0.:1./denom;
                const Real fatfactor = 12. * inv_denom * factor;
                
                // this is to check the transform only
                //const Real fatfactor = norm_factor;
                in_out[linidx] *= fatfactor;
            }
    }
    
    void _fftw2cubX(Real * out, TGrid& grid, const size_t nx, const size_t ny, const size_t ny_hat) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
        const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
                    
                    assert(src_index>=0 && src_index<nx*ny_hat*2); // linking error with openmp
                    TStreamer::setTmp(&out[src_index], b(ix,iy));
                }
        }
    }
    
    void _fftw2cubY(Real * out, TGrid& grid, const size_t nx, const size_t ny) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
        const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            const size_t offset = BlockType::sizeY*info.index[1] + ny*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = TStreamer::channels*(offset + iy + ny*ix);
                    
                    assert(src_index>=0 && src_index<nx*ny); // linking error with openmp
                    TStreamer::operate(&out[src_index], b(ix,iy));
                }
        }
    }
    
    void _FFTx(TGrid& grid)
    {
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("CUB2FFTW");
        _cub2fftwX(grid, datax, gsize[0], gsize[1], gsize[1]/2+1);
        profiler.pop_stop();
        
        profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
        fftw_execute(fwdx);
#else // _SP_COMP_
        fftwf_execute(fwdx);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        const Real norm_factor = 1./(gsize[0]*gsize[1]);
        const Real h = grid.getBlocksInfo().front().h_gridpoint;
        assert(1./gsize[0]==h);
        
        profiler.push_start("SOLVE");
        _solveX((mycomplex *)datax, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
        profiler.pop_stop();
        
        profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
        fftw_execute(bwdx);
#else // _SP_COMP_
        fftwf_execute(bwdx);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        profiler.push_start("FFTW2CUB");
        _fftw2cubX(datax, grid, gsize[0], gsize[1], gsize[1]/2+1);
        profiler.pop_stop();
    }
    
    void _DCTy(TGrid& grid)
    {
        // dimensions of blocks and whole grid
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("CUB2FFTW");
        _cub2fftwY(grid, datay, gsize[0], gsize[1]);
        profiler.pop_stop();
        
        profiler.push_start("DCT FORWARD");
#ifndef _SP_COMP_
        fftw_execute(fwdy);
#else // _SP_COMP_
        fftwf_execute(fwdy);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        const Real norm_factor = .5/(gsize[0]*gsize[1]);
        const Real h = grid.getBlocksInfo().front().h_gridpoint;
        assert(1./gsize[0]==h);
        
        profiler.push_start("SOLVE");
        _solveY(datay, gsize[0], gsize[1], norm_factor, h);
        profiler.pop_stop();
        
        profiler.push_start("DCT INVERSE");
#ifndef _SP_COMP_
        fftw_execute(bwdy);
#else // _SP_COMP_
        fftwf_execute(bwdy);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        profiler.push_start("FFTW2CUB");
        _fftw2cubY(datay, grid, gsize[0], gsize[1]);
        profiler.pop_stop();
    }
    
public:
    
    PoissonSolverScalarFFTW_DCT(const int desired_threads, TGrid& grid): FFTWBase(desired_threads), initialized(false)
    {
        // dimensions of blocks and whole grid
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("SETUP");
        _setup(datax, datay, gsize[0], gsize[1]);
        profiler.pop_stop();
    }
    
    void solve(TGrid& grid, const bool spectral=false)
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                    b(ix,iy).tmp = 0;
        }
        
        _FFTx(grid);
        _DCTy(grid);
        /*
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                    b(ix,iy).divU += b(ix,iy).tmp;
        }
        
        //profiler.printSummary();
    }
    
    void dispose()
    {
        if (initialized)
        {
            initialized = false;
            
#ifndef _SP_COMP_
            fftw_destroy_plan(fwdx);
            fftw_destroy_plan(bwdx);
            
            fftw_free(datax);
            
            fftw_destroy_plan(fwdy);
            fftw_destroy_plan(bwdy);
            
            fftw_free(datay);
#else // _SP_COMP_
            fftwf_destroy_plan(fwdx);
            fftwf_destroy_plan(bwdx);
            
            fftwf_free(datax);
            
            fftwf_destroy_plan(fwdy);
            fftwf_destroy_plan(bwdy);
            
            fftwf_free(datay);
#endif // _SP_COMP_
            FFTWBase::dispose();
        }
    }
};
/*/
/*
// 2D - separated 1D FFT, 1D DCT
template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_DCT : FFTWBase
{
    Profiler profiler;
    
protected:
    typedef typename TGrid::BlockType BlockType;
    
    bool initialized;
    
    myplan fwdx, fwdy, bwdx, bwdy;
    Real * datax, * datay; // rhs in _setup, out in cub2fftw and fftw2cub
    
protected:
    
    virtual void _setup(Real *& rhsx, Real *& rhsy, const size_t nx, const size_t ny)
    {
        if (!initialized)
        {
            initialized = true;
            
            //  Forward X:  FFT
            //  Backward X: FFT
            //  Forward Y:  DST-III
            //  Backward Y: DST-II
#ifndef _SP_COMP_
            abort();
#else
            rhsx = fftwf_alloc_real(2*(nx/2+1)); // comes out of r2c/c2r (Section 2.3 of fftw3.pdf)
            
            fwdx = fftwf_plan_dft_r2c_1d(nx, rhsx, (mycomplex *)rhsx, FFTW_MEASURE);
            bwdx = fftwf_plan_dft_c2r_1d(nx, (mycomplex *)rhsx, rhsx, FFTW_MEASURE);
            
            rhsy = fftwf_alloc_real(ny);
            
            fwdy = fftwf_plan_r2r_1d(ny, rhsy, rhsy, FFTW_REDFT01, FFTW_MEASURE);
            bwdy = fftwf_plan_r2r_1d(ny, rhsy, rhsy, FFTW_REDFT10, FFTW_MEASURE);
            //fwdy = fftwf_plan_r2r_1d(ny, rhsy, rhsy, FFTW_REDFT00, FFTW_MEASURE);
            //bwdy = fftwf_plan_r2r_1d(ny, rhsy, rhsy, FFTW_REDFT00, FFTW_MEASURE);
#endif // _SP_COMP_
        }
    }
    
    void _cub2fftwX(TGrid& grid, Real * out, const size_t nx_hat, const int by, const int ly) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            if (info.index[1]==by)
            {
                // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
                const size_t offset = BlockType::sizeX*info.index[0];
                
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index = TStreamer::channels*(offset + ix);
                    
                    assert(dest_index>=0 && dest_index<nx_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(ix,ly), &out[dest_index]);
                }
            }
        }
    }
    
    void _cub2fftwY(TGrid& grid, Real * out, const size_t ny, const int bx, const int lx)
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            if (info.index[0]==bx)
            {
                // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
                const size_t offset = BlockType::sizeY*info.index[1];
                
                for(int iy=0; iy<BlockType::sizeY; iy++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index = TStreamer::channels*(offset + iy);
                    
                    assert(dest_index>=0 && dest_index<ny); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(lx,iy), &out[dest_index]);
                }
            }
        }
 
        
        // nothing like this should be needed for 0-dirichlet
        const int size = ny;
        const int by = grid.getBlocksPerDimension(1)*32;
        double px = (bx+lx)/(double)size;
        double py = 1;
        double x = 4*px*M_PI;
        double y = 3*py*M_PI_2;
        const double offsetY = -9*M_PI_2*M_PI_2 * cos(y) + -64*M_PI_2*M_PI_2 * sin(x);
        //const double offsetY = out[ny-1]; // this value needs to be the one at the boundary! i.e. out[ny]
        //double offsetY = 0;
        
//#pragma omp parallel for
//        for (int j=0; j<ny; j++)
//            offsetY += out[j];
//        offsetY /= nY;
//#pragma omp parallel for
//        for (int j=0; j<ny; j++)
//            out[j] -= offsetY;
 
    }
    
    virtual void _solveX(mycomplex * in_out, const size_t nx, const size_t nx_hat, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int i=0; i<nx_hat; ++i)
        {
            const int linidx = i;
            
            // based on the 5 point stencil in 1D (h^4 error)
            const Real denom = 32.*cos(2.*M_PI*i/nx) - 2.*cos(4.*M_PI*i/nx) - 30.;
            const Real inv_denom = (denom==0)? 0.:1./denom;
            const Real fatfactor = 12. * inv_denom * factor;
            
            //const Real denom = 2.*cos(2.*M_PI*i/nx) - 2.;
            //const Real inv_denom = (denom==0)? 0.:1./denom;
            //const Real fatfactor = inv_denom * factor;
            
            // this is to check the transform only
            //const Real fatfactor = norm_factor;
            
            in_out[linidx][0] *= fatfactor;
            in_out[linidx][1] *= fatfactor;
        }
        
        //this is sparta!
        in_out[0][0] = in_out[0][1] = 0;
    }
    
    virtual void _solveY(Real * in_out, const size_t ny, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int j = 0; j<ny; ++j)
        {
            const int linidx = j;
            
            // based on the 5 point stencil in 1D (h^4 error)
            const Real denom = - 2.*cos(2*M_PI*(j+.5)/ny) + 32.*cos(M_PI*(j+.5)/ny) - 30.;
            const Real inv_denom = (denom==0)? 0.:1./denom;
            const Real fatfactor = 12. * inv_denom * factor;
        
            //const Real denom = 2.*cos(M_PI*(j+.5)/(ny-1)) - 2.;
            //const Real inv_denom = (denom==0)? 0.:1./denom;
            //const Real fatfactor = inv_denom * factor;
            
            // this is to check the transform only
            //const Real fatfactor = norm_factor;
            in_out[linidx] *= fatfactor;
        }
    }
    
    void _fftw2cubX(Real * out, TGrid& grid, const size_t nx_hat, const int by, const int ly) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            if (info.index[1]==by)
            {
                const size_t offset = BlockType::sizeX*info.index[0];
                
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = TStreamer::channels*(offset + ix);
                    
                    assert(src_index>=0 && src_index<nx_hat*2); // linking error with openmp
                    TStreamer::setTmp(&out[src_index], b(ix,ly));
                }
            }
        }
    }
    
    void _fftw2cubY(Real * out, TGrid& grid, const size_t ny, const int bx, const int lx) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            if (info.index[0]==bx)
            {
                const size_t offset = BlockType::sizeY*info.index[1];
                
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        const size_t src_index = TStreamer::channels*(offset+iy);
                        
                        assert(src_index>=0 && src_index<ny); // linking error with openmp
                        TStreamer::operate(&out[src_index], b(lx,iy));
                    }
            }
        }
    }
    
    void _FFTx(TGrid& grid)
    {
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        for (int by=0; by<grid.getBlocksPerDimension(1); by++)
            for (int ly=0; ly<BlockType::sizeY; ly++)
            {
                profiler.push_start("CUB2FFTW");
                _cub2fftwX(grid, datax, gsize[0]/2+1, by, ly);
                profiler.pop_stop();
                
                profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
                fftw_execute(fwdx);
#else // _SP_COMP_
                fftwf_execute(fwdx);
#endif // _SP_COMP_
                profiler.pop_stop();
                
                const Real norm_factor = 1./gsize[0];
                const Real h = grid.getBlocksInfo().front().h_gridpoint;
                assert(1./gsize[0]==h);
                
                profiler.push_start("SOLVE");
                _solveX((mycomplex *)datax, gsize[0], gsize[0]/2+1, norm_factor, h);
                profiler.pop_stop();
                
                profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
                fftw_execute(bwdx);
#else // _SP_COMP_
                fftwf_execute(bwdx);
#endif // _SP_COMP_
                profiler.pop_stop();
                
                profiler.push_start("FFTW2CUB");
                _fftw2cubX(datax, grid, gsize[0]/2+1, by, ly);
                profiler.pop_stop();
            }
    }
    
    void _DCTy(TGrid& grid)
    {
        // dimensions of blocks and whole grid
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        for (int bx=0; bx<grid.getBlocksPerDimension(0); bx++)
            for (int lx=0; lx<BlockType::sizeX; lx++)
            {
                profiler.push_start("CUB2FFTW");
                _cub2fftwY(grid, datay, gsize[1], bx, lx);
                profiler.pop_stop();
                
                profiler.push_start("DCT FORWARD");
#ifndef _SP_COMP_
                fftw_execute(fwdy);
#else // _SP_COMP_
                fftwf_execute(fwdy);
#endif // _SP_COMP_
                profiler.pop_stop();
                
                const Real norm_factor = .5/gsize[1];
                const Real h = grid.getBlocksInfo().front().h_gridpoint;
                assert(1./gsize[1]==h);
                
                profiler.push_start("SOLVE");
                _solveY(datay, gsize[1], norm_factor, h);
                profiler.pop_stop();
 
                profiler.push_start("DCT INVERSE");
#ifndef _SP_COMP_
                fftw_execute(bwdy);
#else // _SP_COMP_
                fftwf_execute(bwdy);
#endif // _SP_COMP_
                profiler.pop_stop();
 
                profiler.push_start("FFTW2CUB");
                _fftw2cubY(datay, grid, gsize[1], bx, lx);
                profiler.pop_stop();
            }
    }
    
public:
    
    PoissonSolverScalarFFTW_DCT(const int desired_threads, TGrid& grid): FFTWBase(desired_threads), initialized(false)
    {
        // dimensions of blocks and whole grid
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("SETUP");
        _setup(datax, datay, gsize[0], gsize[1]);
        profiler.pop_stop();
    }
    
    void solve(TGrid& grid, const bool spectral=false)
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                    b(ix,iy).tmp = 0;
        }
        
        _FFTx(grid);
        _DCTy(grid);
 
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                    b(ix,iy).divU += b(ix,iy).tmp;
        }
        
        //profiler.printSummary();
    }
    
    void dispose()
    {
        if (initialized)
        {
            initialized = false;
            
#ifndef _SP_COMP_
            fftw_destroy_plan(fwdx);
            fftw_destroy_plan(bwdx);
            
            fftw_free(datax);
            
            fftw_destroy_plan(fwdy);
            fftw_destroy_plan(bwdy);
            
            fftw_free(datay);
#else // _SP_COMP_
            fftwf_destroy_plan(fwdx);
            fftwf_destroy_plan(bwdx);
            
            fftwf_free(datax);
            
            fftwf_destroy_plan(fwdy);
            fftwf_destroy_plan(bwdy);
            
            fftwf_free(datay);
#endif // _SP_COMP_
            FFTWBase::dispose();
        }
    }
};
//*/


/*/
 // 1D
 template<typename TGrid, typename TStreamer>
 class PoissonSolverScalarFFTW_DCT : FFTWBase
 {
 Profiler profiler;
 
 protected:
 typedef typename TGrid::BlockType BlockType;
 
 bool initialized;
 
 myplan fwd, bwd;
 Real * data; // rhs in _setup, out in cub2fftw and fftw2cub
 
 protected:
 
 virtual void _setup(Real *& rhs, const size_t ny)
 {
 if (!initialized)
 {
 initialized = true;
 
 //  Forward X:  FFT
 //  Backward X: FFT
 //  Forward Y:  DST-III
 //  Backward Y: DST-II
 #ifndef _SP_COMP_
 rhs = fftw_alloc_real(ny); // why does it need more space than nx*ny?
 abort();
 // mirroring + DFT should accomplish the same
 fwd = fftw_plan_r2r_1d(ny, rhs, rhs, FFTW_REDFT01, FFTW_MEASURE);
 bwd = fftw_plan_r2r_1d(ny, rhs, rhs, FFTW_REDFT10, FFTW_MEASURE); // is the inverse correct? taken from http://www.manpagez.com/info/fftw3/fftw3-3.1.2/fftw3_35.php#SEC42
 #else // _SP_COMP_
 #ifdef _DCT_
 rhs = fftwf_alloc_real(ny);
 
 fwd = fftwf_plan_r2r_1d(ny, rhs, rhs, FFTW_REDFT01, FFTW_MEASURE);
 bwd = fftwf_plan_r2r_1d(ny, rhs, rhs, FFTW_REDFT10, FFTW_MEASURE);
 #else
 rhs = fftwf_alloc_real(4*(ny-1));
 
 fwd = fftwf_plan_r2r_1d(4*(ny-1), rhs, rhs, FFTW_R2HC, FFTW_MEASURE);
 bwd = fftwf_plan_r2r_1d(4*(ny-1), rhs, rhs, FFTW_HC2R, FFTW_MEASURE);
 #endif
 #endif // _SP_COMP_
 }
 }
 
 void _cub2fftw(TGrid& grid, Real * out, const size_t ny) const
 {
 vector<BlockInfo> infos = grid.getBlocksInfo();
 
 const size_t N = infos.size();
 
 #pragma omp parallel for
 for(int i=0; i<N; ++i)
 {
 const BlockInfo info = infos[i];
 BlockType& b = *(BlockType*)infos[i].ptrBlock;
 
 // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
 const size_t offset = BlockType::sizeY*info.index[1];
 
 for(int iy=0; iy<BlockType::sizeY; iy++)
 for(int ix=0; ix<BlockType::sizeX; ix++)
 {
 // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
 const size_t dest_index1 = TStreamer::channels*(           (offset + iy));
 #ifndef _DCT_
 const size_t dest_index2 = TStreamer::channels*(2*(ny-1) - (offset + iy));
 const size_t dest_index3 = TStreamer::channels*(2*(ny-1) + (offset + iy));
 const size_t dest_index4 = TStreamer::channels*(4*(ny-1) - (offset + iy));
 #endif
 
 //assert(dest_index>=0 && dest_index<ny); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
 TStreamer::operate(b(ix,iy), &out[dest_index1]);
 #ifndef _DCT_
 TStreamer::operate(b(ix,iy), &out[dest_index2]);
 TStreamer::operate(b(ix,iy), &out[dest_index3]);
 if (dest_index4 < 4*(ny-1))
 TStreamer::operate(b(ix,iy), &out[dest_index4]);
 #endif
 }
 }
 }
 
 virtual void _solve(Real * in_out, const size_t ny, const Real norm_factor, const Real h)
 {
 if (TStreamer::channels != 1)
 {
 cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
 abort();
 }
 
 const Real h2 = h*h;
 const Real factor = h2*norm_factor;
 
 #pragma omp parallel for
 #ifdef _DCT_
 for(int j = 0; j<ny; ++j)
 #else
 for(int j = 0; j<4*(ny-1); ++j)
 #endif
 {
 const int linidx = j;
 #ifdef _DCT_
 assert(linidx >=0 && linidx<ny); // linking error with openmp
 const Real denom = 32*cos(M_PI*(j+.5)/ny) - 2*cos(2*M_PI*(j+.5)/ny) - 30;
 const Real inv_denom = (denom==0)? 0.:1./denom;
 const Real fatfactor = 12. * inv_denom * factor;
 #else
 assert(linidx >=0 && linidx<4*(ny-1)); // linking error with openmp
 const Real denom = 32.*cos(2.*M_PI*j/ny) - 2.*cos(4.*M_PI*j/ny) - 30.;
 const Real inv_denom = (denom==0)? 0.:1./denom;
 const Real fatfactor = 12. * inv_denom * factor;
 #endif
 // based on the 5 point stencil in 1D (h^4 error)
 
 // based on the 3 point stencil in 1D (h^2 error)
 //const Real denom = 2.*(cos(2.*M_PI*j/ny)) - 2.;
 //const Real inv_denom = (denom==0)? 0.:1./denom;
 //const Real fatfactor = inv_denom * factor;
 
 // this is to check the transform only
 //const Real fatfactor = norm_factor;
 
 in_out[linidx] *= fatfactor;
 }
 
 //this is sparta!
 //in_out[0] = 0;
 }
 
 virtual void _solveSpectral(Real * in_out, const size_t ny, const Real norm_factor, const Real h)
 {
 abort();
 if (TStreamer::channels != 1)
 {
 cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
 abort();
 }
 
 const Real waveFactY = 2.0*M_PI/(ny*h);
 #pragma omp parallel for
 for(int j = 0; j<ny*4-4; ++j)
 {
 const int linidx = j;
 assert(linidx >=0 && linidx<ny); // linking error with openmp
 
 // wave number
 const int ky = (j <= ny/2) ? j : -(ny-j);
 const Real rky = ky*waveFactY;
 const Real kinv = (ky==0) ? 0.0 : -1.0/(rky*rky);
 in_out[linidx] *= kinv*norm_factor;
 }
 
 //this is sparta!
 //in_out[0] = 0;
 }
 
 void _fftw2cub(Real * out, TGrid& grid, const size_t ny) const
 {
 vector<BlockInfo> infos = grid.getBlocksInfo();
 
 const size_t N = infos.size();
 
 const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
 
 #pragma omp parallel for
 for(int i=0; i<N; ++i)
 {
 const BlockInfo info = infos[i];
 BlockType& b = *(BlockType*)infos[i].ptrBlock;
 
 const size_t offset = BlockType::sizeY*info.index[1];
 
 for(int iy=0; iy<BlockType::sizeY; iy++)
 for(int ix=0; ix<BlockType::sizeX; ix++)
 {
 const size_t src_index = TStreamer::channels*(offset + iy);
 
 assert(src_index>=0 && src_index<ny); // linking error with openmp
 TStreamer::operate(&out[src_index], b(ix,iy));
 }
 }
 }
 
 public:
 
 PoissonSolverScalarFFTW_DCT(const int desired_threads): FFTWBase(desired_threads), initialized(false) {	}
 
 void solve(TGrid& grid, const bool spectral=false)
 {
 // dimensions of blocks and whole grid
 const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
 const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
 
 // allocate memory and prepare plans - in an optimized step, can this not be moved out of solve and into the constructor?
 profiler.push_start("SETUP");
 _setup(data, gsize[1]);
 profiler.pop_stop();
 
 profiler.push_start("CUB2FFTW");
 _cub2fftw(grid, data, gsize[1]);
 profiler.pop_stop();
 
 profiler.push_start("FFTW FORWARD");
 #ifndef _SP_COMP_
 fftw_execute(fwd);
 #else // _SP_COMP_
 fftwf_execute(fwd);
 #endif // _SP_COMP_
 profiler.pop_stop();
 
 #ifdef _DCT_
 const Real norm_factor = .5/(gsize[1]);
 #else
 const Real norm_factor = .25/(gsize[1]);
 #endif
 //const Real norm_factor = .5/gsize[0];//1./(gsize[0]*gsize[1]);
 const Real h = grid.getBlocksInfo().front().h_gridpoint;
 assert(1./gsize[0]==h);
 
 profiler.push_start("SOLVE");
 if(spectral)
 _solveSpectral(data, gsize[1], norm_factor, h);
 else
 _solve(data, gsize[1], norm_factor, h);
 
 profiler.pop_stop();
 
 profiler.push_start("FFTW INVERSE");
 #ifndef _SP_COMP_
 fftw_execute(bwd);
 #else // _SP_COMP_
 fftwf_execute(bwd);
 #endif // _SP_COMP_
 profiler.pop_stop();
 
 profiler.push_start("FFTW2CUB");
 _fftw2cub(data, grid, gsize[1]);
 profiler.pop_stop();
 
 profiler.printSummary();
 }
 
 void dispose()
 {
 if (initialized)
 {
 initialized = false;
 
 #ifndef _SP_COMP_
 fftw_destroy_plan(fwd);
 fftw_destroy_plan(bwd);
 
 fftw_free(data);
 #else // _SP_COMP_
 fftwf_destroy_plan(fwd);
 fftwf_destroy_plan(bwd);
 
 fftwf_free(data);
 #endif // _SP_COMP_
 FFTWBase::dispose();
 }
 }
 };
 //*/

