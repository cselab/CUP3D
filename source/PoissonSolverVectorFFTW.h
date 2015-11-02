#pragma once

#include "PoissonSolverScalarFFTW.h"

template<typename TGrid, typename TStreamer>
class PoissonSolverVectorFFTW : public PoissonSolverScalarFFTW< TGrid, TStreamer >
{
protected:
	void _setup(Real *& rhs, const size_t nx, const size_t ny, const size_t nz)
	{
	  if (!(this->initialized))
		{
			this->initialized = true;
			
			const int dim[3] = {nx, ny ,nz}; // need int instead of ptrdiff for nonMPI

#ifndef _SP_COMP_
            rhs = fftw_alloc_real(TStreamer::channels*2*nx*ny*(nz/2+1));
            // first entry is rank = dimension in space --> hardcoded to 3
			this->fwd = fftw_plan_many_dft_r2c(3, dim, TStreamer::channels, rhs, NULL, TStreamer::channels, 1, (mycomplex *)rhs, NULL, TStreamer::channels, 1, FFTW_MEASURE); 
			this->bwd = fftw_plan_many_dft_c2r(3, dim, TStreamer::channels, (mycomplex *)rhs, NULL, TStreamer::channels, 1, rhs, NULL, TStreamer::channels, 1, FFTW_MEASURE);
#else // _SP_COMP_
            rhs = fftwf_alloc_real(TStreamer::channels*2*nx*ny*(nz/2+1));
			this->fwd = fftwf_plan_many_dft_r2c(3, dim, TStreamer::channels, rhs, NULL, TStreamer::channels, 1, (mycomplex *)rhs, NULL, TStreamer::channels, 1, FFTW_MEASURE); 
			this->bwd = fftwf_plan_many_dft_c2r(3, dim, TStreamer::channels, (mycomplex *)rhs, NULL, TStreamer::channels, 1, rhs, NULL, TStreamer::channels, 1, FFTW_MEASURE);
#endif // _SP_COMP_
		}
	}
	
	void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const Real norm_factor, const Real h)
	{
	  if (TStreamer::channels == 1)
	    {
	      //cout << "PoissonSolverVectorFFTW::PoissonSolverVectorFFTW(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be greater than 1. Aborting\n";
	      //abort();
         // do nothing, this also works    
	    }

		const Real h2 = h*h;
		const Real factor = h2*norm_factor; 
#pragma omp parallel for
		for(int i=0; i<nx; ++i)        
            for(int j = 0; j<ny; ++j)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (i*ny+j)*nz_hat + k );
					
					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp
					
					const Real denom = 32.*(cos(2.*M_PI*i/nx)+cos(2.*M_PI*j/ny)+cos(2.0*M_PI*k/nz))-2.*(cos(4.*M_PI*i/nx)+cos(4.*M_PI*j/ny)+cos(4.*M_PI*k/nz))-90.;
                    
					const Real inv_denom = (denom==0)? 0.:1./denom;
					const Real fatfactor = 12. * inv_denom * factor;
					
					for(int c=0; c<TStreamer::channels; ++c)
					{
						in_out[linidx+c][0] *= fatfactor;
						in_out[linidx+c][1] *= fatfactor;	
					}
				}
        
        //this is sparta!
        for(int c=0; c<TStreamer::channels; ++c)
          in_out[c][0] = in_out[c][1] = 0;
	}

    void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const Real norm_factor, const Real h)
	{
        if (TStreamer::channels == 1)
	    {
            //cout << "PoissonSolverVectorFFTW::PoissonSolverVectorFFTW(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be greater than 1. Aborting\n";
            //abort();
            // do nothing, this also works    
	    }
		
        const Real waveFactX = 2.0*M_PI/(nx*h);
        const Real waveFactY = 2.0*M_PI/(ny*h);
        const Real waveFactZ = 2.0*M_PI/(nz*h);
#pragma omp parallel for
		for(int i=0; i<nx; ++i)        
            for(int j = 0; j<ny; ++j)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (i*ny+j)*nz_hat + k );					
					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp
                    
					// wave number
                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (j <= ny/2) ? j : -(ny-j);
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const Real rkx = kx*waveFactX;
                    const Real rky = ky*waveFactY;
                    const Real rkz = kz*waveFactZ;                    
                    const Real kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky+rkz*rkz);
					for(int c=0; c<TStreamer::channels; ++c)
					{
						in_out[linidx+c][0] *= kinv*norm_factor;
						in_out[linidx+c][1] *= kinv*norm_factor;	
					}                    
				}
		
        //this is sparta!
        for(int c=0; c<TStreamer::channels; ++c)
            in_out[c][0] = in_out[c][1] = 0;
	}
    
    void _solveReproject(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const Real norm_factor, const Real h)
	{
        if (TStreamer::channels != 3)
	    {
            cout << "PoissonSolverVectorFFTW::PoissonSolverVectorFFTW(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 3. Aborting\n";
            abort();
	    }
		
        const Real waveFactX = 2.0*M_PI/(nx*h);
        const Real waveFactY = 2.0*M_PI/(ny*h);
        const Real waveFactZ = 2.0*M_PI/(nz*h);
#pragma omp parallel for
		for(int i=0; i<nx; ++i)        
            for(int j = 0; j<ny; ++j)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (i*ny+j)*nz_hat + k );					
					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp
                    
					// wave number
                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (j <= ny/2) ? j : -(ny-j);
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const Real rkx = kx*waveFactX;
                    const Real rky = ky*waveFactY;
                    const Real rkz = kz*waveFactZ;                    
                    const Real kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : 1.0/(rkx*rkx+rky*rky+rkz*rkz);
                    const Real wR[3] = {in_out[linidx][0],in_out[linidx+1][0],in_out[linidx+2][0]};
                    const Real wC[3] = {in_out[linidx][1],in_out[linidx+1][1],in_out[linidx+2][1]};
                    const Real wdotkR = (wR[0]*rkx + wR[1]*rky + wR[2]*rkz)*kinv;
                    const Real wdotkC = (wC[0]*rkx + wC[1]*rky + wC[2]*rkz)*kinv;
                    const Real rk[3] = {rkx,rky,rkz};
					for(int c=0; c<TStreamer::channels; ++c)
					{
						in_out[linidx+c][0] = (wR[c] - wdotkR*rk[c])*norm_factor;
						in_out[linidx+c][1] = (wC[c] - wdotkC*rk[c])*norm_factor;
					}                    
				}
		
        //this is sparta!
        for(int c=0; c<TStreamer::channels; ++c)
            in_out[c][0] = in_out[c][1] = 0;
	}    

    void _solveVelocity(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const Real norm_factor, const Real h)
	{
        if (TStreamer::channels != 3)
	    {
            cout << "PoissonSolverVectorFFTW::PoissonSolverVectorFFTW(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 3. Aborting\n";
            abort();
	    }
		
        const Real waveFactX = 2.0*M_PI/(nx*h);
        const Real waveFactY = 2.0*M_PI/(ny*h);
        const Real waveFactZ = 2.0*M_PI/(nz*h);
#pragma omp parallel for
		for(int i=0; i<nx; ++i)        
            for(int j = 0; j<ny; ++j)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (i*ny+j)*nz_hat + k );					
					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp
                    
					// wave number
                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (j <= ny/2) ? j : -(ny-j);
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const Real rkx = kx*waveFactX;
                    const Real rky = ky*waveFactY;
                    const Real rkz = kz*waveFactZ;                    
                    const Real kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : 1.0/(rkx*rkx+rky*rky+rkz*rkz);
                    const Real wR[3] = {in_out[linidx][0],in_out[linidx+1][0],in_out[linidx+2][0]};
                    const Real wC[3] = {in_out[linidx][1],in_out[linidx+1][1],in_out[linidx+2][1]};

                    in_out[linidx+0][0] = (-rky*wC[2] + rkz*wC[1])*kinv*norm_factor;
                    in_out[linidx+0][1] = ( rky*wR[2] - rkz*wR[1])*kinv*norm_factor;
                    in_out[linidx+1][0] = (-rkz*wC[0] + rkx*wC[2])*kinv*norm_factor;
                    in_out[linidx+1][1] = ( rkz*wR[0] - rkx*wR[2])*kinv*norm_factor;
                    in_out[linidx+2][0] = (-rkx*wC[1] + rky*wC[0])*kinv*norm_factor;
                    in_out[linidx+2][1] = ( rkx*wR[1] - rky*wR[0])*kinv*norm_factor;
				}
		
        //this is sparta!
        for(int c=0; c<TStreamer::channels; ++c)
            in_out[c][0] = in_out[c][1] = 0;
	}
    
public:	
 PoissonSolverVectorFFTW(const int desired_threads) : PoissonSolverScalarFFTW< TGrid, TStreamer >(desired_threads) 
{
  this->initialized = false; 
}
    
    void solveVelocity(TGrid& grid) 
	{
		const int bs[3] = {TGrid::BlockType::sizeX, TGrid::BlockType::sizeY, TGrid::BlockType::sizeZ};
		const size_t gsize[3] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1], grid.getBlocksPerDimension(2)*bs[2]};		
        
		_setup(this->data, gsize[0], gsize[1], gsize[2]);
        
		this->_cub2fftw(grid, this->data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
		
#ifndef _SP_COMP_
		fftw_execute(this->fwd);
#else // _SP_COMP_
		fftwf_execute(this->fwd);
#endif // _SP_COMP_
		
		const Real norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
		const Real h = grid.getBlocksInfo().front().h_gridpoint;
        
        _solveVelocity((mycomplex *)this->data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);        
        
#ifndef _SP_COMP_
		fftw_execute(this->bwd);
#else // _SP_COMP_
		fftwf_execute(this->bwd);
#endif // _SP_COMP_
		
		this->_fftw2cub(this->data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
        
		//profiler.printSummary();
	} 
    
    void reproject(TGrid& grid) 
	{
		const int bs[3] = {TGrid::BlockType::sizeX, TGrid::BlockType::sizeY, TGrid::BlockType::sizeZ};
		const size_t gsize[3] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1], grid.getBlocksPerDimension(2)*bs[2]};		
        
		_setup(this->data, gsize[0], gsize[1], gsize[2]);
        
		this->_cub2fftw(grid, this->data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
		
#ifndef _SP_COMP_
		fftw_execute(this->fwd);
#else // _SP_COMP_
		fftwf_execute(this->fwd);
#endif // _SP_COMP_
		
		const Real norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
		const Real h = grid.getBlocksInfo().front().h_gridpoint;
        
        _solveReproject((mycomplex *)this->data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);        
        
#ifndef _SP_COMP_
		fftw_execute(this->bwd);
#else // _SP_COMP_
		fftwf_execute(this->bwd);
#endif // _SP_COMP_
		
		this->_fftw2cub(this->data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
        
		//profiler.printSummary();
	}     
};
