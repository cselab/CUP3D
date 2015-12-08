/*
 *  Layer.h
 *  VM2D
 *
 *  Created by Diego Rossinelli on 2/9/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once 
#include <math.h>
#include <string>
#include <vector>
using namespace std;
#include <assert.h>

#include "common.h"

struct Layer
{
	const int sizeX;
	const int sizeY;
	const int sizeZ;
	const int nDim;

	Real * data;
	
	Layer(const int sizeX, const int sizeY, const int sizeZ, const int nDim) : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), nDim(nDim)
	{
		data = new Real[nDim*sizeX*sizeY*sizeZ];
	}
	
	~Layer()
	{
		delete [] data;
	}

	inline Real& operator()(int ix=0, int iy=0, int iz=0, int dim=0)
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);
		
		return data[dim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
	}

	inline Real read(int ix=0, int iy=0, int iz=0, int dim=0) const
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);
		
		return data[dim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
	}

	const Layer& operator=(const Real val)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iz = 0; iz<sizeZ; iz++)
				for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] = val;
		
		return *this;
	}

	const Layer& operator=(const Layer& l)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iz = 0; iz<sizeZ; iz++)
			for(int iy = 0; iy<sizeY; iy++)
				for(int ix = 0; ix<sizeX; ix++)
					data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] = l.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
		
		return *this;
	}

	template<int dim>
	void clear(Real val)
	{
		for(int iz = 0; iz<sizeZ; iz++)
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[dim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] = val;
	}
	
	double getH0() const
	{
		return 1./(double)sizeX;
	}
	
	double getH1() const
	{
		return 1./(double)sizeY;
	}
	
	double getH2() const
	{
		return 1./(double)sizeZ;
	}

	const vector<double> operator -(const Layer& layer)
	{
		vector<double> result;
		
		//compute linf distance
		{
			double LInf_diff = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iz = 0; iz<sizeZ; iz++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						LInf_diff = max(LInf_diff, (double)fabs(data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] - layer.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix]));
			
			result.push_back(LInf_diff);
		}
		
		//compute linf distance
		{
			double L2error = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iz = 0; iz<sizeZ; iz++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						L2error += pow(data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] - layer.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix], 2);
			
			result.push_back(sqrt((double)L2error/(sizeY*sizeX)));
		}
		
		return result;
	}
	
	
	void difference(const Layer& input, Layer& difference)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iz = 0; iz<sizeZ; iz++)
			for(int iy = 0; iy<sizeY; iy++)
				for(int ix = 0; ix<sizeX; ix++)
					difference.data[idim*sizeX*sizeY + iy*sizeX + ix]= data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] - input.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
	}

	template<int iDim>
	Real * getPlane() 
	{
		return (Real*)&data[iDim*sizeZ*sizeX*sizeY];
	}
};


#include <xmmintrin.h>
template <typename LayerT> 
LayerT * allocate()
{
	void * data = _mm_malloc(sizeof(LayerT), 64);
	return new (data) LayerT;
}


