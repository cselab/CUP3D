/*
 *  BoundaryConditions.h
 *  MPCFnode
 *
 *  Created by Christian Conti on 11/2/15.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "BlockLab.h"
#include <Matrix3D.h>

template<typename TBlock, typename TElement, template<typename X> class allocator=std::allocator>
class BoundaryCondition
{
protected:
	
	int s[3], e[3];
	int stencilStart[3], stencilEnd[3];
	Matrix3D<TElement, true, allocator> * cacheBlock;
	
	// this guy includes corners!
	// what happens if they are not allocated?
	template<int dir, int side>
	void _setup()
	{
		s[0] =	dir==0 ? (side==0 ? stencilStart[0]: TBlock::sizeX) : stencilStart[0];
		s[1] =	dir==1 ? (side==0 ? stencilStart[1]: TBlock::sizeY) : stencilStart[1];
		s[2] =	dir==2 ? (side==0 ? stencilStart[2]: TBlock::sizeZ) : stencilStart[2];
		
		e[0] =	dir==0 ? (side==0 ? 0: TBlock::sizeX + stencilEnd[0]-1) : TBlock::sizeX +  stencilEnd[0]-1;
		e[1] =	dir==1 ? (side==0 ? 0: TBlock::sizeY + stencilEnd[1]-1) : TBlock::sizeY +  stencilEnd[1]-1;
		e[2] =	dir==2 ? (side==0 ? 0: TBlock::sizeZ + stencilEnd[2]-1) : TBlock::sizeZ +  stencilEnd[2]-1;
	}
public:
	
	BoundaryCondition(const int ss[3], const int se[3], Matrix3D<TElement, true, allocator> * cacheBlock):
	cacheBlock(cacheBlock)
	{
		s[0]=s[1]=s[2]=0;
		e[0]=e[1]=e[2]=0;
		
		stencilStart[0] = ss[0];
		stencilStart[1] = ss[1];
		stencilStart[2] = ss[2];
		
		stencilEnd[0] = se[0];
		stencilEnd[1] = se[1];
		stencilEnd[2] = se[2];
	}
	
	TElement& operator()(int ix, int iy, int iz)
	{
		return cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
	}
	
	template<int dir, int side>
	void applyBC_mixedBottom(const TElement& p)
	{
		assert(dir==1);
        assert(side==0);
		
		_setup<dir,side>();
		
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
        		{
                	//(*this)(ix,iy,iz).rho  = p.rho;
                	//(*this)(ix,iy,iz).tmp  = p.tmp;
                	(*this)(ix,iy,iz).chi  = 0;
                
                	// dirichlet BC
					(*this)(ix,iy,iz).u = 2*p.u - (*this)(ix, -iy-1,iz).u;
					(*this)(ix,iy,iz).v = 2*p.v - (*this)(ix, -iy-1,iz).v;
					(*this)(ix,iy,iz).w = 2*p.w - (*this)(ix, -iy-1,iz).w;
					(*this)(ix,iy,iz).tmpU = 2*p.tmpU - (*this)(ix, -iy-1,iz).tmpU;
					(*this)(ix,iy,iz).tmpV = 2*p.tmpV - (*this)(ix, -iy-1,iz).tmpV;
					(*this)(ix,iy,iz).tmpW = 2*p.tmpW - (*this)(ix, -iy-1,iz).tmpW;
					
                	// Neumann BC
                	(*this)(ix,iy,iz).p    = (*this)(ix, -iy-1,iz).p;
                	(*this)(ix,iy,iz).pOld = (*this)(ix, -iy-1,iz).pOld;
                	//(*this)(ix,iy,iz).divU = (*this)(ix, -iy-1,iz).divU;
				}
	}
	
	template<int dir, int side>
	void applyBC_mixedTop(const TElement& p)
	{
		assert(dir==1);
        assert(side==1);
		
		_setup<dir,side>();
		
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
           		{
                	//(*this)(ix,iy,iz).rho  = (*this)(ix, TBlock::sizeY-1,iz).rho;//p.rho;
               		//(*this)(ix,iy,iz).tmp  = (*this)(ix, TBlock::sizeY-1,iz).tmp;//p.rho;
                	(*this)(ix,iy,iz).chi  = 0;
                
                	// dirichlet BC
					(*this)(ix,iy,iz).p    = 2*p.p    - (*this)(ix,2*TBlock::sizeY-1-iy,iz).p;
                	(*this)(ix,iy,iz).pOld = 2*p.pOld - (*this)(ix,2*TBlock::sizeY-1-iy,iz).pOld;
                	//(*this)(ix,iy,iz).divU = 2*p.divU - (*this)(ix,2*TBlock::sizeY-1-iy,iz).divU; // needed because we compute gradP on this!
                
                	// Neumann BC
					(*this)(ix,iy,iz).u = (*this)(ix, 2*TBlock::sizeY-1-iy,iz).u;
					(*this)(ix,iy,iz).v = (*this)(ix, 2*TBlock::sizeY-1-iy,iz).v;
					(*this)(ix,iy,iz).w = (*this)(ix, 2*TBlock::sizeY-1-iy,iz).w;
					(*this)(ix,iy,iz).tmpU = (*this)(ix, 2*TBlock::sizeY-1-iy,iz).tmpU;
					(*this)(ix,iy,iz).tmpV = (*this)(ix, 2*TBlock::sizeY-1-iy,iz).tmpV;
					(*this)(ix,iy,iz).tmpW = (*this)(ix, 2*TBlock::sizeY-1-iy,iz).tmpW;
				}
	}
	
	template<int dir, int side>
	void applyBC_vortex(const BlockInfo info)
	{
		_setup<dir,side>();
		
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
				// what to put here?
				Real p[3];
				info.pos(p, ix, iy,iz);
				
				p[0] = p[0]*2.-1.;
				p[1] = p[1]*2.-1.;
                
                const Real r = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
				const Real invR = 1./r;
                
                //(*this)(ix,iy,iz).rho = r;
				(*this)(ix,iy,iz).u   =   sin(p[1])*cos(r*M_PI/2)*invR;//-p[1];//
				(*this)(ix,iy,iz).v   =  -sin(p[0])*cos(r*M_PI/2)*invR;// p[0];//
				(*this)(ix,iy,iz).w   =  -sin(p[0])*cos(r*M_PI/2)*invR;// p[0];//
				(*this)(ix,iy,iz).chi = 0;
				}
	}
};
