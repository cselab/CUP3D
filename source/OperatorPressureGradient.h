//
//  OperatorPressureGradient.h
//  CubismUP_3D
//
//	Add a gradient of pressure, from gradient to 0
//
//  Created by Christian Conti on 5/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorPressureGradient_h
#define CubismUP_3D_OperatorPressureGradient_h

class OperatorPressureGradient : public GenericOperator
{
private:
	double dt;
	Real gradient[3];
	
public:
	OperatorPressureGradient(Real gradient[3], double dt) : dt(dt), gradient{gradient[0],gradient[1],gradient[2]}
	{
	}
	
	~OperatorPressureGradient() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				// assumes [0,1[^2 domain
				
				double p[3];
				info.pos(p,ix,iy,iz);
				
				block(ix,iy,iz).u += dt*gradient[0]*(1-p[0]);
				block(ix,iy,iz).v += dt*gradient[1]*(1-p[1]);
				block(ix,iy,iz).w += dt*gradient[2]*(1-p[2]);
			}
	}
};


#endif
