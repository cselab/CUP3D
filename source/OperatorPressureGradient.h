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
	Real gradient[2];
	
public:
	OperatorPressureGradient(Real gradient[2], double dt) : dt(dt), gradient{gradient[0],gradient[1]}
	{
	}
	
	~OperatorPressureGradient() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				// assumes [0,1[^2 domain
				
				double p[2];
				info.pos(p,ix,iy);
				
				block(ix,iy).u += dt*gradient[0]*(1-p[0]);
				block(ix,iy).v += dt*gradient[1]*(1-p[1]);
			}
	}
};


#endif
