//
//  Sim_Bubble.h
//  CubismUP_3D
//
//  Created by Christian Conti on 4/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__Sim_Bubble__
#define __CubismUP_3D__Sim_Bubble__

#include "Simulation_MP.h"

class Sim_Bubble : public Simulation_MP
{
protected:
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
public:
	Sim_Bubble(const int argc, const char ** argv);
	~Sim_Bubble();
	
	void init();
	void simulate();
};

#endif /* defined(__CubismUP_3D__Sim_Bubble__) */
