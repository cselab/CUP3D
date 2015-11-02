//
//  Sim_FSI_Fixed.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__Sim_FSI_Fixed__
#define __CubismUP_3D__Sim_FSI_Fixed__

#include "Simulation_FSI.h"

class Sim_FSI_Fixed : public Simulation_FSI
{
protected:
    double uinf, re, nu;
	double dtCFL, dtLCFL, dtFourier;
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
public:
	Sim_FSI_Fixed(const int argc, const char ** argv);
	virtual ~Sim_FSI_Fixed();
	
	void init();
    void simulate();
};


#endif /* defined(__CubismUP_3D__Sim_FSI_Fixed__) */
