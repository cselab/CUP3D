//
//  Sim_FSI_Rotation.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/10/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef Sim_FSI_Rotation_h
#define Sim_FSI_Rotation_h

#include "Simulation_FSI.h"

class Sim_FSI_Rotation : public Simulation_FSI
{
protected:
	Real omegaBody[3], uBody[3];
	double dtBody, dtCFL, dtLCFL, dtFourier;
	double re, nu;
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
public:
	Sim_FSI_Rotation(const int argc, const char ** argv);
	virtual ~Sim_FSI_Rotation();
	
	void init();
	void simulate();
};

#endif /* Sim_FSI_Rotation_h */
