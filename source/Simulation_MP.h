//
//  Simulation_MP.h
//  CubismUP_3D
//
//  Created by Christian Conti on 4/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__Simulation_MP__
#define __CubismUP_3D__Simulation_MP__

#include "Simulation_Fluid.h"

class Simulation_MP: public Simulation_Fluid
{
protected:
	double dtCFL, dtLCFL, dtFourier;
	double re, nu;
	double minRho, rhoS;
	bool bSplit;
	
	Real gravity[3];
	
	void _outputSettings(ostream& outStream)
	{
		outStream << "Multiphase\n";
		outStream << "re " << re << endl;
		outStream << "nu " << nu << endl;
		outStream << "minRho " << minRho << endl;
		outStream << "rhoS " << rhoS << endl;
		
		Simulation_Fluid::_outputSettings(outStream);
	}
	
	void _inputSettings(istream& inStream)
	{
		string variableName;
		
		inStream >> variableName;
		if (variableName != "Multiphase")
		{
			cout << "Error in deserialization - Simulation_Multiphase\n";
			abort();
		}
		
		// read data
		inStream >> variableName;
		assert(variableName=="re");
		inStream >> re;
		inStream >> variableName;
		assert(variableName=="nu");
		inStream >> nu;
		inStream >> variableName;
		assert(variableName=="minRho");
		inStream >> minRho;
		inStream >> variableName;
		assert(variableName=="rhoS");
		inStream >> rhoS;
		
		Simulation_Fluid::_inputSettings(inStream);
	}
	
public:
	Simulation_MP(const int argc, const char ** argv) : Simulation_Fluid(argc, argv), gravity{0,-9.81,0}, dtCFL(0), dtFourier(0), re(0), nu(0), minRho(0), rhoS(1), bSplit(false)
	{
		
	}
	
	~Simulation_MP()
	{
	}
	
	void init()
	{
		Simulation_Fluid::init();
		
		if (!bRestart)
		{
			// simulation settings
			bSplit = parser("-split").asBool(false);
			nu = parser("-nu").asDouble(1e-2);
			rhoS = parser("-rhoS").asDouble(1);
			minRho = min((Real)1.,(Real)rhoS);
			
			if (rank==0)
				if (bSplit)
					cout << "Using split method with constant coefficients Poisson solver\n";
				else
					cout << "Solving full variable coefficient Poisson equation for pressure\n";
		}
	}
	
	void simulate()
	{
		
	}
};

#endif /* defined(__CubismUP_3D__Sim_Multiphase__) */
