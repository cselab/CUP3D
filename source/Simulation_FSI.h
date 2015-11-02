//
//  Simulation_FSI.h
//  CubismUP_3D
//
//	Base class for Fluid-Structure Interaction (FSI) simulations from which any FSI simulation case should inherit
//	Contains the base structure and interface that any FSI simulation class should have
//	Inherits from Simulation_Fluid
//	Assumes use of Penalization to handle rigid body
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_Simulation_FSI_h
#define CubismUP_3D_Simulation_FSI_h

#include "Simulation_Fluid.h"
#include "Shape.h"

class Simulation_FSI : public Simulation_Fluid
{
protected:
	// penalization parameter
	Real lambda, dlm;
	
	// body
	Shape * shape;
	
	
	void _outputSettings(ostream& outStream)
	{
		outStream << "Simulation_FSI\n";
		outStream << "lambda " << lambda << endl;
		outStream << "dlm " << dlm << endl;
		shape->outputSettings(outStream);
		
		Simulation_Fluid::_outputSettings(outStream);
	}
	
	virtual void _inputSettings(istream& inStream)
	{
		string variableName;
		
		inStream >> variableName;
		if (variableName != "Simulation_FSI")
		{
			cout << "Error in deserialization - Simulation_FSI\n";
			abort();
		}
		
		// read data
		inStream >> variableName;
		assert(variableName=="lambda");
		inStream >> lambda;
		inStream >> variableName;
		assert(variableName=="dlm");
		inStream >> dlm;
		
		inStream >> variableName;
		Real center[3];
		bool bPeriodic[3] = {false,false,false};
		Real rhoS;
		Real mollChi, mollRho;
		if (variableName=="Sphere")
		{
			Real radius; // need one more angle
			
			inStream >> variableName;
			assert(variableName=="radius");
			inStream >> radius;
			inStream >> variableName;
			assert(variableName=="centerX");
			inStream >> center[0];
			inStream >> variableName;
			assert(variableName=="centerY");
			inStream >> center[1];
			inStream >> variableName;
			assert(variableName=="centerZ");
			inStream >> center[2];
			inStream >> variableName;
			assert(variableName=="orientation");
			inStream >> variableName;
			inStream >> variableName;
			assert(variableName=="rhoS");
			inStream >> rhoS;
			inStream >> variableName;
			assert(variableName=="mollChi");
			inStream >> mollChi;
			inStream >> variableName;
			assert(variableName=="mollRho");
			inStream >> mollRho;
            
            vector<BlockInfo> vInfo = grid->getBlocksInfo();
			const Real domainSize[3] = { FluidBlock::sizeX * grid->getBlocksPerDimension(0) * vInfo[0].h_gridpoint,
										 FluidBlock::sizeY * grid->getBlocksPerDimension(1) * vInfo[0].h_gridpoint,
										 FluidBlock::sizeZ * grid->getBlocksPerDimension(2) * vInfo[0].h_gridpoint};
			shape = new Sphere(center, radius, rhoS, mollChi, mollRho, bPeriodic, domainSize);
		}
		else if (variableName=="Ellipsoid")
		{
			Real semiAxis[3];
			Real angle; // need one more angle
			
			
			inStream >> variableName;
			assert(variableName=="semiAxisX");
			inStream >> semiAxis[0];
			inStream >> variableName;
			assert(variableName=="semiAxisY");
			inStream >> semiAxis[1];
			inStream >> variableName;
			assert(variableName=="semiAxisZ");
			inStream >> semiAxis[2];
			inStream >> variableName;
			assert(variableName=="centerX");
			inStream >> center[0];
			inStream >> variableName;
			assert(variableName=="centerY");
			inStream >> center[1];
			inStream >> variableName;
			assert(variableName=="centerZ");
			inStream >> center[2];
			inStream >> variableName;
			assert(variableName=="orientation");
			inStream >> angle;
			inStream >> variableName;
			assert(variableName=="rhoS");
			inStream >> rhoS;
			inStream >> variableName;
			assert(variableName=="mollChi");
			inStream >> mollChi;
			inStream >> variableName;
			assert(variableName=="mollRho");
            inStream >> mollRho;
            
            vector<BlockInfo> vInfo = grid->getBlocksInfo();
			const Real domainSize[3] = { FluidBlock::sizeX * grid->getBlocksPerDimension(0) * vInfo[0].h_gridpoint,
										 FluidBlock::sizeY * grid->getBlocksPerDimension(1) * vInfo[0].h_gridpoint,
										 FluidBlock::sizeZ * grid->getBlocksPerDimension(2) * vInfo[0].h_gridpoint};
			shape = new Ellipsoid(center, semiAxis, angle, rhoS, mollChi, mollRho, bPeriodic, domainSize);
		}
		else
		{
			cout << "Error - this shape is not currently implemented! Aborting now\n";
			abort();
		}
		
		Simulation_Fluid::_inputSettings(inStream);
	}
	
public:
	Simulation_FSI(const int argc, const char ** argv) : Simulation_Fluid(argc,argv)
	{
	}
	
	virtual void init()
	{
		Simulation_Fluid::init();
		
		if (!bRestart)
		{
			lambda = parser("-lambda").asDouble(1e5);
			dlm = parser("-lambda").asDouble(1.);
			
			double rhoS = parser("-rhoS").asDouble(1);
			Real centerOfMass[3] = {0,0,0};
			bool bPeriodic[3] = {false,false,false};
            
            vector<BlockInfo> vInfo = grid->getBlocksInfo();
            const Real domainSize[3] = { FluidBlock::sizeX * grid->getBlocksPerDimension(0) * vInfo[0].h_gridpoint,
										 FluidBlock::sizeY * grid->getBlocksPerDimension(1) * vInfo[0].h_gridpoint,
										 FluidBlock::sizeZ * grid->getBlocksPerDimension(2) * vInfo[0].h_gridpoint};
			
			string shapeType = parser("-shape").asString("disk");
            const int eps = 2;
			if (shapeType=="disk")
			{
				Real radius = parser("-radius").asDouble(0.1);
				shape = new Sphere(centerOfMass, radius, rhoS, eps, eps, bPeriodic, domainSize);
			}
			else if (shapeType=="ellipse")
			{
				Real semiAxis[3] = {parser("-semiAxisX").asDouble(0.1),parser("-semiAxisY").asDouble(0.2),parser("-semiAxisZ").asDouble(0.05)};
				Real angle = parser("-angle").asDouble(0.0);
				shape = new Ellipsoid(centerOfMass, semiAxis, angle, rhoS, eps, eps, bPeriodic, domainSize);
			}
			else
			{
				cout << "Error - this shape is not currently implemented! Aborting now\n";
				abort();
			}
		}
		
		// nothing needs to be done on restart
	}
	
	virtual ~Simulation_FSI()
	{
		delete shape;
	}
};

#endif
