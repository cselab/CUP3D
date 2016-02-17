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
	string shapeType;
	Shape * shape;
	
	
	void _outputSettings(ostream& outStream)
	{
		outStream << "Simulation_FSI\n";
		outStream << "lambda " << lambda << endl;
		outStream << "dlm " << dlm << endl;
		outStream << "shapeType " << shapeType << endl;
		
		shape->serialize(outStream);
		
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
		assert(variableName=="shapeType");
		inStream >> shapeType;
		
		if (shapeType=="samara")
			shape = GeometryMesh::deserialize(inStream);
		else
		{
			cout << "Deserialization for this shape not implemented yet\n";
			abort();
		}
		
		/*
		inStream >> variableName;
		Real center[3];
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
			assert(variableName=="rhoS");
			inStream >> rhoS;
			inStream >> variableName;
			assert(variableName=="mollChi");
			inStream >> mollChi;
			inStream >> variableName;
			assert(variableName=="mollRho");
			inStream >> mollRho;
			shape = new Sphere(center, radius, rhoS, mollChi, mollRho);
		}
		else
		{
			cout << "Error - this shape is not currently implemented! Aborting now\n";
			abort();
		}
		*/
		
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
			dlm = parser("-dlm").asDouble(1.);
			
			double rhoS = parser("-rhoS").asDouble(1);
			const Real aspectRatio = 1;
			if (rank==0) cout << "WARNING - Aspect ratio for correct positioning of sphere not implemented yet\n";
			Real center[3] = {parser("-xpos").asDouble(.5*aspectRatio),parser("-ypos").asDouble(.85),parser("-zpos").asDouble(.5*aspectRatio)};
			
			string shapeType = parser("-shape").asString("sphere");
            const int eps = 2;
			if (shapeType=="sphere")
			{
				Real radius = parser("-radius").asDouble(0.1);
				shape = new Sphere(center, radius, rhoS, eps, eps);
			}
			else if (shapeType=="samara" || shapeType=="triangle")
			{
				/*
				const Real center[3] = {.5,.5,.5};
				const Real rhoS = 1;
				const Real moll = 2;
				const int gridsize = 1024;
				const Real scale = .12;
				const Real tx = .07;
				const Real ty = .18;
				const Real tz = .08;
				const Real isosurface = parser("-isosurface").asDouble(.004);
				Geometry::Quaternion q1(cos(.5*M_PI), 0, 0, sin(.5*M_PI));
				//Geometry::Quaternion q2(1, 0, 0, 0);
				Geometry::Quaternion q2(cos(45./360.*M_PI), sin(45./360.*M_PI), 0, 0);
				Geometry::Quaternion q = q1*q2;
				
				double qm = q.magnitude();
				q.w /= qm;
				q.x /= qm;
				q.y /= qm;
				q.z /= qm;
				const string filename = "/cluster/home/infk/cconti/CubismUP_3D/launch/geometries/Samara_v3.obj";
				shape = new GeometryMesh(filename, gridsize, isosurface, center, rhoS, moll, moll, scale, tx, ty, tz, q);
				 */
			}
			else
			{
				cout << "Error - this shape is not currently implemented! Aborting now\n";
				abort();
			}
		}
		else
		{
			// shape stuff
			cout << "Restart - Simulation_FSI\n";
			//abort();
		}
		
		// nothing needs to be done on restart
	}
	
	virtual ~Simulation_FSI()
	{
		delete shape;
	}
};

#endif
