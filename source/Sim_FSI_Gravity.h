//
//  Sim_FSI_Gravity.h
//  CubismUP_3D
//
//	Class for the simulation of gravity driven FSI
//
//	Makefile arguments
//		particles:		if active, use particles (remeshing kernel should be chosen by commenting/uncommenting in the code), otherwise use finite differences
//		poisson:		Poisson solver, split-fftw/hypre
//		multiphase:		activate/deactivate advection of densities, false should be used
//		bc:				boundary conditions, for this setting, mixed should be used
//		precision:		single/double, default: double
//		config:			configuration debug/production, default: debug
//		nthreads:		#threads, default: 48
//		bs:				Cubism block size, default: 32
//		vertexcentered:	always set to false for this setting
//
//	Command line arguments
//		-bpdx:			mandatory, #blocks in x direction
//		-bpdy:			mandatory, #blocks in y direction
//		-bpdz:			mandatory, #blocks in y direction
//		-restart:		needed to restart from a previous checkpoint
//		-nsteps:		simulation ends with nsteps iterations. If 0, ignored
//		-tend:			simulation ends at time tend. If 0, ignored
//		-fdump:			#timesteps between grid outputs. If 0, ignored
//		-tdump:			time interval for grid output. If 0, ignored
//		-file:			location and base name of output files
//		-CFL:			CFL constant used for CFL condition and for diffusion
//		-LCFL:			particle Lagrangian-CFL (not required when using finite differences)
//		-verbose:		activate verbosity
//		-serialization:	mandatory if restarting, path to serialized checkpoint data
//		-lambda:		Penalization parameter, default: 1e5
//		-rhoS:			solid density
//		-shape:			select shape to be used: disk/ellipse, default: disk
//		-radius:		if shape==disk, select radius, default: 0.1
//		-semiAxisX:		if shape==ellipse, select X axis, default: 0.1
//		-semiAxisY:		if shape==ellipse, select Y axis, default: 0.2
//		-semiAxisZ:		if shape==ellipse, select Z axis, default: 0.05
//		-angle:			if shape==ellipse, select orientation
//		-split:			use split technique for the pressure solver (only used by Hypre implementation)
//		-nu:			mu for diffusion operator (which is divided by the density to recover nu)
//		-ypos:			position in the vertical direction used to place the shape in the initial conditions
//
//	Output
//		diagnostics.dat
//			 0 - step
//			 1 - time
//			 2 - dt
//			 3 - bpdx
//			 4 - lambda
//			 5 - cD
//			 6 - Re
//			 7 - x[0]
//			 8 - x[1]
//			 9 - x[2]
//			10 - u[0]
//			11 - u[1]
//			12 - u[2]
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__Sim_FSI_Gravity__
#define __CubismUP_3D__Sim_FSI_Gravity__

#include "Simulation_FSI.h"

class Sim_FSI_Gravity : public Simulation_FSI
{
protected:
	Real uBody[3], uBodyOld[3];
	Real aBody[3];
	double dtCFL, dtLCFL, dtFourier, dtBody;
	Real re, nu;
	Real minRho;
	bool bSplit;
	double maxU;
	
	Real gravity[3];
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
	// should this stuff be moved? - serialize method will do that
	void _dumpSettings(ostream& outStream);
	
public:
	Sim_FSI_Gravity(const int argc, const char ** argv);
	~Sim_FSI_Gravity();
	
	void init();
	void simulate();
};

#endif /* defined(__CubismUP_3D__Sim_FSI_Gravity__) */
