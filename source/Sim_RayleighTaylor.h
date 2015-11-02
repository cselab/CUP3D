//
//  Sim_RayleighTaylor.h
//  CubismUP_3D
//
//	Class for the simulation of Rayleigh-Taylor Instability
//		used to validate against Herrmann
//			fixed timestep
//			fixed nondimentionalization of time
//
//	Makefile arguments
//		particles:		if active, use particles (remeshing kernel should be chosen by commenting/uncommenting in the code), otherwise use finite differences
//		bc:				boundary conditions, for this setting, mixed should be used
//		poisson:		Poisson solver, for this setting, hypre should be used
//		multiphase:		activate/deactivate advection of densities, for this setting, true should be used
//		precision:		single/double, default: double
//		config:			configuration debug/production, default: debug
//		nthreads:		#threads, default: 48
//		bs:				Cubism block size, default: 32
//		vertexcentered:	always set to false for this setting
//
//	Command line arguments
//		-bpdx:			mandatory, #blocks in x direction
//		-bpdy:			mandatory, #blocks in y direction
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
//		-split:			use split technique for the pressure solver (only used by Hypre implementation)
//		-nu:			mu for diffusion operator (which is divided by the density to recover nu)
//		-rhoS:			heavy phase density
//
//  Created by Christian Conti on 4/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__Sim_RayleighTaylor__
#define __CubismUP_3D__Sim_RayleighTaylor__

#include "Simulation_MP.h"

class Sim_RayleighTaylor : public Simulation_MP
{
protected:
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
public:
	Sim_RayleighTaylor(const int argc, const char ** argv);
	~Sim_RayleighTaylor();
	
	void init();
	void simulate();
};

#endif /* defined(__CubismUP_3D__Sim_RayleighTaylor__) */
