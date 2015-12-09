//
//  Simulation_Fluid.h
//  CubismUP_3D
//
//	Base class for fluid simulations from which any fluid simulation case should inherit
//	Contains the base structure and interface that any fluid simulation class should have
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_Simulation_Fluid_h
#define CubismUP_3D_Simulation_Fluid_h

#include "Definitions.h"
#include "ProcessOperatorsOMP.h"
#include "GenericCoordinator.h"
#include "GenericOperator.h"

#include "CoordinatorVorticity.h"

//#include "SerializerIO_WaveletCompression.h"

#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

class Simulation_Fluid
{
protected:
	ArgumentParser parser;
	Profiler profiler;
	
	// Serialization
	bool bPing; // needed for ping-pong scheme
	string path4serialization;
	bool bRestart;
	
	// MPI stuff - required for Hypre
	int rank, nprocs;
	
	vector<GenericCoordinator *> pipeline;
	
	// grid
	int bpdx, bpdy, bpdz;
	FluidGrid * grid;
	
	// simulation status
	int step, nsteps;
	double dt, time, endTime;
	
	// simulation settings
	double CFL, LCFL;
	
	// verbose
	bool verbose;
	
	// output
	int dumpFreq;
	double dumpTime;
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
//	SerializerIO_WaveletCompression_SimpleBlocking<FluidGrid, FluidVPStreamer> waveletdumper_grid;
//	SerializerIO_WaveletCompression_SimpleBlocking<FluidGrid, TmpVPStreamer> waveletdumper_vorticity;
	
	virtual void _diagnostics() = 0;
	virtual void _ic() = 0;
	virtual double _nonDimensionalTime() = 0;
	
	virtual void _dump(double & nextDumpTime)
	{
#ifndef NDEBUG
		if (rank==0)
		{
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			const int N = vInfo.size();
			
#pragma omp parallel for schedule(static)
			for(int i=0; i<N; i++)
			{
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				
				for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						if (std::isnan(b(ix,iy,iz).rho) ||
							std::isnan(b(ix,iy,iz).u) ||
							std::isnan(b(ix,iy,iz).v) ||
							std::isnan(b(ix,iy,iz).chi) ||
							std::isnan(b(ix,iy,iz).p) ||
							std::isnan(b(ix,iy,iz).pOld))
							cout << "dump" << endl;
						
						if (b(ix,iy,iz).rho <= 0)
							cout << "dump " << b(ix,iy,iz).rho << "\t" << info.index[0] << " " << info.index[1] << " " << info.index[2] << " " << ix << " " << iy << " " << iz << endl;
						
						assert(b(ix,iy,iz).rho > 0);
						assert(!std::isnan(b(ix,iy,iz).rho));
						assert(!std::isnan(b(ix,iy,iz).u));
						assert(!std::isnan(b(ix,iy,iz).v));
						assert(!std::isnan(b(ix,iy,iz).chi));
						assert(!std::isnan(b(ix,iy,iz).p));
						assert(!std::isnan(b(ix,iy,iz).pOld));
						assert(!std::isnan(b(ix,iy,iz).tmpU));
						assert(!std::isnan(b(ix,iy,iz).tmpV));
						assert(!std::isnan(b(ix,iy,iz).tmp));
						assert(!std::isnan(b(ix,iy,iz).divU));
					}
			}
		}
#endif
		
		const int sizeX = bpdx * FluidBlock::sizeX;
		const int sizeY = bpdy * FluidBlock::sizeY;
		const int sizeZ = bpdz * FluidBlock::sizeZ;
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		if(rank==0 && (dumpFreq>0 && step % dumpFreq == 0) || (dumpTime>0 && abs(nextDumpTime-_nonDimensionalTime()) < 10*std::numeric_limits<Real>::epsilon()))
		{
			nextDumpTime += dumpTime;
			
#ifndef _USE_HDF_
			stringstream ss;
			ss << path2file << "-" << step << ".vti";
			cout << ss.str() << endl;
			
			dumper.Write(*grid, ss.str());
#else
			CoordinatorVorticity<Lab> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << "-" << std::setfill('0') << std::setw(6) << step;
			cout << ss.str() << endl;
			DumpHDF5<FluidGrid, StreamerHDF5>(*grid, step, ss.str());
#endif
			_serialize();
			
			/*
			// VP
			std::stringstream streamer;
			streamer << path2file << "-datawavelet";
			streamer.setf(ios::dec | ios::right);
			streamer.width(6);
			streamer.fill('0');
			streamer << step;
				
			double vpeps = parser("-vpeps").asDouble(1e-3);
			int wavelet_type = parser("-wtype").asInt(1);
				
			waveletdumper_grid.verbose();
			waveletdumper_grid.set_wtype_write(wavelet_type);
				
			waveletdumper_grid.set_threshold (vpeps);
			if (vpchannels.find('0') != std::string::npos)
				waveletdumper_grid.Write<0>(grid, streamer.str());
				
			waveletdumper_grid.set_threshold (vpeps);
			if (vpchannels.find('1') != std::string::npos)
				waveletdumper_grid.Write<1>(grid, streamer.str());
				
			waveletdumper_grid.set_threshold (vpeps);
			if (vpchannels.find('2') != std::string::npos)
				waveletdumper_grid.Write<2>(grid, streamer.str());
				
			waveletdumper_grid.set_threshold (vpeps);
			if (vpchannels.find('3') != std::string::npos)
				waveletdumper_grid.Write<3>(grid, streamer.str());
				
			waveletdumper_grid.set_threshold (vpeps);
			if (vpchannels.find('4') != std::string::npos)
				waveletdumper_grid.Write<4>(grid, streamer.str());
				
			waveletdumper_grid.set_threshold (vpeps);
			if (vpchannels.find('5') != std::string::npos)
				waveletdumper_grid.Write<5>(grid, streamer.str());
				
			waveletdumper_vorticity.verbose();
			waveletdumper_vorticity.set_wtype_write(wavelet_type);
			waveletdumper_vorticity.set_threshold (vpeps);
			if (vpchannels.find('w') != std::string::npos || vpchannels.find('W') != std::string::npos)
				waveletdumper_vorticity.Write<5>(grid, streamer.str());
				
			waveletdumper_sos.verbose();
			waveletdumper_sos.set_wtype_write(wavelet_type);
			waveletdumper_sos.set_threshold (10*vpeps);
			if (vpchannels.find('c') != std::string::npos)
				waveletdumper_sos.Write<0>(grid, streamer.str());
				
			waveletdumper_mach.verbose();
			waveletdumper_mach.set_wtype_write(wavelet_type);
			waveletdumper_mach.set_threshold (vpeps);
			if (vpchannels.find('M') != std::string::npos)
				waveletdumper_mach.Write<0>(grid, streamer.str());
				
			waveletdumper_velocity_magnitude.verbose();
			waveletdumper_velocity_magnitude.set_wtype_write(wavelet_type);
			waveletdumper_velocity_magnitude.set_threshold (vpeps);
			if (vpchannels.find('m') != std::string::npos)
				waveletdumper_velocity_magnitude.Write<0>(grid, streamer.str());
			 */
		}
	}
	
	virtual void _outputSettings(ostream& outStream)
	{
		outStream << "Simulation_Fluid\n";
		
		outStream << "step " << step << endl;
		outStream << "nsteps " << nsteps << endl;
		outStream << "time " << time << endl;
		outStream << "endTime " << endTime << endl;
		
		outStream << "verbose " << verbose << endl;
		
		outStream << "CFL " << CFL << endl;
		outStream << "LCFL " << LCFL << endl;
		
		outStream << "dumpFreq " << dumpFreq << endl;
		outStream << "dumpTime " << dumpTime << endl;
		outStream << "path2file " << path2file << endl;
		
		outStream << "Grid " << bpdx << " " << bpdy << " " << bpdz << endl;
	}
	
	virtual void _inputSettings(istream& inStream)
	{
		string variableName;
		
		inStream >> variableName;
		if (variableName != "Simulation_Fluid")
		{
			cout << "Error in deserialization - Simulation_Fluid\n";
			abort();
		}
		
		// read data
		inStream >> variableName;
		assert(variableName=="step");
		inStream >> step;
		inStream >> variableName;
		assert(variableName=="nsteps");
		inStream >> nsteps;
		inStream >> variableName;
		assert(variableName=="time");
		inStream >> time;
		inStream >> variableName;
		assert(variableName=="endTime");
		inStream >> endTime;
		inStream >> variableName;
		assert(variableName=="verbose");
		inStream >> verbose;
		inStream >> variableName;
		assert(variableName=="CFL");
		inStream >> CFL;
		inStream >> variableName;
		assert(variableName=="LCFL");
		inStream >> LCFL;
		inStream >> variableName;
		assert(variableName=="dumpFreq");
		inStream >> dumpFreq;
		inStream >> variableName;
		assert(variableName=="dumpTime");
		inStream >> dumpTime;
		inStream >> variableName;
		assert(variableName=="path2file");
		inStream >> path2file;
		inStream >> variableName;
		assert(variableName=="Grid");
		inStream >> bpdx;
		inStream >> bpdy;
		inStream >> bpdz;
	}
	
	void _serialize()
	{
		if (rank==0)
		{
			stringstream ss;
			ss << path4serialization << "Serialized-" << bPing << ".dat";
			cout << ss.str() << endl;
			
			stringstream serializedGrid;
			serializedGrid << "SerializedGrid-" << bPing << ".grid";
			DumpZBin<FluidGrid, StreamerSerialization>(*grid, serializedGrid.str(), path4serialization);
			
			ofstream file;
			file.open(ss.str());
			
			if (file.is_open())
			{
				_outputSettings(file);
				
				file.close();
			}
			
			bPing = !bPing;
		}
	}
	
	void _deserialize()
	{
		stringstream ss0, ss1, ss;
		struct stat st0, st1;
		ss0 << path4serialization << "Serialized-0.dat";
		ss1 << path4serialization << "Serialized-1.dat";
		stat(ss0.str().c_str(), &st0);
		stat(ss1.str().c_str(), &st1);
		
		
		// direct comparison of the two quantities leads to segfault
		bPing = st0.st_mtime<st1.st_mtime ? false : true;
		ss << (!bPing ? ss0.str() : ss1.str());
		
		ifstream file;
		file.open(ss.str());
		
		if (file.is_open())
		{
			_inputSettings(file);
			
			file.close();
		}
		
		grid = new FluidGrid(bpdx,bpdy,bpdz);
		assert(grid != NULL);
		
		if (rank==0)
		{
			stringstream serializedGrid;
			serializedGrid << "SerializedGrid-" << bPing << ".grid";
			ReadZBin<FluidGrid, StreamerSerialization>(*grid, serializedGrid.str(), path4serialization);
		}
	}
	
public:
	Simulation_Fluid(const int argc, const char ** argv) : parser(argc,argv), step(0), time(0), dt(0), rank(0), nprocs(1), bPing(false)
	{
	}
	
	virtual ~Simulation_Fluid()
	{
		delete grid;
		
		while(!pipeline.empty())
		{
			GenericCoordinator * g = pipeline.back();
			pipeline.pop_back();
			delete g;
		}
	}
	
	virtual void init()
	{
		bRestart = parser("-restart").asBool(false);
		if (rank==0)
			cout << "bRestart is " << bRestart << endl;
		
		if (!bRestart)
		{
			// initialize grid
			parser.set_strict_mode();
			bpdx = parser("-bpdx").asInt();
			bpdy = parser("-bpdy").asInt();
			bpdz = parser("-bpdz").asInt();
			grid = new FluidGrid(bpdx,bpdy,bpdz);
			assert(grid != NULL);
			
			// simulation ending parameters
			parser.unset_strict_mode();
			nsteps = parser("-nsteps").asInt(0);		// nsteps==0   means that this stopping criteria is not active
			endTime = parser("-tend").asDouble(0);		// endTime==0  means that this stopping criteria is not active
			
			// output parameters
			dumpFreq = parser("-fdump").asDouble(0);	// dumpFreq==0 means that this dumping frequency (in #steps) is not active
			dumpTime = parser("-tdump").asDouble(0);	// dumpTime==0 means that this dumping frequency (in time)   is not active
			path2file = parser("-file").asString("../data/Simulation_Fluid");
			path4serialization = parser("-serialization").asString(path2file);
			
			CFL = parser("-CFL").asDouble(.25);
			LCFL = parser("-LCFL").asDouble(.1);
			
			verbose = parser("-verbose").asBool(false);
		}
		else
		{
			if (rank==0)
				cout << "Deserializing...";
			
			parser.set_strict_mode();
			path4serialization = parser("-serialization").asString();
			parser.unset_strict_mode();
			_deserialize();
			
			// evenutally read new endTime, nsteps
			nsteps = parser("-nsteps").asInt(nsteps);
			endTime = parser("tend").asDouble(endTime);
			
			if (rank==0)
			{
				cout << " done - parameters:\n";
				_outputSettings(cout);
			}
			
			if (rank==0)
			{
				double d = _nonDimensionalTime();
				_dump(d);
			}
			
#ifdef _MULTIGRID_
			MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
		}
	}
	
	virtual void simulate() = 0;
};

#endif
