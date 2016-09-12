//
//  Simulation_Fluid.h
//  CubismUP_2D
//
//	Base class for fluid simulations from which any fluid simulation case should inherit
//	Contains the base structure and interface that any fluid simulation class should have
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Simulation.h"
	
void Simulation::_ic()
{
    CoordinatorIC coordIC(grid);
    profiler.push_start(coordIC.getName());
    coordIC(0);
}
	
void Simulation::setupGrid()
{
	parser.set_strict_mode();
	nprocsx = parser("-nprocsx").asInt();
	nprocsy = parser("-nprocsy").asInt();
	nprocsz = parser("-nprocsz").asInt();
	bpdx = parser("-bpdx").asInt();
	bpdy = parser("-bpdy").asInt();
	bpdz = parser("-bpdz").asInt();

	grid = new FluidGridMPI(nprocsx, nprocsy, nprocsz, bpdx, bpdy, bpdz);
	assert(grid != NULL);

    vInfo = grid->getBlocksInfo();
}

void Simulation::parseArguments()
{
    nu = parser("-nu").asDouble();
    
    parser.unset_strict_mode();
    theta = 0.5;
    bRestart = parser("-restart").asBool(false);
    dumpFreq = parser("-fdump").asDouble(0);	// dumpFreq==0 means that this dumping frequency (in #steps) is not active
    dumpTime = parser("-tdump").asDouble(0.25);	// dumpTime==0 means that this dumping frequency (in time)   is not active
    saveFreq = parser("-fsave").asDouble(0);	// dumpFreq==0 means that this dumping frequency (in #steps) is not active
    saveTime = parser("-tsave").asDouble(1.0);	// dumpTime==0 means that this dumping frequency (in time)   is not active
    nsteps = parser("-nsteps").asInt(0);		// nsteps==0   means that this stopping criteria is not active
    endTime = parser("-tend").asDouble(8);		// endTime==0  means that this stopping criteria is not active
    
    path2file = parser("-file").asString("./paternoster");
    path4serialization = parser("-serialization").asString("./");
    lambda = parser("-lambda").asDouble(.25);
    CFL = parser("-CFL").asDouble(.25);
    LCFL = parser("-LCFL").asDouble(.1);
    uinf[0] = parser("-uinf").asDouble(0.0);
    length = parser("-length").asDouble(0.0);
    printf("Fluid kinematic viscosity: %20.20e (length scale = %20.20e)\n", nu, length);
    verbose = parser("-verbose").asBool(false);
    
    parser.save_options();
}
    
void Simulation::setupObstacles()
{
    IF2D_ObstacleFactory obstacleFactory(grid, uinf);
    obstacle_vector = new IF2D_ObstacleVector(grid, obstacleFactory.create(parser));
    parser.unset_strict_mode();
}
    
void Simulation::setupOperators()
{
    pipeline.clear();
    pipeline.push_back(new CoordinatorComputeShape(grid, &obstacle_vector, &step, &time, uinf));
    pipeline.push_back(new CoordinatorAdvection<Lab>(uinf, grid));
    pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
    pipeline.push_back(new CoordinatorPenalization(grid, &obstacle_vector, &lambda, uinf));
    pipeline.push_back(new CoordinatorPressure<Lab>(grid, &obstacle_vector));

    pipeline.push_back(new CoordinatorComputeForces(grid, &obstacle_vector, &step, &time, &nu, &bDump, uinf));
    pipeline.push_back(new CoordinatorComputeDiagnostics(grid, &obstacle_vector, &step, &time, &lambda, uinf));
    
    cout << "Coordinator/Operator ordering:\n";
    for (int c=0; c<pipeline.size(); c++) cout << "\t" << pipeline[c]->getName() << endl;
}

void Simulation::_dump(double & nextDumpTime)
{
#ifndef NDEBUG
    if (rank==0) {
        vector<BlockInfo> vInfo = grid->getBlocksInfo();
        const int N = vInfo.size();
        
        #pragma omp parallel for schedule(static)
        for(int i=0; i<N; i++) {
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;

			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
            for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
                if (std::isnan(b(ix,iy,iz).u) ||
                    std::isnan(b(ix,iy,iz).v) ||
                    std::isnan(b(ix,iy,iz).w) ||
                    std::isnan(b(ix,iy,iz).chi) ||
                    std::isnan(b(ix,iy,iz).p) ||
                    std::isnan(b(ix,iy,iz).pOld))
                    cout << "dump" << endl;
                assert(!std::isnan(b(ix,iy,iz).u));
                assert(!std::isnan(b(ix,iy,iz).v));
                assert(!std::isnan(b(ix,iy,iz).w));
                assert(!std::isnan(b(ix,iy,iz).chi));
                assert(!std::isnan(b(ix,iy,iz).p));
                assert(!std::isnan(b(ix,iy,iz).pOld));
                assert(!std::isnan(b(ix,iy,iz).tmpU));
                assert(!std::isnan(b(ix,iy,iz).tmpV));
                assert(!std::isnan(b(ix,iy,iz).tmpW));
                assert(!std::isnan(b(ix,iy,iz).tmp));
            }
        }
    }
#endif
    
    const int sizeX = bpdx * FluidBlock::sizeX;
    const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    if((dumpFreq>0 && step%dumpFreq==0)||(dumpTime>0 && time>nextDumpTime)) {
        nextDumpTime += dumpTime;
        
#ifdef _USE_HDF_
			CoordinatorVorticity<LabMPI> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << "-" << std::setfill('0') << std::setw(6) << step;
			if (rank==0) cout << ss.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
			// if (rank==0) _serialize(); TODO: ask
			MPI_Barrier(MPI_COMM_WORLD);
			/* VP
				 std::stringstream streamer;
				 streamer << path2file << "-datawavelet";
				 streamer.setf(ios::dec | ios::right);
				 streamer.width(6);
				 streamer.fill('0');
				 streamer << step;

				 double vpeps = parser("-vpeps").asDouble(1e-5);
				 int wavelet_type = parser("-wtype").asInt(1);

				 waveletdumper_grid.verbose();
				 waveletdumper_grid.set_wtype_write(wavelet_type);

				 waveletdumper_grid.set_threshold (vpeps);
				 waveletdumper_grid.Write<0>(grid, streamer.str());
				 waveletdumper_grid.Write<1>(grid, streamer.str());
				 waveletdumper_grid.Write<2>(grid, streamer.str());
				 waveletdumper_grid.Write<3>(grid, streamer.str());
				 waveletdumper_grid.Write<4>(grid, streamer.str());
				 waveletdumper_grid.Write<5>(grid, streamer.str());
					/*
				 waveletdumper_vorticity.verbose();
				 waveletdumper_vorticity.set_wtype_write(wavelet_type);
				 waveletdumper_vorticity.set_threshold (vpeps);
				 if (vpchannels.find('w') != std::string::npos || vpchannels.find('W') != std::string::npos)
					waveletdumper_vorticity.Write<5>(grid, streamer.str());

				 waveletdumper_velocity_magnitude.verbose();
				 waveletdumper_velocity_magnitude.set_wtype_write(wavelet_type);
				 waveletdumper_velocity_magnitude.set_threshold (vpeps);
				 if (vpchannels.find('m') != std::string::npos)
					waveletdumper_velocity_magnitude.Write<0>(grid, streamer.str());
			*/
    }
}
    
void Simulation::_selectDT()
{
    const double maxU = findMaxUOMP(vInfo,*grid,uinf);
    dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
    dtCFL     = CFL*vInfo[0].h_gridpoint/abs(maxU);
    dt = min(dtCFL,dtFourier);
#ifdef _PARTICLES_
    //const double maxA = findMaxAOMP<Lab>(vInfo,*grid);
    //dtLCFL = maxA==0 ? 1e5 : LCFL/abs(maxA);
    //dt = min(dt,dtLCFL);
#endif
    
    printf("maxU %f dtF %f dtC %f dt %f\n",maxU,dtFourier,dtCFL,dt);
    //if (dumpTime>0) dt = min(dt,nextDumpTime-time);
    //if (saveTime>0) dt = min(dt,nextSaveTime-time);
    //if (endTime>0)  dt = min(dt,endTime-time);
    
    if ( step<1000 ) {
        const double dt_max = 0.01;
        const double dt_min = 1e-6;
        const double dt_ramp = dt_min + step*(dt_max - dt_min)/1000.;
        dt = min(dt,dt_ramp);
        printf("dt_ramp %f dt %f\n",dt_ramp,dt);
    }

    //if (verbose) cout << "dt (Fourier, CFL): " << dtFourier << " " << dtCFL << endl;
}

void Simulation::_serialize(double & nextSaveTime)
{
	/*
	 * MUSINGS TODO:
	 * 	- Original version was ping scheme:
	 * 		advantage is that you only ever have 2 restart files
	 * 	- IF2(3)D uses some restart frequency and keeps a growing list of restarts
	 * 		advantage is that is publication friendly
	 * 		disadvantage is that you save two times per save
	 *
	 * 	I'd use the original version since these files weigh a ton. So TODO: make it work.
	 */
	if((saveFreq>0 && step%saveFreq==0)||(saveTime>0 && time>nextSaveTime)) {
		nextSaveTime += saveTime;
		if (rank==0)
		{ //only rank 0 saves step id and obstacles
			printf("****SERIALIZING****\n");
			/*
			{ //write status
				string restart_status;
				char buf[500];
				sprintf(buf, "/restart.status");
				restart_status = path4serialization+string(buf);

				FILE * f = fopen(restart_status.c_str(), "w");
				assert(f != NULL);
				fprintf(f, "time: %20.20e\n", time);
				fprintf(f, "stepid: %d\n", (int)step);
				fclose(f);

				printf( "time: %20.20e\n", time);
				printf( "stepid: %d\n", (int)step);
			}
			*/
			{ //write numbered status (extra safety measure)
				string numbered_status;
				char buf[500];
				sprintf(buf, "restart_%07d.status", (int)step);
				numbered_status = path4serialization+string(buf);

				FILE * f = fopen(numbered_status.c_str(), "w");
				assert(f != NULL);
				fprintf(f, "time: %20.20e\n", time);
				fprintf(f, "stepid: %d\n", (int)step);
				fclose(f);
			}
			obstacle_vector->save(step,time,path4serialization);
		}

		{
			//DumpZBin<FluidGrid, StreamerSerialization>(*grid, "restart", path4serialization);
		}
		{
			bPing = !bPing;
			string numbered_filename;
			char buf[500];
			sprintf(buf, "restart_%07d", (int)step);
			numbered_filename = string(buf);
			DumpZBin_MPI<FluidGridMPI, StreamerSerialization>(*grid, numbered_filename.c_str(), path4serialization);
		}
		printf("****SERIALIZING DONE****\n");
	}
}
	
void Simulation::_deserialize()
{
	{
		string restartfile = path4serialization+"restart.status";
		FILE * f = fopen(restartfile.c_str(), "r");
		assert(f != NULL);
		float val = -1;
		fscanf(f, "time: %e\n", &val);
		assert(val>=0);
		time=val;
		int step_id_fake = -1;
		fscanf(f, "stepid: %d\n", &step_id_fake);
		assert(step_id_fake >= 0);
		step = step_id_fake;
		fclose(f);
	}
	{
		stringstream ssHDF;
		ssHDF << path4serialization << "restart" << step;
		if (rank==0) cout << ssHDF.str() << endl;
		ReadHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, ssHDF.str());
	}
	{
		string restartfile = path4serialization+"restart";
		obstacle_vector->restart(time,restartfile);
	}

    printf("DESERIALIZATION: time is %f and step id is %d\n", time, (int)step);
}
    
void Simulation::init()
{
    parseArguments();
    setupGrid();
    setupObstacles();
    setupOperators();

    if(bRestart) _deserialize();
    else _ic();
    /*
    if (bRestart)
    {
        parser.set_strict_mode();
        parser.unset_strict_mode();
        _deserialize();
        
        // evenutally read new endTime, nsteps
        nsteps = parser("-nsteps").asInt(nsteps);
        endTime = parser("tend").asDouble(endTime);
        
        if (rank==0)
        {
            double d = _nonDimensionalTime();
            _dump(d);
        }
        
#ifdef _MULTIGRID_
        MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
    }
     */
    
    MPI_Barrier(MPI_COMM_WORLD);
}
	
void Simulation::simulate()
{
    const int sizeX = bpdx * FluidBlock::sizeX;
    const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
    
    double nextDumpTime = time;
    double nextSaveTime = time + saveTime;
    double maxU = max(uinf[2],max(uinf[0],uinf[1]));
    double maxA = 0;
    
    while (true)
    {
        profiler.push_start("DT");
        _selectDT();
        profiler.pop_stop();
        bDump = (dumpTime>0 && time+dt>nextDumpTime);
        
        if(bDLM)
        	lambda = 1.0/dt;

        for (int c=0; c<pipeline.size(); c++) {
            profiler.push_start(pipeline[c]->getName());
            (*pipeline[c])(dt);
            profiler.pop_stop();
        }

        profiler.push_start("Dump");
        {
            _dump(nextDumpTime);
        }
        profiler.pop_stop();


        profiler.push_start("Save");
        {
        	_serialize(nextSaveTime);
        }
        profiler.pop_stop();

        time += dt;
        step++;
        printf("Step %d time %f\n",step,time);
        
        if (step % 1000 == 0) profiler.printSummary();
        if ((endTime>0 && time>endTime) || (nsteps!=0 && step>=nsteps)) {
        	cout<<"Finished at time "<<time<<" (target end time "<<endTime<<") in "<<step<<" step of "<<nsteps<<endl;
			exit(0);
        }
    }
}
