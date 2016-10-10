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
#include <HDF5Dumper_MPI.h>
	
void Simulation::_ic()
{
    CoordinatorIC coordIC(grid);
    profiler.push_start(coordIC.getName());
    coordIC(0);
	profiler.pop_stop();
}
	
void Simulation::setupGrid()
{
	parser.set_strict_mode();
	bpdx = parser("-bpdx").asInt();
	bpdy = parser("-bpdy").asInt();
	bpdz = parser("-bpdz").asInt();
	nprocsx = parser("-nprocsx").asInt();
	parser.unset_strict_mode();
	nprocsy = parser("-nprocsy").asInt(1);
	nprocsz = parser("-nprocsz").asInt(1);
	assert(bpdx%nprocsx==0 && bpdy%nprocsy==0 && bpdz%nprocsz==0);
	if(rank==0) {
		printf("Blocks per dimension: [%d %d %d]\n",bpdx,bpdy,bpdz);
		printf("Nranks per dimension: [%d %d %d]\n",nprocsx,nprocsy,nprocsz);
		printf("Resident blocks per dim: [%d %d %d]\n",bpdx/nprocsx,bpdy/nprocsy,bpdz/nprocsz);
	}
	bpdx /= nprocsx;
	bpdy /= nprocsy;
	bpdz /= nprocsz;
	grid = new FluidGridMPI(nprocsx, nprocsy, nprocsz, bpdx, bpdy, bpdz);
	assert(grid != NULL);
   vInfo = grid->getBlocksInfo();
}

void Simulation::parseArguments()
{
    nu = parser("-nu").asDouble();
    assert(nu>=0);
    parser.unset_strict_mode();
    bRestart = parser("-restart").asBool(false);
    b2Ddump = parser("-2Ddump").asDouble(true);
    bDLM = parser("-use-dlm").asBool(false);
    dumpFreq = parser("-fdump").asDouble(0);	// dumpFreq==0 means that this dumping frequency (in #steps) is not active
    dumpTime = parser("-tdump").asDouble(0.05);	// dumpTime==0 means that this dumping frequency (in time)   is not active
    saveFreq = parser("-fsave").asDouble(0);	// dumpFreq==0 means that this dumping frequency (in #steps) is not active
    saveTime = parser("-tsave").asDouble(10.0);	// dumpTime==0 means that this dumping frequency (in time)   is not active
    nsteps = parser("-nsteps").asInt(0);		// nsteps==0   means that this stopping criteria is not active
    endTime = parser("-tend").asDouble(0);		// endTime==0  means that this stopping criteria is not active

    path2file = parser("-file").asString("./paternoster");
    path4serialization = parser("-serialization").asString("./");
    lambda = parser("-lambda").asDouble(.25);
    CFL = parser("-CFL").asDouble(.25);
    LCFL = parser("-LCFL").asDouble(.1);
    uinf[0] = parser("-uinfx").asDouble(0.0);
    uinf[1] = parser("-uinfy").asDouble(0.0);
    uinf[2] = parser("-uinfz").asDouble(0.0);
    length = parser("-length").asDouble(0.0);
    verbose = parser("-verbose").asBool(false);
    
    //parser.save_options();
}
    
void Simulation::setupObstacles()
{
    IF3D_ObstacleFactory obstacleFactory(grid, uinf);
    obstacle_vector = new IF3D_ObstacleVector(grid, obstacleFactory.create(parser));
    parser.unset_strict_mode();

    //const int nObst = obstacle_vector->nObstacles();
    //_D.resize(nObst);
    //obstacle_vector->_getData(_D);

    Real maxU = max(uinf[0],max(uinf[1],uinf[2]));
    length = obstacle_vector->getD();
    re = length*max(maxU,length)/nu;
    assert(length>0);

    if(rank==0)
	printf("Fluid kinematic viscosity: %9.9e, Reynolds number: %9.9e (length scale: %9.9e)\n",nu,re,length);
}
    
void Simulation::setupOperators()
{
    pipeline.clear();
    pipeline.push_back(new CoordinatorComputeShape(grid, &obstacle_vector, &step, &time, uinf));
    pipeline.push_back(new CoordinatorAdvection<LabMPI>(uinf, grid));
    //if(nu>0)
    pipeline.push_back(new CoordinatorDiffusion<LabMPI>(nu, grid));
    pipeline.push_back(new CoordinatorPenalization(grid, &obstacle_vector, &lambda, uinf));
    pipeline.push_back(new CoordinatorPressure<LabMPI>(grid, &obstacle_vector));
    pipeline.push_back(new CoordinatorComputeForces(grid, &obstacle_vector, &step, &time, &nu, &bDump, uinf));
    pipeline.push_back(new CoordinatorComputeDiagnostics(grid, &obstacle_vector, &step, &time, &lambda, uinf));
    pipeline.push_back(new CoordinatorFadeOut(grid));
    if(rank==0) {
    	cout << "Coordinator/Operator ordering:\n";
    	for (int c=0; c<pipeline.size(); c++) cout << "\t" << pipeline[c]->getName() << endl;
    }
}

void Simulation::areWeDumping(Real & nextDumpTime)
{
	bDump = (dumpFreq>0 && step%dumpFreq==0) || (dumpTime>0 && time>=nextDumpTime);
	if (bDump) nextDumpTime += dumpTime;
#ifdef _BSMART_
	for (int i=0; i<obstacle_vector->nObstacles(); i++)
	{if (t+DT>=_D[i]->t_next_comm) bDump=true;}
	if (bDump) obstacle_vector->getFieldOfView();
#endif
}

void Simulation::_dump(const string append = string())
{
    //stringstream ss;
    //ss << path2file << append << "-" << std::setfill('0') << std::setw(6) << step;

    stringstream ssR;
    if (append == string())
    	ssR<<path4serialization<<"./restart_"<<std::setfill('0')<<std::setw(9)<<step;
    else
        ssR<<path4serialization<<"./"<<append<<std::setfill('0')<<std::setw(9)<<step;

    if (rank==0) cout << ssR.str() << endl;

    if (rank==0) { //rank 0 saves step id and obstacles
    	obstacle_vector->save(step,time,ssR.str());
    	//safety status in case of crash/timeout during grid save:
    	string statusname = ssR.str()+".status";
    	FILE * f = fopen(statusname.c_str(), "w");
    	assert(f != NULL);
    	fprintf(f, "time: %20.20e\n", time);
    	fprintf(f, "stepid: %d\n", (int)step);
    	//fprintf(f, "ping: %d\n", (int)bPing);
    	fclose(f);
    }

#if defined(_USE_HDF_)
    if(b2Ddump) {
      stringstream ssF;
      if (append == string())
         ssF<<path4serialization<<"./avemaria_"<<std::setfill('0')<<std::setw(9)<<step;
      else
         ssF<<path4serialization<<"./2D_"<<append<<std::setfill('0')<<std::setw(9)<<step;
    	DumpHDF5flat_MPI(*grid, time, ssF.str());
    }
    DumpHDF5_MPI(*grid, time, ssR.str());
#endif
#if defined(_USE_LZ4_) //TODO: does not compile

    CoordinatorVorticity<LabMPI> coordVorticity(grid);
    coordVorticity(dt);
    MPI_Barrier(MPI_COMM_WORLD);
    Real vpeps = parser("-vpeps").asDouble(1e-5);
    int wavelet_type = parser("-wtype").asInt(1);

    waveletdumper_grid.verbose();
    waveletdumper_grid.set_wtype_write(wavelet_type);
    waveletdumper_grid.set_threshold (vpeps);
    waveletdumper_grid.Write<0>(grid, ss.str());
    waveletdumper_grid.Write<1>(grid, ss.str());
    waveletdumper_grid.Write<2>(grid, ss.str());
    waveletdumper_grid.Write<3>(grid, ss.str());
    waveletdumper_grid.Write<4>(grid, ss.str());
    waveletdumper_grid.Write<5>(grid, ss.str());
    waveletdumper_grid.Write<6>(grid, ss.str());
    waveletdumper_grid.Write<7>(grid, ss.str());
#endif

	if (rank==0) { //saved the grid! Write status to remember most recent ping
		string restart_status = path4serialization+"./restart.status";
		FILE * f = fopen(restart_status.c_str(), "w");
		assert(f != NULL);
		fprintf(f, "time: %20.20e\n", time);
		fprintf(f, "stepid: %d\n", (int)step);
		//fprintf(f, "ping: %d\n", (int)bPing);
		fclose(f);

		printf( "time: %20.20e\n", time);
		printf( "stepid: %d\n", (int)step);
		//printf( "ping: %d\n", (int)bPing);
	}

    CoordinatorDiagnostics coordDiags(grid,time,step);
    coordDiags(dt);
}
    
void Simulation::_selectDT()
{
	double local_maxU = (double)findMaxUOMP(vInfo,*grid,uinf);
	double global_maxU;
   const double h = vInfo[0].h_gridpoint;
    MPI::COMM_WORLD.Allreduce(&local_maxU, &global_maxU, 1, MPI::DOUBLE, MPI::MAX);
    dtFourier = CFL*h*h/nu;
    dtCFL     = CFL*h/abs(global_maxU);
    dt = min(dtCFL,dtFourier);

    if(rank==0) printf("maxU %f dtF %f dtC %f dt %f\n",global_maxU,dtFourier,dtCFL,dt);
    //if (dumpTime>0) dt = min(dt,nextDumpTime-time);
    //if (saveTime>0) dt = min(dt,nextSaveTime-time);
    //if (endTime>0)  dt = min(dt,endTime-time);

    if ( step<100 ) {
        const Real dt_max = 2e-2*CFL;
        const Real dt_min = 2e-4*CFL;
        const Real dt_ramp = dt_min + step*(dt_max - dt_min)/100.;
        if (dt_ramp<dt) {
        	dt = dt_ramp;
        	if(rank==0) printf("Dt bounded by ramp-up: dt_ramp=%f\n",dt_ramp);
        }
    }
}

void Simulation::_serialize(Real & nextSaveTime)
{// FAKE LOLOLOLOLOLOLOLOLOLOLOL (read the dump)
	bool bSaving = (saveFreq>0 && step%saveFreq==0)||(saveTime>0 && time>nextSaveTime);
	if(!bSaving) return;
	nextSaveTime += saveTime;
	//DumpZBin_MPI<FluidGridMPI, StreamerSerialization>(*grid, pingname.c_str(), path4serialization);
}
	
void Simulation::_deserialize()
{
	string restartfile = path4serialization+"./restart.status";
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
	//int ping_fake = -1;
	//fscanf(f, "ping: %d\n", &ping_fake);
	//assert(ping_fake >= 0);
	//bPing = (bool)ping_fake;
	fclose(f);

	stringstream ssR;
	ssR<<path4serialization<<"./restart_"<<std::setfill('0')<<std::setw(9)<<step;
	if (rank==0) cout << "Restarting from " << ssR.str() << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	ReadHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, ssR.str());
	obstacle_vector->restart(time,ssR.str());

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
    
    MPI_Barrier(MPI_COMM_WORLD);
}
	
void Simulation::simulate()
{
    Real nextDumpTime = time;
    Real nextSaveTime = time + saveTime;
    
    while (true)
    {
        profiler.push_start("DT");
        _selectDT();
        areWeDumping(nextDumpTime);
        profiler.pop_stop();

#ifdef _BSMART_
        profiler.push_start("LEARN");
        bool bDoOver = false;
        const int nO = obstacle_vector->nObstacles();
        for(int i=1; i<nO; i++) {
        	bDoOver = _D[i]->checkFail(_D[0]->Xrel, _D[0]->Yrel, _D[0]->thExp, length);
        	if (bDoOver) {
        		if (i>0) _D[i]->finalizePos(_D[0]->Xrel,  _D[0]->Yrel,  _D[0]->thExp, _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1.0);
        		_D[i]->info = 2;
        		obstacle_vector->execute(comm, i, t);
        		return;
        	}
        }
        for(int i=0; i<nO; i++) if(t>=_D[i]->t_next_comm) {
        	if (i>0) _D[i]->finalize(_D[0]->Xrel,  _D[0]->Yrel,  _D[0]->thExp, _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1.0);
        	obstacle_vector->execute(comm, i, t);
        }
        if (bDoOver) exit(0);
        profiler.pop_stop();
#endif
        if(bDLM) lambda = 1.0/dt;

        for (int c=0; c<pipeline.size(); c++)
        {
        	//cout << pipeline[c]->getName() << " start" << endl;
            profiler.push_start(pipeline[c]->getName());
            (*pipeline[c])(dt);
            profiler.pop_stop();
            //if(time>.0025)  _dump(pipeline[c]->getName());
        }
        //if(time>.0025) abort();
        step++;
        time += dt;
        if(rank==0)
        printf("Step %d time %f\n",step,time);

        profiler.push_start("Dump");
        if(bDump)
        {
            _dump();
        }
        profiler.pop_stop();

        profiler.push_start("Save");
        {
        	_serialize(nextSaveTime);
        }
        profiler.pop_stop();
        
        if (step % 100 == 0) profiler.printSummary();
        if ((endTime>0 && time>endTime) || (nsteps!=0 && step>=nsteps))
        {
        	cout<<"Finished at time "<<time<<" in "<<step<<" step of "<<nsteps<<endl;
			exit(0);
        }
    }
}
