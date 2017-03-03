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
#include "ProcessOperatorsOMP.h"
#include <chrono>

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

	bpdx /= nprocsx;
	bpdy /= nprocsy;
	bpdz /= nprocsz;
	grid = new FluidGridMPI(nprocsx, nprocsy, nprocsz, bpdx, bpdy, bpdz, 1, app_comm);
	dump = new  DumpGridMPI(nprocsx, nprocsy, nprocsz, bpdx, bpdy, bpdz, 1, app_comm);
	assert(grid != NULL);
  vInfo = grid->getBlocksInfo();

  MPI_Comm_rank(grid->getCartComm(),&rank);
  MPI_Comm_size(grid->getCartComm(),&nprocs);
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  const int nthreads = omp_get_max_threads();
  printf("Rank %d (of %d) with %d threads on host Hostname: %s\n",
          rank, nprocs, nthreads, hostname);
  if (communicator not_eq nullptr) //Yo dawg I heard you like communicators.
    communicator->comm_MPI = grid->getCartComm();
  fflush(0);
	if(rank==0) {
		printf("Blocks per dimension: [%d %d %d]\n",bpdx,bpdy,bpdz);
		printf("Nranks per dimension: [%d %d %d]\n",nprocsx,nprocsy,nprocsz);
		printf("Resident blocks per dim: [%d %d %d]\n",
          bpdx/nprocsx,bpdy/nprocsy,bpdz/nprocsz);
	}
}

void Simulation::parseArguments()
{
    nu = parser("-nu").asDouble();
    assert(nu>=0);
    parser.unset_strict_mode();
    bRestart = parser("-restart").asBool(false);
    b2Ddump = parser("-2Ddump").asDouble(true);
    bDLM = parser("-use-dlm").asBool(false);
    dumpFreq = parser("-fdump").asDouble(0);	// dumpFreq==0 means dump freq (in #steps) is not active
    dumpTime = parser("-tdump").asDouble(0.05);	// dumpTime==0 means dump freq (in time)   is not active
    saveFreq = parser("-fsave").asDouble(0);	// dumpFreq==0 means dump freq (in #steps) is not active
    saveTime = parser("-tsave").asDouble(10.0);	// dumpTime==0 means dump freq (in time)   is not active
    nsteps = parser("-nsteps").asInt(0);		// nsteps==0   means stopping criteria is not active
    endTime = parser("-tend").asDouble(0);		// endTime==0  means stopping criteria is not active

    path2file = parser("-file").asString("./paternoster");
    path4serialization = parser("-serialization").asString("./");
    maxClockDuration = parser("-Wtime").asDouble(1e30);
    lambda = parser("-lambda").asDouble(1e4);
    CFL = parser("-CFL").asDouble(.1);
    uinf[0] = parser("-uinfx").asDouble(0.0);
    uinf[1] = parser("-uinfy").asDouble(0.0);
    uinf[2] = parser("-uinfz").asDouble(0.0);
    length = parser("-length").asDouble();
    verbose = parser("-verbose").asBool(false);
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
    if(nu>0)
    	pipeline.push_back(new CoordinatorDiffusion<LabMPI>(nu, grid));
    pipeline.push_back(new CoordinatorPenalization(grid, &obstacle_vector, &lambda, uinf));
    pipeline.push_back(new CoordinatorPressure<LabMPI>(grid, &obstacle_vector));
    pipeline.push_back(new CoordinatorComputeForces(grid, &obstacle_vector, &step, &time, &nu, &bDump, uinf));
    pipeline.push_back(new CoordinatorComputeDiagnostics(grid, &obstacle_vector, &step, &time, &lambda, uinf));
    //#ifndef _OPEN_BC_
    pipeline.push_back(new CoordinatorFadeOut(grid));
    //#endif
    if(rank==0) {
    	cout << "Coordinator/Operator ordering:\n";
    	for (int c=0; c<pipeline.size(); c++) cout << "\t" << pipeline[c]->getName() << endl;
    }
}

void Simulation::areWeDumping(Real & nextDumpTime)
{
	bDump = (dumpFreq>0 && step+1%dumpFreq==0) || (dumpTime>0 && time>=nextDumpTime);
	if (bDump) nextDumpTime += dumpTime;

  #ifdef __SMARTIES_
  //const auto _D = obstacle_vector->_getData();
	//for (int i=0; i<obstacle_vector->nObstacles(); i++)
	//   if (time+dt>=_D[i]->t_next_comm) bDump=true;
	//if (bDump) obstacle_vector->getFieldOfView();
  #endif
}

void Simulation::_dump(const string append = string())
{
    stringstream ssR;
    if (append == string())
    	ssR<<"restart_";
    else
      ssR<<append;
      ssR << std::setfill('0')<<std::setw(9)<<step;
    if (rank==0) cout << ssR.str() << endl;

    if (rank==0) { //rank 0 saves step id and obstacles
    	obstacle_vector->save(step,time,path4serialization+"/"+ssR.str());
    	//safety status in case of crash/timeout during grid save:
    	string statusname = path4serialization+"/"+ssR.str()+".status";
    	FILE * f = fopen(statusname.c_str(), "w");
    	assert(f != NULL);
    	fprintf(f, "time: %20.20e\n", time);
    	fprintf(f, "stepid: %d\n", (int)step);
    	fprintf(f, "uinfx: %20.20e\n", uinf[0]);
    	fprintf(f, "uinfy: %20.20e\n", uinf[1]);
    	fprintf(f, "uinfz: %20.20e\n", uinf[2]);
    	fclose(f);
    }

    #if defined(_USE_HDF_)

      if(b2Ddump) {
      stringstream ssF;
      if (append == string())
       ssF<<"avemaria_"<<std::setfill('0')<<std::setw(9)<<step;
      else
       ssF<<"2D_"<<append<<std::setfill('0')<<std::setw(9)<<step;
      DumpHDF5flat_MPI(*grid, time, ssF.str(),path4serialization);
      }
      /*
      copyDumpGrid(*grid, *dump);
      if(dumper not_eq nullptr){
        bool warned = false;
        while(not dumper->joinable()) {
          MPI_Barrier(grid->getCartComm());
          if (!warned && rank==0) {
            printf("Waiting for previous dump.\n");
            fflush(0);
            warned = true;
          }
        }
        dumper->join();
        if(dumper not_eq nullptr)
          delete dumper;
      }
      //why do they even let me code?
      string string_path = path4serialization;
      string string_file = ssR.str();
    	DumpGridMPI* ptr_dump = dump;
    	Real dumptime = time;
      dumper = new std::thread( [ptr_dump, dumptime, string_file, string_path] () {
        DumpHDF5_MPI_Vector<DumpGridMPI, StreamerHDF5Dump>(*ptr_dump, dumptime, string_file, string_path);
        DumpHDF5_MPI_Channel<DumpGridMPI,StreamerHDF5Dump,3>(*ptr_dump, dumptime, string_file, string_path);
      });
      */
      DumpHDF5_MPI_Vector<FluidGridMPI,StreamerHDF5>(*grid, time, ssR.str(),path4serialization);
    #endif

    if (rank==0) { //saved the grid! Write status to remember most recent ping
      string restart_status = path4serialization+"/restart.status";
      FILE * f = fopen(restart_status.c_str(), "w");
      assert(f != NULL);
      fprintf(f, "time: %20.20e\n", time);
      fprintf(f, "stepid: %d\n", (int)step);
      fprintf(f, "uinfx: %20.20e\n", uinf[0]);
      fprintf(f, "uinfy: %20.20e\n", uinf[1]);
      fprintf(f, "uinfz: %20.20e\n", uinf[2]);
      fclose(f);
      printf("time:  %20.20e\n", time);
      printf("stepid: %d\n", (int)step);
      printf("uinfx: %20.20e\n", uinf[0]);
      printf("uinfy: %20.20e\n", uinf[1]);
      printf("uinfz: %20.20e\n", uinf[2]);
    }

    //CoordinatorDiagnostics coordDiags(grid,time,step);
    //coordDiags(dt);
}

void Simulation::_selectDT()
{
	double local_maxU = (double)findMaxUOMP(vInfo,*grid,uinf);
	double global_maxU;
   const double h = vInfo[0].h_gridpoint;

 		MPI_Allreduce(&local_maxU, &global_maxU, 1, MPI::DOUBLE, MPI::MAX, grid->getCartComm());
    const Real  _CFL = ( step<100 ) ? (step/100.+0.001) * CFL : CFL;
    dtFourier = _CFL*h*h/nu;
    dtCFL     = _CFL*h/abs(global_maxU);
    dt = min(dtCFL,dtFourier);

    if(rank==0) printf("maxU %f dtF %f dtC %f dt %f\n",global_maxU,dtFourier,dtCFL,dt);
    //if (dumpTime>0) dt = min(dt,nextDumpTime-time);
    //if (saveTime>0) dt = min(dt,nextSaveTime-time);
    //if (endTime>0)  dt = min(dt,endTime-time);

    if ( step<100 ) {
        const Real dt_max = 2e-2*CFL;
        const Real dt_min = 2e-6*CFL;
        const Real dt_ramp = dt_min + step*(dt_max - dt_min)/100.;
        if (dt_ramp<dt) {
        	dt = dt_ramp;
        	if(rank==0) printf("Dt bounded by ramp-up: dt_ramp=%f\n",dt_ramp);
        }
    }
}

void Simulation::_serialize(Real & nextSaveTime)
{// FAKE LOL (read the dump) kept it here if I ever decide to implement compressed dumps
	bool bSaving = (saveFreq>0 && step%saveFreq==0)||(saveTime>0 && time>nextSaveTime);
	if(!bSaving) return;
	nextSaveTime += saveTime;
	//DumpZBin_MPI<FluidGridMPI,StreamerSerialization>(*grid,pingname.c_str(),path4serialization);
}

void Simulation::_deserialize()
{
	string restartfile = path4serialization+"/restart.status";
	FILE * f = fopen(restartfile.c_str(), "r");
  if (f == NULL) {
    printf("Could not restart... starting a new sim.\n");
    return;
  }
	assert(f != NULL);

	float val = -1;
	fscanf(f, "time: %e\n", &val);
	assert(val>=0);
	time=val;

	int step_id_fake = -1;
	fscanf(f, "stepid: %d\n", &step_id_fake);
	assert(step_id_fake >= 0);
	step = step_id_fake;
  int ret = 0;
	ret = fscanf(f, "uinfx: %e\n", &val);
  if (ret) uinf[0] = val;
	ret = fscanf(f, "uinfy: %e\n", &val);
  if (ret) uinf[1] = val;
	ret = fscanf(f, "uinfz: %e\n", &val);
  if (ret) uinf[2] = val;
	fclose(f);

	stringstream ssR;
	ssR<<"restart_";
	ssR<<std::setfill('0')<<std::setw(9)<<step;
	if (rank==0) cout << "Restarting from " << ssR.str() << endl;
	MPI_Barrier(grid->getCartComm());
  #if 1
	ReadHDF5_MPI_Vector<FluidGridMPI, StreamerHDF5>(*grid, ssR.str(),path4serialization);
  #else
	ReadHDF5_MPI_Channel<FluidGridMPI, StreamerHDF5, 0>(*grid, ssR.str(),path4serialization);
	ReadHDF5_MPI_Channel<FluidGridMPI, StreamerHDF5, 1>(*grid, ssR.str(),path4serialization);
	ReadHDF5_MPI_Channel<FluidGridMPI, StreamerHDF5, 2>(*grid, ssR.str(),path4serialization);
  #endif
	obstacle_vector->restart(time,path4serialization+"/"+ssR.str());

  printf("DESERIALIZATION: time is %f and step id is %d\n", time, (int)step);
	//bDump=true;
	//_dump("restarted");
	//if (time>0.01) MPI_Abort(grid->getCartComm(), 0);
}

void Simulation::init()
{
    parseArguments();
    setupGrid();
    setupObstacles();
    setupOperators();

    if(bRestart) _deserialize();
    else _ic();

    MPI_Barrier(grid->getCartComm());
}

int Simulation::_3Fish_RLstep(const int nO)
{
  #ifndef __2Leads_
  MPI_Abort(grid->getCartComm(), 0);
  #endif

  const auto _D = obstacle_vector->_getData();
  assert(_D.size() == nO && nO==3);
  bool bDoOver = false;

  bDoOver = _D[2]->checkFail(0.5*(_D[0]->Xrel +_D[1]->Xrel),
                             0.5*(_D[0]->Yrel +_D[1]->Yrel),
                             0.5*(_D[0]->thExp+_D[1]->thExp),
                             length);
  if (bDoOver) {
    _D[2]->finalizePos( 0.5*(_D[0]->Xrel +_D[1]->Xrel),
                        0.5*(_D[0]->Yrel +_D[1]->Yrel),
                        0.5*(_D[0]->thExp+_D[1]->thExp),
                        0.5*(_D[0]->vxExp+_D[1]->vxExp),
                        0.5*(_D[0]->vyExp+_D[1]->vyExp),
                        0.5*(_D[0]->avExp+_D[1]->avExp),
                        length, 1.);
    _D[2]->info = 2;
    obstacle_vector->execute(communicator, 2, time, 0);
    printf("Sim terminated\n");
    return 1;
  }

  if(time >= _D[2]->t_next_comm) {
    _D[2]->finalize(0.5*(_D[0]->Xrel +_D[1]->Xrel),
                    0.5*(_D[0]->Yrel +_D[1]->Yrel),
                    0.5*(_D[0]->thExp+_D[1]->thExp),
                    0.5*(_D[0]->vxExp+_D[1]->vxExp),
                    0.5*(_D[0]->vyExp+_D[1]->vyExp),
                    0.5*(_D[0]->avExp+_D[1]->avExp),
                    length, 1.0);
    _D[2]->finalizePos(0.5*(_D[0]->Xrel +_D[1]->Xrel),
                       0.5*(_D[0]->Yrel +_D[1]->Yrel),
                       0.5*(_D[0]->thExp+_D[1]->thExp),
                       0.5*(_D[0]->vxExp+_D[1]->vxExp),
                       0.5*(_D[0]->vyExp+_D[1]->vyExp),
                       0.5*(_D[0]->avExp+_D[1]->avExp),
                       length, 1.0);
    obstacle_vector->execute(communicator, 2, time, 0);
   }
   return 0;
}
/*
int Simulation::_2Fish_RLstep(const int nO)
{
  const auto _D = obstacle_vector->_getData();
  assert(_D.size() == nO && nO==2);
  bool bDoOver = false;

  for(int i=1; i<nO; i++) {
    bDoOver = _D[i]->checkFail(_D[0]->Xrel, _D[0]->Yrel, _D[0]->thExp, length);
    if (bDoOver) {
      _D[i]->finalizePos(_D[0]->Xrel,  _D[0]->Yrel,  _D[0]->thExp,
                         _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1.);
      _D[i]->info = 2;
      obstacle_vector->execute(communicator, i, time, i-1);
      printf("Sim terminated\n");
      return 1;
    }
  }

  for(int i=0; i<nO; i++) if(time >= _D[i]->t_next_comm) {
    if(i>0) {
       _D[i]->finalize(   _D[0]->Xrel, _D[0]->Yrel, _D[0]->thExp,
                          _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1);
       _D[i]->finalizePos(_D[0]->Xrel, _D[0]->Yrel, _D[0]->thExp,
                          _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1);
     }
    obstacle_vector->execute(communicator, i, time, i-1);
  }
  return 0;
}
*/
int Simulation::_2Fish_RLstep(const int nO)
{
  const auto _D = obstacle_vector->_getData();
  assert(_D.size() == nO && nO<3);
  bool bDoOver = false;

  const Real theta_frame =  _D[0]->thExp + 0.15;
  for(int i=1; i<nO; i++) {
    bDoOver = _D[i]->checkFail(_D[0]->Xrel, _D[0]->Yrel, theta_frame, length);
    if (bDoOver) {
      _D[i]->finalizePos(_D[0]->Xrel,  _D[0]->Yrel, theta_frame,
                 _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1.);
      _D[i]->info = 2;
      obstacle_vector->execute(communicator, i, time, i-1);
      printf("Sim terminated\n");
      return 1;
    }
  }

  for(int i=0; i<nO; i++) if(time >= _D[i]->t_next_comm) {
    if(i>0) {
        _D[i]->finalize(   _D[0]->Xrel, _D[0]->Yrel, theta_frame,
                           _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1);
        _D[i]->finalizePos(_D[0]->Xrel, _D[0]->Yrel, theta_frame,
                           _D[0]->vxExp, _D[0]->vyExp, _D[0]->avExp, length, 1);
    }
    obstacle_vector->execute(communicator, i, time, i-1);
  }
  return 0;
}

void Simulation::simulate()
{
	using namespace std::chrono;
    Real nextDumpTime = dumpTime>0 ? std::ceil(time/dumpTime)*dumpTime : 1e3; //next = if time=0 then 0, if restart time+dumpTime
    if(!rank) printf("First dump at %g (freq %g)\n",nextDumpTime,dumpTime);
    Real nextSaveTime = time + saveTime;
    time_point<high_resolution_clock> last_save, this_save, start_sim;
    start_sim = high_resolution_clock::now();
    last_save = high_resolution_clock::now();

    while (true)
    {
        #ifdef __SMARTIES_
				//if (communicator not_eq nullptr)
				{
					profiler.push_start("RL");
					bool bDoOver = false;
          const int nO = obstacle_vector->nObstacles();
          if (nO == 3) bDoOver = _3Fish_RLstep(nO);
          else         bDoOver = _2Fish_RLstep(nO);
					if (bDoOver) exit(0);
					profiler.pop_stop();
				}
        #endif

        profiler.push_start("DT");
        _selectDT();
        if(bDLM) lambda = 1.0/dt;
        areWeDumping(nextDumpTime);
        profiler.pop_stop();

        for (int c=0; c<pipeline.size(); c++)
        {
        	//cout<<pipeline[c]->getName()<<" start"<<endl;
            profiler.push_start(pipeline[c]->getName());
            (*pipeline[c])(dt);
            profiler.pop_stop();
            //if(time>=0)  _dump(pipeline[c]->getName());
        }
        //if(step>10) abort();
        step++;
        time+=dt;
        if(rank==0)
        printf("Step %d time %f\n",step,time);

        profiler.push_start("Dump");
        if(bDump)
        {
            _dump();

            #ifdef _USE_ZLIB_
              profiler.pop_stop();
              profiler.push_start("Zlib");
              MPI_Barrier(grid->getCartComm());
              const int wtype_write = parser("-wtype_write").asInt(1);
              double threshold = parser("-wthreshold").asDouble(1e-4);

              std::stringstream streamer;
              streamer<<path4serialization;
              streamer<<"/";
              streamer<<"chi_channel_";
              streamer<<std::setfill('0')<<std::setw(9)<<step;
              if(!rank)
                printf("Saving chi to %s, threshold %g\n",
                        streamer.str().c_str(),threshold);

              waveletdumper_grid.verbose();
              waveletdumper_grid.set_threshold(threshold);
              waveletdumper_grid.set_wtype_write(wtype_write);
              waveletdumper_grid.Write<3>(*grid, streamer.str());
            #else
              stringstream ssR;
              ssR<<"restart_"<< std::setfill('0')<<std::setw(9)<<step;
              DumpHDF5_MPI_Channel<FluidGridMPI,StreamerHDF5,3>(*grid, time, ssR.str(), path4serialization);
            #endif

            this_save = high_resolution_clock::now();
            const Real dClock = duration<Real>(this_save-last_save).count();
            const Real totClock = duration<Real>(this_save-start_sim).count();
            if(maxClockDuration < totClock + dClock*1.1) {
            	if(rank==0) {
            	cout<<"Save and exit at time "<<time<<" in "<<step<<" step of "<<nsteps<<endl;
            	cout<<"Not enough allocated clock time ("<<maxClockDuration<<") to reach next dump point"<<endl;
               cout<<"(time since last dump = "<<dClock<<" total duration of the sim = "<<totClock<<")\n"<<endl;
            	}
            	exit(0);
            }
            last_save = high_resolution_clock::now();
        }
        profiler.pop_stop();

        profiler.push_start("Save");
        {
        	_serialize(nextSaveTime);
        }
        profiler.pop_stop();

        if (step % 50 == 0 && !rank) profiler.printSummary();
        if ((endTime>0 && time>endTime) || (nsteps!=0 && step>=nsteps))
        {
        	if(rank==0)
        	cout<<"Finished at time "<<time<<" in "<<step<<" step of "<<nsteps<<endl;
			//exit(0);
		return;
        }
    }
}
