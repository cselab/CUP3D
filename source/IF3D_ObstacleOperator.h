//
//  IF3D_ObstacleOperator.h
//  IF3D_ROCKS
//
//  Created by Wim van Rees on 04/10/14.
//
//


#ifndef IF3D_ROCKS_IF3D_ObstacleOperator_h
#define IF3D_ROCKS_IF3D_ObstacleOperator_h

#include "Definitions.h"
#include "IF3D_ObstacleLibrary.h"
#include "IF2D_FactoryFileLineParser.h"
#include <fstream>

// forward declaration of derived class for visitor

class IF3D_ObstacleOperator;
class IF3D_ObstacleVector;

struct ObstacleVisitor
{
    virtual void visit(IF3D_ObstacleOperator * obstacle) {}
    virtual void visit(IF3D_ObstacleVector * obstacle) {}
};

class IF3D_ObstacleOperator
{
protected:
	FluidGridMPI * grid;
    surfacePoints surfData;
    vector<BlockInfo> vInfo;
    std::map<int,ObstacleBlock*> obstacleBlocks;
    
    Real quaternion[4]; //representing orientation
    Real position[3], transVel[3], angVel[3], volume, J[6]; // moment of inertia
    Real mass, force[3], torque[3]; //from diagnostics
    Real totChi, surfForce[3], drag, thrust, Pout, PoutBnd, defPower, defPowerBnd, Pthrust, Pdrag, EffPDef, EffPDefBnd; //from compute forces
    double transVel_correction[3], angVel_correction[3], length;
    //double ext_X, ext_Y, ext_Z;
    double *pX, *pY, *pZ; //dumb dum-dumb dum-duuuumb
    int Npts, rank;
    bool bFixToPlanar;


    virtual void _parseArguments(ArgumentParser & parser);
    virtual void _writeComputedVelToFile(const int step_id, const Real t, const double * uInf);
    virtual void _writeDiagForcesToFile(const int step_id, const Real t);
    virtual void _makeDefVelocitiesMomentumFree(const double CoM[3]);
    void _finalizeAngVel(Real (&AV)[3], const Real (&J)[6], const Real& gam0, const Real& gam1, const Real& gam2);

public:
    int obstacleID;
    bool bFixFrameOfRef;
    IF3D_ObstacleOperator(FluidGridMPI * grid, ArgumentParser& parser) :
    	grid(grid), obstacleID(0), quaternion{1.,0.,0.,0.}, transVel{0.,0.,0.}, angVel{0.,0.,0.},
		volume(0.0), J{0.,0.,0.,0.,0.,0.}, pZ(nullptr), pY(nullptr), pX(nullptr), Npts(0)
    {
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        vInfo = grid->getBlocksInfo();
        /*
        const double extent = grid->maxextent;
        const unsigned int maxbpd = max(grid->NX*FluidBlock::sizeX,
        							max(grid->NY*FluidBlock::sizeY,
        								grid->NZ*FluidBlock::sizeZ));
        const double scale[3] = {
        		(double)(grid->NX*FluidBlock::sizeX)/(double)maxbpd,
        		(double)(grid->NY*FluidBlock::sizeY)/(double)maxbpd,
        		(double)(grid->NZ*FluidBlock::sizeZ)/(double)maxbpd
        };
        ext_X = scale[0]*extent;
        ext_Y = scale[1]*extent;
        ext_Z = scale[2]*extent;
        */
        _parseArguments(parser);
    }

    IF3D_ObstacleOperator(FluidGridMPI * grid):
    grid(grid), obstacleID(0), quaternion{1.,0.,0.,0.}, transVel{0.,0.,0.}, angVel{0.,0.,0.},
	volume(0.0), J{0.,0.,0.,0.,0.,0.}, pZ(nullptr), pY(nullptr), pX(nullptr), Npts(0)
	{
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    	vInfo = grid->getBlocksInfo();
    	/*
        const double extent = grid->maxextent;
        const unsigned int maxbpd = max(grid->NX*FluidBlock::sizeX,
        							max(grid->NY*FluidBlock::sizeY,
        								grid->NZ*FluidBlock::sizeZ));
        const double scale[3] = {
        		(double)(grid->NX*FluidBlock::sizeX)/(double)maxbpd,
        		(double)(grid->NY*FluidBlock::sizeY)/(double)maxbpd,
        		(double)(grid->NZ*FluidBlock::sizeZ)/(double)maxbpd
        };
        ext_X = scale[0]*extent;
        ext_Y = scale[1]*extent;
        ext_Z = scale[2]*extent;
        */
	}
    
    virtual void Accept(ObstacleVisitor * visitor);

    virtual Real getD() const {return length;}

    virtual void computeDiagnostics(const int stepID, const double time, const double* Uinf, const double lambda) ;
    virtual void computeVelocities(const double* Uinf);
    virtual void computeForces(const int stepID, const double time, const double* Uinf, const double NU, const bool bDump);
    virtual void update(const int step_id, const double t, const double dt, const double* Uinf);
    virtual void save(const int step_id, const double t, std::string filename = std::string());
    virtual void restart(const double t, std::string filename = std::string());

    void execute(Communicator * comm, const int iAgent, const double time);
    
    // some non-pure methods
    virtual void create(const int step_id,const double time, const double dt, const double *Uinf) { }
    
    //methods that work for all obstacles
    virtual std::map<int,ObstacleBlock*> getObstacleBlocks() const
    {
        return obstacleBlocks;
    }

    virtual void getObstacleBlocks(std::map<int,ObstacleBlock*>*& obstblock_ptr)
    {
        obstblock_ptr = &obstacleBlocks;
    }
    
    virtual void characteristic_function();

    virtual std::vector<int> intersectingBlockIDs(const int buffer) const;

    virtual ~IF3D_ObstacleOperator()
    {
        for(auto & entry : obstacleBlocks) {
            if(entry.second != nullptr) {
                delete entry.second;
                entry.second = nullptr;
            }
        }
        obstacleBlocks.clear();
    }
    
    virtual void getTranslationVelocity(Real (&UT)[3]) const;
    virtual void getAngularVelocity(Real (&W)[3]) const;
    virtual void getCenterOfMass(Real (&CM)[3]) const;
    virtual void setTranslationVelocity(Real UT[3]);
    virtual void setAngularVelocity(const Real W[3]);
    virtual double getForceX() const;
    virtual double getForceY() const;

    //THIS IS WHY I CAN NEVER HAVE NICE THINGS!!
	template <typename Kernel>
	void compute(vector<Kernel> kernels)
	{
#if 0
		SynchronizerMPI& Synch = grid->sync(kernels[0]);
		const int nthreads = omp_get_max_threads();
		LabMPI * labs = new LabMPI[nthreads];
		vector<BlockInfo> avail0, avail1;
		for(int i=0; i<nthreads; ++i)
			labs[i].prepare(*grid, Synch);

		static int rounds = -1;
		static int one_less = 1;
		if (rounds == -1) {
			char *s = getenv("MYROUNDS");
			if (s != NULL) rounds = atoi(s);
			else 		   rounds = 0;

			char *s2 = getenv("USEMAXTHREADS");
			if (s2 != NULL) one_less = !atoi(s2);
		}

		MPI::COMM_WORLD.Barrier();

		avail0 = Synch.avail_inner();
		const int Ninner = avail0.size();
		BlockInfo * ary0 = &avail0.front();

		int nthreads_first;
		if (one_less) nthreads_first = nthreads-1;
		else 		  nthreads_first = nthreads;

		if (nthreads_first == 0) nthreads_first = 1;
		int Ninner_first = (nthreads_first)*rounds;
		if (Ninner_first > Ninner) Ninner_first = Ninner;
		int Ninner_rest = Ninner - Ninner_first;

#pragma omp parallel num_threads(nthreads_first)
		{
			const int tid = omp_get_thread_num();
			Kernel& kernel=kernels[tid];
			LabMPI& lab = labs[tid];

#pragma omp for schedule(dynamic,1)
			for(int i=0; i<Ninner_first; i++) {
				lab.load(ary0[i], 0);
				kernel(lab, ary0[i], *(FluidBlock*)ary0[i].ptrBlock);
			}
		}

		avail1 = Synch.avail_halo();
		const int Nhalo = avail1.size();
		BlockInfo * ary1 = &avail1.front();

#pragma omp parallel num_threads(nthreads)
		{
			const int tid = omp_get_thread_num();
			Kernel& kernel=kernels[tid];
			LabMPI& lab = labs[tid];

#pragma omp for schedule(dynamic,1)
			for(int i=-Ninner_rest; i<Nhalo; i++) {
				if (i < 0) {
					int ii = i + Ninner;
					lab.load(ary0[ii], 0);
					kernel(lab, ary0[ii], *(FluidBlock*)ary0[ii].ptrBlock);
				} else {
					lab.load(ary1[i], 0);
					kernel(lab, ary1[i], *(FluidBlock*)ary1[i].ptrBlock);
				}
			}
		}

		if(labs!=NULL) {
			delete[] labs;
			labs=NULL;
		}

		MPI::COMM_WORLD.Barrier();
#else
		SynchronizerMPI& Synch = grid->sync(kernels[0]);
		vector<BlockInfo> avail0, avail1;

		const int nthreads = omp_get_max_threads();
		LabMPI * labs = new LabMPI[nthreads];
		for(int i = 0; i < nthreads; ++i)
			labs[i].prepare(*grid, Synch);

		MPI::COMM_WORLD.Barrier();
		avail0 = Synch.avail_inner();
		const int Ninner = avail0.size();

#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			Kernel& kernel=kernels[tid];
			LabMPI& lab = labs[tid];

#pragma omp for schedule(dynamic,1)
			for(int i=0; i<Ninner; i++) {
				BlockInfo info = avail0[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				lab.load(info, 0);
				kernel(lab, info, b); // why is this using the local blockInfo? or is it global? is dh correct?
			}
		}

		avail1 = Synch.avail_halo();
		const int Nhalo = avail1.size();

#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			Kernel& kernel=kernels[tid];
			LabMPI& lab = labs[tid];

#pragma omp for schedule(dynamic,1)
			for(int i=0; i<Nhalo; i++) {
					BlockInfo info = avail1[i];
					FluidBlock& b = *(FluidBlock*)info.ptrBlock;
					lab.load(info, 0);
					kernel(lab, info, b);
			}
		}

		if(labs!=NULL) {
			delete [] labs;
			labs=NULL;
		}

		MPI::COMM_WORLD.Barrier();
#endif
	}
};

#endif
