//
//  IF3D_ObstacleOperator.h
//  IF3D_ROCKS
//
//  Created by Wim van Rees on 04/10/14.
//
//


#ifndef IF3D_ROCKS_IF3D_ObstacleOperator_h
#define IF3D_ROCKS_IF3D_ObstacleOperator_h

#include "common.h"
#include "Definitions.h"
//#include "IF3D_ObstacleLibrary.h"
#include "IF3D_FactoryFileLineParser.h"
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
    FluidGrid * grid;
    surfacePoints surfData;
    vector<BlockInfo> vInfo;
    std::map<int,ObstacleBlock*> obstacleBlocks;
    
    Real quaternion[4]; //representing orientation
    Real position[3], transVel[3], angVel[3], volume, J[6]; // moment of inertia
    Real mass, force[3], torque[3]; //from diagnostics
    Real totChi, surfForce[3], drag, thrust, Pout, PoutBnd, defPower, defPowerBnd, Pthrust, Pdrag, EffPDef, EffPDefBnd; //from compute forces
    double transVel_correction[3], angVel_correction[3], length;
    double *pX, *pY, *pZ; //dumb dum-dumb dum-duuuumb
    int Npts;
    bool bFixToPlanar;


    virtual void _parseArguments(ArgumentParser & parser);
    virtual void _writeComputedVelToFile(const int step_id, const Real t);
    virtual void _writeDiagForcesToFile(const int step_id, const Real t);
    virtual void _makeDefVelocitiesMomentumFree(const double CoM[3]);
    void _finalizeAngVel(Real& AV[3], const Real& J[6], const Real& gam0, const Real& gam1, const Real& gam2);

public:
    int obstacleID;
    IF3D_ObstacleOperator(FluidGrid * grid, ArgumentParser& parser) :
    	grid(grid), obstacleID(0), transVel{0.,0.,0.}, angVel{0.,0.,0.}, volume(0.0),
		J{0.,0.,0.,0.,0.,0.}, pZ(nullptr), pY(nullptr), pX(nullptr), Npts(0)
    {
        vInfo = grid->getBlocksInfo();
        _parseArguments(parser);
    }

    IF3D_ObstacleOperator(FluidGrid * grid):
    grid(grid), obstacleID(0), transVel{0.,0.,0.}, angVel{0.,0.,0.}, volume(0.0),
	J{0.,0.,0.,0.,0.,0.}, pZ(nullptr), pY(nullptr), pX(nullptr), Npts(0)
	{
    	vInfo = grid->getBlocksInfo();
	}
    
    void Accept(ObstacleVisitor * visitor)
    {
    	visitor->visit(this);
    }

    Real getD() const {return length;}

    virtual void computeDiagnostics(const int stepID, const double time, const double* Uinf, const double lambda) ;
    virtual void computeVelocities(const double* Uinf);
    virtual void computeForces(const int stepID, const double time, const double* Uinf, const double NU, const bool bDump);
    virtual void update(const int step_id, const double t, const double dt);
    virtual void save(const int step_id, const double t, std::string filename = std::string());
    virtual void restart(const double t, std::string filename = std::string());
    
    // some non-pure methods
    virtual void create(const int step_id,const double time, const double dt) { }
    
    //methods that work for all obstacles
    const std::map<int,ObstacleBlock*> getObstacleBlocks() const
    {
        return obstacleBlocks;
    }

    void getObstacleBlocks(std::map<int,ObstacleBlock*>*& obstblock_ptr)
    {
        obstblock_ptr = &obstacleBlocks;
    }
    
    virtual void characteristic_function();

    virtual std::vector<int> intersectingBlockIDs(const int buffer);

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
    
    virtual void getTranslationVelocity(Real UT[3]) const;
    virtual void getAngularVelocity(Real & W[3]) const;
    virtual void getCenterOfMass(Real CM[3]) const;
    virtual void setTranslationVelocity(Real UT[3]);
    virtual void setAngularVelocity(const Real W[3]);
    double getForceX() const;
    double getForceY() const;
};

#endif
