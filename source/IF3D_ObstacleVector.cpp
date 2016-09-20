//
//  IF3D_ObstacleVector.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 8/23/13.
//
//

#include "IF3D_ObstacleVector.h"
#include <string>

void IF3D_ObstacleVector::characteristic_function()
{
    for(const auto & obstacle_ptr : obstacles)
        obstacle_ptr->characteristic_function();
}

void IF3D_ObstacleVector::update(const int step_id, const double time, const double dt, const double* Uinf)
{
    for(const auto & obstacle_ptr : obstacles)
        obstacle_ptr->update(step_id,time,dt,Uinf);
}

void IF3D_ObstacleVector::create(const int step_id,const double time, const double dt, const double *Uinf)
{
	for(const auto & obstacle_ptr : obstacles)
		obstacle_ptr->create(step_id,time,dt,Uinf);
}

std::vector<int> IF3D_ObstacleVector::intersectingBlockIDs(const int buffer) const
{
    std::set<int> IDcollection;
    for(const auto & obstacle_ptr : obstacles) {
        std::vector<int> myIDs = obstacle_ptr->intersectingBlockIDs(buffer);
        IDcollection.insert(myIDs.begin(), myIDs.end()); // it's a set, so only unique values are inserted
    }
    return std::vector<int>(IDcollection.begin(), IDcollection.end());
}

void IF3D_ObstacleVector::computeDiagnostics(const int stepID, const double time, const double* Uinf, const double lambda)
{
	for(const auto & obstacle_ptr : obstacles)
		obstacle_ptr->computeDiagnostics(stepID,time,Uinf,lambda);
}

void IF3D_ObstacleVector::computeVelocities(const double* Uinf)
{
	for(const auto & obstacle_ptr : obstacles)
		obstacle_ptr->computeVelocities(Uinf);
}

void IF3D_ObstacleVector::computeForces(const int stepID, const double time, const double* Uinf, const double NU, const bool bDump)
{
	for(const auto & obstacle_ptr : obstacles)
		obstacle_ptr->computeForces(stepID,time,Uinf,NU,bDump);
}

IF3D_ObstacleVector::~IF3D_ObstacleVector()
{
    for(const auto & obstacle_ptr : obstacles)
        delete obstacle_ptr;
    obstacles.clear();
}

void IF3D_ObstacleVector::save(const int step_id, const double t, std::string filename)
{
    std::string fname = (filename == std::string() ? "IF3D_ObstacleVector_restart" : filename);
    int cntr = 0;
    for(const auto & obstacle_ptr : obstacles) {
        obstacle_ptr->save(step_id, t, fname+"_"+std::to_string(cntr));
        cntr++;
    }
}

void IF3D_ObstacleVector::restart(const double t, std::string filename)
{
    std::string fname = (filename == std::string() ? "IF3D_ObstacleVector_restart" : filename);
    int cntr = 0;
    for(const auto & obstacle_ptr : obstacles) {
        obstacle_ptr->restart(t, fname+"_"+std::to_string(cntr));
        cntr++;
    }
}

void IF3D_ObstacleVector::Accept(ObstacleVisitor * visitor)
{
	for(int i=0;i<obstacles.size();++i)
		obstacles[i]->Accept(visitor);
}

void IF3D_ObstacleVector::execute(Communicator * comm, const int iAgent, const double time)
{
   // obstacles[iAgent]->execute(comm, iAgent, time);
}
