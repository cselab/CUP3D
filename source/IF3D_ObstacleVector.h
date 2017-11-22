//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef IncompressibleFluids3D_IF3D_ObstacleVector_h
#define IncompressibleFluids3D_IF3D_ObstacleVector_h

#include "IF3D_ObstacleOperator.h"

class IF3D_ObstacleVector : public IF3D_ObstacleOperator
{
 protected:
    std::vector<IF3D_ObstacleOperator*> obstacles;

 public:
    IF3D_ObstacleVector(FluidGridMPI* g) : IF3D_ObstacleOperator(g) {}
    IF3D_ObstacleVector(FluidGridMPI* g, std::vector<IF3D_ObstacleOperator*> obstacles_in)
    : IF3D_ObstacleOperator(g), obstacles(obstacles_in)
    {}
    ~IF3D_ObstacleVector();

    void characteristic_function() override;
    int nObstacles() const {return obstacles.size();}
    void computeVelocities(const Real* Uinf) override;
    void update(const int step_id, const double t, const double dt, const Real* Uinf) override;
    void restart(const double t, std::string filename = std::string()) override;
    void save(const int step_id, const double t, std::string filename = std::string()) override;
    std::vector<int> intersectingBlockIDs(const int buffer) const override;

    void computeDiagnostics(const int stepID, const double time, const Real* Uinf, const double lambda) override;
    void computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump) override;

    void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
    void Accept(ObstacleVisitor * visitor) override;

    void getFieldOfView(const double lengthscale);
    vector<std::array<int, 2>> collidingObstacles();
    std::vector<StateReward*> _getData();

    void execute(const int iAgent, const double time, const vector<double> action) override;

    void interpolateOnSkin(const double time, const int step, const int iobst = -1, bool dumpWake=false);

    void addObstacle(IF3D_ObstacleOperator * obstacle)
    {
        obstacles.push_back(obstacle);
    }

    const std::vector<IF3D_ObstacleOperator*> & getObstacleVector() const
    {
        return obstacles;
    }

    Real getD() const override;
};

#endif
