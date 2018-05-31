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
    //IF3D_ObstacleVector(FluidGridMPI* g) : IF3D_ObstacleOperator(g) {}
    IF3D_ObstacleVector(FluidGridMPI* g, std::vector<IF3D_ObstacleOperator*> obstacles_in)
    : IF3D_ObstacleOperator(g), obstacles(obstacles_in)
    {
      // sort obstacles to make sure that all those that are blocking
      // when chi is compute (ie they do towers from a sdf on the grid
      // rather than computing the sdf on the fly) are last
      // a is put before b if b is blocking and a is not
      const auto isAbeforeB = [&] (
        const IF3D_ObstacleOperator* a, const IF3D_ObstacleOperator* b)
        { return b->isMPIBarrierOnChiCompute and not
                 a->isMPIBarrierOnChiCompute; };
      std::stable_sort(obstacles.begin(), obstacles.end(), isAbeforeB);
    }
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

    void interpolateOnSkin(const double time, const int step, bool dumpWake=false) override;

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
