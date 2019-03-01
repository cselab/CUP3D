//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef IncompressibleFluids3D_IF3D_ObstacleVector_h
#define IncompressibleFluids3D_IF3D_ObstacleVector_h

#include "obstacles/IF3D_ObstacleOperator.h"

class IF3D_ObstacleVector : public IF3D_ObstacleOperator
{
 protected:
    std::vector<IF3D_ObstacleOperator*> obstacles;

 public:
    IF3D_ObstacleVector(SimulationData&s) : IF3D_ObstacleOperator(s) {}
    IF3D_ObstacleVector(SimulationData&s, std::vector<IF3D_ObstacleOperator*> o)
    : IF3D_ObstacleOperator(s), obstacles(o)
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
    void computeVelocities(const double dt, const Real lambda) override;
    void update(const int step_id, const double t, const double dt, const Real* Uinf) override;
    void restart(const double t, std::string filename = std::string()) override;
    void save(const int step_id, const double t, std::string filename = std::string()) override;
    std::vector<int> intersectingBlockIDs(const int buffer) const override;

    void computeForces(const int stepID, const double time, const double dt, const double NU, const bool bDump) override;

    void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
    void Accept(ObstacleVisitor * visitor) override;

    std::vector<std::array<int, 2>> collidingObstacles();

    void addObstacle(IF3D_ObstacleOperator * obstacle)
    {
        obstacles.push_back(obstacle);
    }

    const std::vector<IF3D_ObstacleOperator*> & getObstacleVector() const
    {
        return obstacles;
    }

    Real getD() const override;


    #ifdef RL_LAYER
      std::vector<StateReward*> _getData();

      void getFieldOfView(const double lengthscale);

      void execute(const int iAgent, const double time, const std::vector<double> action) override;

      void interpolateOnSkin(const double time, const int step, bool dumpWake=false) override;
    #endif
};

#endif
