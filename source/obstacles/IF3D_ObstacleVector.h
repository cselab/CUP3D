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
    : IF3D_ObstacleOperator(s), obstacles(o) {}
    ~IF3D_ObstacleVector();

    void characteristic_function() override;
    int nObstacles() const {return obstacles.size();}
    void computeVelocities() override;
    void update() override;
    void restart(std::string filename = std::string()) override;
    void save(std::string filename = std::string()) override;
    std::vector<int> intersectingBlockIDs(const int buffer) const override;

    void computeForces() override;

    void create() override;
    void finalize() override;
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

    std::vector<std::vector<ObstacleBlock*>*> getAllObstacleBlocks() const
    {
      const size_t Nobs = obstacles.size();
      std::vector<std::vector<ObstacleBlock*>*> ret(Nobs, nullptr);
      for(size_t i=0; i<Nobs; i++) ret[i]= obstacles[i]->getObstacleBlocksPtr();
      return ret;
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
