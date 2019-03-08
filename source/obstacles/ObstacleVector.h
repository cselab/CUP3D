//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_ObstacleVector_h
#define CubismUP_3D_ObstacleVector_h

#include "obstacles/Obstacle.h"
#include <memory>
#include <utility>

CubismUP_3D_NAMESPACE_BEGIN

class ObstacleVector : public Obstacle
{
 public:
    typedef std::vector<std::shared_ptr<Obstacle>> VectorType;

    ObstacleVector(SimulationData&s) : Obstacle(s) {}
    ObstacleVector(SimulationData&s, VectorType o)
        : Obstacle(s), obstacles(std::move(o)) {}

    int nObstacles() const {return obstacles.size();}
    void computeVelocities() override;
    void update() override;
    void restart(std::string filename = std::string()) override;
    void save(std::string filename = std::string()) override;

    void computeForces() override;

    void create() override;
    void finalize() override;
    void Accept(ObstacleVisitor * visitor) override;

    std::vector<std::array<int, 2>> collidingObstacles();

    void addObstacle(std::shared_ptr<Obstacle> obstacle)
    {
        obstacles.emplace_back(std::move(obstacle));
    }

    const VectorType& getObstacleVector() const
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

 protected:
    VectorType obstacles;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ObstacleVector_h
