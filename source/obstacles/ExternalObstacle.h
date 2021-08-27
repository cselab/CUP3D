//
//  Cubism3D
//  Copyright (c) 2021 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#ifndef CubismUP_3D_ExternalObstacle_h
#define CubismUP_3D_ExternalObstacle_h

#include "Obstacle.h"
#include "extra/triangleMeshSDF.hpp"

CubismUP_3D_NAMESPACE_BEGIN

class ExternalObstacle : public Obstacle
{
  std::string path;
  std::vector<Vector3<Real>> coordinates;
  std::vector<Vector3<int>> indices;
  Real maxSize;

public:
  ExternalObstacle(SimulationData&s,cubism::ArgumentParser&p);

  void create() override;
  void finalize() override;
  void computeVelocities() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Sphere_h
