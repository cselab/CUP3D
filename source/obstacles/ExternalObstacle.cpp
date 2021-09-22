//
//  Cubism3D
//  Copyright (c) 2021 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "ExternalObstacle.h"
#include "extra/ObstacleLibrary.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "extra/happly.h"
#pragma GCC diagnostic pop
#include <filesystem>

#include <Cubism/ArgumentParser.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace ExternalObstacleObstacle
{
struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real safety;
  const double position[3], quaternion[4];
  Mesh mesh;
  std::vector<Vector3<Real>> randomNormals;
  const Real w=quaternion[0], x=quaternion[1], y=quaternion[2], z=quaternion[3];
  const Real Rmatrix[3][3] = {
      {1-2*(y*y+z*z),   2*(x*y+z*w),   2*(x*z-y*w)},
      {  2*(x*y-z*w), 1-2*(x*x+z*z),   2*(y*z+x*w)},
      {  2*(x*z+y*w),   2*(y*z-x*w), 1-2*(x*x+y*y)}
  };
  std::array<std::array<Real,2>,3> box;

  FillBlocks(Real _h, const double p[3], const double q[4], Mesh &_mesh, std::vector<Vector3<Real>> _randomNomals): 
  safety((2+SURFDH)*_h), position{p[0],p[1],p[2]}, quaternion{q[0],q[1],q[2],q[3]}, mesh(_mesh), randomNormals(_randomNomals)
  {
    mesh.rotate( Rmatrix, position );
    // Compute maximal extents
    Real MIN = std::numeric_limits<Real>::min();
    Real MAX = std::numeric_limits<Real>::max();
    Vector3<Real> min = { MAX, MAX, MAX };
    Vector3<Real> max = { MIN, MIN, MIN };
    for( const auto& point: mesh.x_ )
    for( size_t i = 0; i<3; i++ )
    {
      if( point[i] < min[i] ) min[i] = point[i];
      if( point[i] > max[i] ) max[i] = point[i];
    }
    Vector3<Real> diff = max-min;
    box = {{
     {(Real)position[0]-(diff[0]+safety), (Real)position[0]+(diff[0]+safety)},
     {(Real)position[1]-(diff[1]+safety), (Real)position[1]+(diff[1]+safety)},
     {(Real)position[2]-(diff[2]+safety), (Real)position[2]+(diff[2]+safety)}
    }};
  }

  inline bool isTouching(const BlockInfo&, const FluidBlock&b) const
  {
    // check intersection of block with obstacle bounding box
    const Real intersect[3][2] = {
      {std::max(b.min_pos[0], box[0][0]), std::min(b.max_pos[0], box[0][1]) },
      {std::max(b.min_pos[1], box[1][0]), std::min(b.max_pos[1], box[1][1]) },
      {std::max(b.min_pos[2], box[2][0]), std::min(b.max_pos[2], box[2][1]) }
    };
    bool bBoxIntersect = intersect[0][1]-intersect[0][0]>0 && intersect[1][1]-intersect[1][0]>0 && intersect[2][1]-intersect[2][0]>0;

    return bBoxIntersect;
  }

  inline Real signedDistance(const Real px, const Real py, const Real pz) const
  {
    Vector3<Real> p = { px, py, pz };
    // compute signed distance function
    Real dist = mesh.nonConvexSDF( p, randomNormals ); // pos inside, neg outside
    return dist;
  }
};
}

ExternalObstacle::ExternalObstacle(SimulationData& s, ArgumentParser& p)
    : Obstacle(s, p)
{
  path = p("-externalObstaclePath").asString();
  // reading coordinates / indices from file
  if( std::filesystem::exists(path) )
  {
    if( sim.rank == 0 )
        std::cout << "[ExternalObstacle] Reading mesh from " << path << std::endl;
    // Construct the data object by reading from file
    happly::PLYData plyIn(path);

    // Get mesh-style data from the object
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<int>> fInd = plyIn.getFaceIndices<int>();

    // Compute maximal extent
    Real MIN = std::numeric_limits<Real>::min();
    Real MAX = std::numeric_limits<Real>::max();
    Vector3<Real> min = { MAX, MAX, MAX };
    Vector3<Real> max = { MIN, MIN, MIN };
    for(const auto& point: vPos)
    for( size_t i = 0; i<3; i++ )
    {
      if( point[i] < min[i] ) min[i] = point[i];
      if( point[i] > max[i] ) max[i] = point[i];
    }
    Vector3<Real> diff = max-min;
    Real maxSize = std::max({diff[0], diff[1], diff[2]});
    Real scalingFac = length / maxSize;
    if( sim.rank == 0 )
      std::cout << "[ExternalObstacle] Largest extent = " << maxSize << ", target length = " << length << ", scaling factor = " << scalingFac << std::endl;

    // Initialize vectors of Vector3 required by triangleMeshSDF
    for(const auto& point: vPos)
    {
      Vector3<Real> pt = { scalingFac*point[0], scalingFac*point[1], scalingFac*point[2] };
      // if( sim.rank == 0 ) {
      //   std::cout << "reading point [" << point[0] << ", " << point[1] << ", " << point[2] << "]" << std::endl;
      //   std::cout << "after scaling [" << pt[0] << ", " << pt[1] << ", " << pt[2] << "]" << std::endl;
      // }
      coordinates.push_back(pt);
    }

    for(const auto& indx: fInd)
    {
      Vector3<int> id = { indx[0], indx[1], indx[2] };
      // if( sim.rank == 0 )
      //   std::cout << "reading point indices [" << id[0] << ", " << id[1] << ", " << id[2] << "]" << std::endl;
      indices.push_back(id);
    }
    if( sim.rank == 0 )
      std::cout << "[ExternalObstacle] Read grid with nPoints = " << coordinates.size() << ", nTriangles = " << indices.size() << std::endl;
  }
  else{
    fprintf(stderr, "[ExternalObstacle] ERROR: Unable to find %s file\n", path.c_str());
    fflush(0); abort();
  }
  // create 10 random vectors to determine if point is inside obstacle
  gen = std::mt19937();
  normalDistribution = std::normal_distribution<Real>(0.0, 1.0);
  for( size_t i = 0; i<10; i++ ) {
    Real normRandomNormal = 0.0;
    Vector3<Real> randomNormal;
    while ( approximatelyEqual(normRandomNormal, 0.0, 1e-7) ) {
      randomNormal = { normalDistribution(gen), 
                       normalDistribution(gen), 
                       normalDistribution(gen) };
      normRandomNormal = norm( randomNormal );
    }
    randomNormal = ( 1/normRandomNormal )*randomNormal;
    randomNormals.push_back(randomNormal);
  }
}

void ExternalObstacle::create()
{
  const Real h = sim.hmin;
  Mesh mesh(coordinates, indices);
  ExternalObstacleObstacle::FillBlocks K(h, position, quaternion, mesh, randomNormals);
  create_base<ExternalObstacleObstacle::FillBlocks>(K);
}

void ExternalObstacle::finalize()
{
  // this method allows any computation that requires the char function
  // to be computed. E.g. compute the effective center of mass or removing
  // momenta from udef
}


void ExternalObstacle::computeVelocities()
{
  Obstacle::computeVelocities();
}

CubismUP_3D_NAMESPACE_END
