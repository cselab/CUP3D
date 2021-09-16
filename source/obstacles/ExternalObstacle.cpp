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
  Mesh mesh;
  const Real safety;
  const Real maxSize;
  const double position[3], quaternion[4];
  const Real w=quaternion[0], x=quaternion[1], y=quaternion[2], z=quaternion[3];
  const Real Rmatrix[3][3] = {
      {1-2*(y*y+z*z),   2*(x*y+z*w),   2*(x*z-y*w)},
      {  2*(x*y-z*w), 1-2*(x*x+z*z),   2*(y*z+x*w)},
      {  2*(x*z+y*w),   2*(y*z-x*w), 1-2*(x*x+y*y)}
  };
  std::array<std::array<Real,2>,3> box;

  FillBlocks(Mesh &_mesh, Real _h, Real _maxSize, const double p[3], const double q[4]):
  mesh(_mesh), safety((2+SURFDH)*_h), maxSize(_maxSize), position{p[0],p[1],p[2]}, quaternion{q[0],q[1],q[2],q[3]} 
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

    // check whether block is close to surface or inside
    if( not bBoxIntersect ) return false;

    // rays start at edges of bounding box
    Real blockEdges[4][2] = { { b.min_pos[0], b.min_pos[1] },
                              { b.min_pos[0], b.max_pos[1] },
                              { b.max_pos[0], b.min_pos[1] },
                              { b.max_pos[0], b.min_pos[1] } };
    bool bInvalid;

    for( const auto& ray: blockEdges ) {
      bInvalid = false;
      size_t numIntersections = 0;
      for( const auto& tri: mesh.tri_ ) {
        // get triangle points
        Vector3<Vector3<Real>> t{ mesh.x_[tri[0]],
                                  mesh.x_[tri[1]],
                                  mesh.x_[tri[2]] };

        // check whether block intersects triangle bounding box
        Real minX = std::min( {t[0][0], t[1][0], t[2][0]} ) - safety;
        Real maxX = std::max( {t[0][0], t[1][0], t[2][0]} ) + safety;

        Real minY = std::min( {t[0][1], t[1][1], t[2][1]} ) - safety;
        Real maxY = std::max( {t[0][1], t[1][1], t[2][1]} ) + safety;

        Real minZ = std::min( {t[0][2], t[1][2], t[2][2]} ) - safety;
        Real maxZ = std::max( {t[0][2], t[1][2], t[2][2]} ) + safety;

        const Real intersectTriangle[3][2] = {
          {std::max(b.min_pos[0], minX), std::min(b.max_pos[0], maxX) },
          {std::max(b.min_pos[1], minY), std::min(b.max_pos[1], maxY) },
          {std::max(b.min_pos[2], minZ), std::min(b.max_pos[2], maxZ) }
        };

        // check if projected triangle intersects projected point
        if( intersectTriangle[0][1]-intersectTriangle[0][0]>0 &&
            intersectTriangle[1][1]-intersectTriangle[1][0]>0 ) {
          // ..if block intersects triangle bounding box
          if( intersectTriangle[2][1]-intersectTriangle[2][0]>0 ) return true;

          // .. if ray intersects corner
          for( size_t i = 0; i<3; i++ )
          if( approximatelyEqual(t[i][0], ray[0]) &&
              approximatelyEqual(t[i][1], ray[1]) ) bInvalid = true;
          if( bInvalid ) break;

          // .. if ray goes through projected triangle
          // from http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
          const Real x1 = t[0][0], x2 = t[1][0], x3 = t[2][0];
          const Real y1 = t[0][1], y2 = t[1][1], y3 = t[2][1];

          Real denominator = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
          Real a = ((y2 - y3)*(ray[0] - x3) + (x3 - x2)*(ray[1] - y3)) / denominator;
          Real b = ((y3 - y1)*(ray[0] - x3) + (x1 - x3)*(ray[1] - y3)) / denominator;
          Real c = 1 - a - b;

          if( 0 <= a && a <= 1 && 0 < b && b <= 1 && 0 <= c && c <= 1 )
            numIntersections++;
          // .. if ray intersects edge
          if( approximatelyEqual(a, 0.0) ||
              approximatelyEqual(b, 0.0) ||
              approximatelyEqual(c, 0.0) ) bInvalid = true;
          if( bInvalid ) break;

        }
      }
      if( numIntersections%2 == 0 && not bInvalid )
        return true;
    }
    if( bInvalid ){
      fprintf(stderr, "ERROR: Could not find valid ray\n");
      fflush(0); abort();
    }
    return false;
  }

  inline Real signedDistance(const Real px, const Real py, const Real pz) const
  {
    // invert translation of obstacle 
    Vector3<Real> p = {px-position[0], py-position[1], pz-position[2]};
    // invert rotation of obstacle
    Vector3<Real> t = {
        Rmatrix[0][0]*p[0] + Rmatrix[0][1]*p[1] + Rmatrix[0][2]*p[2],
        Rmatrix[1][0]*p[0] + Rmatrix[1][1]*p[1] + Rmatrix[1][2]*p[2],
        Rmatrix[2][0]*p[0] + Rmatrix[2][1]*p[1] + Rmatrix[2][2]*p[2]
    };
    // compute signed distance function
    Real dist = -mesh.nonConvexSDF( t, safety ); // pos inside, neg outside
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

    // Initialize vectors of Vector3 required by triangleMeshSDF
    for(const auto& point: vPos)
    {
      Vector3<Real> pt = { scalingFac*point[0], scalingFac*point[1], scalingFac*point[2] };
      coordinates.push_back(pt);
    }

    for(const auto& indx: fInd)
    {
      Vector3<int> id = { indx[0], indx[1], indx[2] };
      indices.push_back(id);
    }
  }
  else{
    fprintf(stderr, "ERROR: Unable to find %s file\n", path.c_str());
    fflush(0); abort();
  }
}

void ExternalObstacle::create()
{
  const Real h = sim.hmin;
  Mesh mesh(coordinates, indices);
  ExternalObstacleObstacle::FillBlocks K(mesh, h, length, position, quaternion);
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
