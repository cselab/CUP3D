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
    // if( not bBoxIntersect ) return false;

    // check whether block is close to surface or inside
    // size_t numIntersections = 0;
    // bool bInvalid;
    // for( size_t i = 0; i<randomNormals.size(); i++ ) { 
    //   bInvalid = false;
    //   for( const auto& tri: mesh.tri_ ) {
    //     // get triangle points
    //     Vector3<Vector3<Real>> t{ mesh.x_[tri[0]],
    //                               mesh.x_[tri[1]],
    //                               mesh.x_[tri[2]] };

    //     // check whether block intersects triangle bounding box
    //     Real minX = std::min( {t[0][0], t[1][0], t[2][0]} ) - safety;
    //     Real maxX = std::max( {t[0][0], t[1][0], t[2][0]} ) + safety;

    //     Real minY = std::min( {t[0][1], t[1][1], t[2][1]} ) - safety;
    //     Real maxY = std::max( {t[0][1], t[1][1], t[2][1]} ) + safety;

    //     Real minZ = std::min( {t[0][2], t[1][2], t[2][2]} ) - safety;
    //     Real maxZ = std::max( {t[0][2], t[1][2], t[2][2]} ) + safety;

    //     const Real intersectTriangle[3][2] = {
    //       {std::max(b.min_pos[0], minX), std::min(b.max_pos[0], maxX) },
    //       {std::max(b.min_pos[1], minY), std::min(b.max_pos[1], maxY) },
    //       {std::max(b.min_pos[2], minZ), std::min(b.max_pos[2], maxZ) }
    //     };

    //     // check if block intersects triangle bounding box
    //     if( intersectTriangle[0][1]-intersectTriangle[0][0]>0 &&
    //         intersectTriangle[1][1]-intersectTriangle[1][0]>0 &&
    //         intersectTriangle[2][1]-intersectTriangle[2][0]>0  )
    //       return true;

    //     // pick lower left corner of block and send ray
    //     Vector3<Real> startPoint = { b.min_pos[0], b.min_pos[1], b.min_pos[2] };
    //     Vector3<Real> intersectionPoint{};
    //     // returns 0 for miss, 1 for hit, -1 for failure
    //     int intersection = rayIntersectsTriangle( startPoint, randomNormals[i], t, intersectionPoint );

    //     // .. if ray is parallel to triangle
    //     if( intersection == -1 ) {
    //       fprintf(stderr, "WARNING: Ray is parallel to triangle: randomNormals[%ld][0]=%f, randomNormals[%ld][1]=%f, randomNormals[%ld][2]=%f, t[0][0]=%f, t[0][1]=%f, t[0][2]=%f, t[1][0]=%f, t[1][1]=%f, t[1][2]=%f, t[2][0]=%f, t[2][1]=%f, t[2][2]=%f, tri[0]=%d, tri[1]=%d, tri[2]=%d\n", i, randomNormals[i][0], i, randomNormals[i][1], i, randomNormals[i][2], t[0][0], t[0][1], t[0][2], t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2], tri[0], tri[1], tri[2]);
    //       bInvalid = true;
    //     }
    //     if( bInvalid ) break;

    //     // .. if ray intersects triangle
    //     if( intersection == 1 ) {
    //       // .. if ray intersects corner
    //       for( size_t j = 0; j<3; j++ )
    //       if( approximatelyEqual(t[j][0], intersectionPoint[0]) &&
    //           approximatelyEqual(t[j][1], intersectionPoint[1]) &&
    //           approximatelyEqual(t[j][2], intersectionPoint[2]) ) {
    //         fprintf(stderr, "WARNING: Ray interesects a corner: intersectionPoint[0]=%f, intersectionPoint[1]=%f, intersectionPoint[2]=%f, t[%ld][0]=%f, t[%ld][1]=%f, t[%ld][2]=%f\n", intersectionPoint[0], intersectionPoint[1], intersectionPoint[2], j, t[j][0], j, t[j][1], j, t[j][2]);
    //         bInvalid = true;
    //       }
    //       if( bInvalid ) break;

    //       // .. if ray intersects edge (use triangle inequality)
    //       for( size_t j = 0; j<3; j++ ){
    //         Vector3<double> vecA= t[j] - intersectionPoint;
    //         Vector3<double> vecB= intersectionPoint - t[(j+1)%3];
    //         Vector3<double> vecC= t[j] - t[(j+1)%3];
    //         Real normA = norm( vecA );
    //         Real normB = norm( vecB );
    //         Real normC = norm( vecC );
    //         if( approximatelyEqual( normA+normB, normC ) ) {
    //           fprintf(stderr, "WARNING: Ray interesects an edge: a=%f, b=%f, c=%f\n", normA, normB, normC);
    //           bInvalid = true;
    //         }
    //       }
    //       if( bInvalid ) break;

    //       // count intersection
    //       numIntersections++;
    //     }
    //   }
    //   if( not bInvalid ) break;
    // }
    // if( bInvalid ) {
    //   fprintf(stderr, "ERROR: Unable to find a valid ray. Aborting..\n");
    //   fflush(0); abort();
    // }
    // if( numIntersections%2 == 1 )
    //   return true;
    // // return to avoid compilation warning
    // return false;
  }

  inline Real signedDistance(const Real px, const Real py, const Real pz) const
  {
    Vector3<Real> p = { px, py, pz };
    // compute signed distance function
    Real dist = -mesh.nonConvexSDF( p, randomNormals ); // pos inside, neg outside
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
    if( sim.rank == 0 )
      std::cout << "nPoints = " << coordinates.size() << ", nTriangles = " << indices.size() << std::endl;
  }
  else{
    fprintf(stderr, "ERROR: Unable to find %s file\n", path.c_str());
    fflush(0); abort();
  }
  // create 10 random vectors to determine if block is inside obstacle
  gen = std::mt19937();
  normalDistribution = std::normal_distribution<Real>(0.0, 1.0);
  for( size_t i = 0; i<10; i++ ) {
    Real normRandomNormal = 0.0;
    Vector3<Real> randomNormal;
    while ( approximatelyEqual(normRandomNormal, 0.0) ) {
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
