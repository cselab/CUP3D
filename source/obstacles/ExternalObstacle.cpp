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
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

#include <Cubism/ArgumentParser.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace ExternalObstacleObstacle
{

inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0); 
}

//struct FillBlocks : FillBlocksBase<FillBlocks>
struct FillBlocks
{
  const Real safety;
  const Real position[3], quaternion[4];
  Mesh mesh;
  std::vector<Vector3<Real>> randomNormals;
  const Real w=quaternion[0], x=quaternion[1], y=quaternion[2], z=quaternion[3];
  const Real Rmatrix[3][3] = {
      {1-2*(y*y+z*z),   2*(x*y+z*w),   2*(x*z-y*w)},
      {  2*(x*y-z*w), 1-2*(x*x+z*z),   2*(y*z+x*w)},
      {  2*(x*z+y*w),   2*(y*z-x*w), 1-2*(x*x+y*y)}
  };
  std::array<std::array<Real,2>,3> box;
  std::vector<std::vector<int>> BlocksToTriangles;
  Real length;

  FillBlocks(Real _h, const Real p[3], const Real q[4], Mesh &_mesh, std::vector<Vector3<Real>> _randomNomals, std::vector<std::vector<int>> _BlocksToTriangles, Real _length): 
  safety((2+SURFDH)*_h), position{p[0],p[1],p[2]}, quaternion{q[0],q[1],q[2],q[3]}, mesh(_mesh), randomNormals(_randomNomals), BlocksToTriangles(_BlocksToTriangles), length(_length)
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

  inline bool isTouching(const BlockInfo& info, const FluidBlock&b) const
  {
    Real MINP[3], MAXP[3];
    info.pos(MINP, 0, 0, 0);
    info.pos(MAXP, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    const Real intersect[3][2] = {
      {std::max(MINP[0], box[0][0]), std::min(MAXP[0], box[0][1]) },
      {std::max(MINP[1], box[1][0]), std::min(MAXP[1], box[1][1]) },
      {std::max(MINP[2], box[2][0]), std::min(MAXP[2], box[2][1]) }
    };
    return intersect[0][1]-intersect[0][0]>0 &&
           intersect[1][1]-intersect[1][0]>0 &&
           intersect[2][1]-intersect[2][0]>0;
  }

  inline Real signedDistance(const Real px, const Real py, const Real pz, const int id) const
  {
    const std::vector<int> & myTriangles = BlocksToTriangles[id];
    if (myTriangles.size() == 0) return -1; //very far

    // Find the closest triangles and the distance to them.
    Vector3<Real> p = { px, py, pz };
    Real minSqrDist = 1e10;
    std::vector<Vector3<Vector3<Real>>> closest{};
    for (size_t index = 0; index < myTriangles.size(); ++index)
    {
      const int i = myTriangles[index];
      Vector3<Vector3<Real>> t{mesh.x_[mesh.tri_[i][0]],mesh.x_[mesh.tri_[i][1]],mesh.x_[mesh.tri_[i][2]],};
      const Real sqrDist = pointTriangleSqrDistance(t[0],t[1],t[2], p);
      if (std::fabs(sqrDist- minSqrDist)< length * 0.001)
      {
        if (sqrDist < minSqrDist) minSqrDist = sqrDist;
        closest.push_back(t);
      }
      else if (sqrDist < minSqrDist)
      {
        minSqrDist = sqrDist;
        closest.clear();
        closest.push_back(t);
      }
    }
    const Real dist = minSqrDist;

    bool trust = true;    
    Real side = -1;
    for (size_t c = 0; c<closest.size(); c++)
    {
      Vector3<Real> n{};
      n = cross(closest[c][1] - closest[c][0], closest[c][2] - closest[c][0]);
      const Real delta0 = n[0]*closest[c][0][0]+n[1]*closest[c][0][1]+n[2]*closest[c][0][2];
      const Real delta1 = n[0]*closest[c][1][0]+n[1]*closest[c][1][1]+n[2]*closest[c][1][2];
      const Real delta2 = n[0]*closest[c][2][0]+n[1]*closest[c][2][1]+n[2]*closest[c][2][2];
      const Real delta_max = std::max({delta0,delta1,delta2});
      const Real delta_min = std::min({delta0,delta1,delta2});
      const Real delta = std::fabs(delta_max) > std::fabs(delta_min) ? delta_max : delta_min;   
      const Real dot_prod = n[0]*p[0]+n[1]*p[1]+n[2]*p[2];  
      const Real newside = -(dot_prod-delta);
      if (c > 0 && newside*side < 0) trust = false;
      side = newside;
      if (!trust) break;
    }
    if (trust) return std::copysign(dist, side);
    else return isInner(p) ? dist : -dist;
  }

  inline bool isInner(Vector3<Real> p) const
  {
    std::vector<size_t> numIntersections(randomNormals.size(),0   );
    std::vector<bool  > validRay        (randomNormals.size(),true);
    for( size_t i = 0; i<randomNormals.size(); i++ )
    {
        numIntersections[i] = 0;

        std::vector<Vector3<int>> my_triangles;
        Vector3<Real> ray_start = p;
        Vector3<Real> ray_end;
        const Real LL = 3.0;
        ray_end = p + LL*randomNormals[i];
        for( const auto& tri: mesh.tri_ )
        {
          Vector3<Vector3<Real>> t{mesh.x_[tri[0]],mesh.x_[tri[1]],mesh.x_[tri[2]]};
          Vector3<Real> C = (1.0/3.0)*(t[0] + t[1] + t[2]);
          Vector3<Real> proj = ProjectToLine(ray_start, ray_end, C);
          const Real d2= (C[0]-proj[0])*(C[0]-proj[0])
                          +(C[1]-proj[1])*(C[1]-proj[1])
                          +(C[2]-proj[2])*(C[2]-proj[2]);
          if (d2 < (0.5*length)*(0.5*length)) my_triangles.push_back(tri);
        }

        for( const auto& tri: my_triangles )
        {
          Vector3<Vector3<Real>> t{mesh.x_[tri[0]],mesh.x_[tri[1]],mesh.x_[tri[2]]};

          // send ray
          // returns 0 for miss, 1 for hit, -1 for parallel triangle, and -2 for line intersection
          Vector3<Real> intersectionPoint{};
          const int intersection = rayIntersectsTriangle( p, randomNormals[i], t, intersectionPoint );

          // check whether ray is invalid (corner or edge intersection)
          if( intersection >= 0 )
          {
            for( size_t j = 0; j<3; j++ )
            {
              if( std::fabs(t[j][0] - intersectionPoint[0]) < 1e-10 && 
                  std::fabs(t[j][1] - intersectionPoint[1]) < 1e-10 &&
                  std::fabs(t[j][2] - intersectionPoint[2]) < 1e-10 ) validRay[i] = false;
              Vector3<Real> vecA= t[(j+1)%3] - intersectionPoint;
              Vector3<Real> vecB= intersectionPoint - t[j];
              Vector3<Real> vecC= t[(j+1)%3] - t[j];
              Real normA = norm( vecA );
              Real normB = norm( vecB );
              Real normC = norm( vecC );
              if( std::fabs( normA+normB - normC) < 1e-10) validRay[i] = false;
            }
            numIntersections[i] += intersection;
          }
        }
        if( validRay[i] && numIntersections[i]%2 == 0) return false; //definetely outside
        if( validRay[i] && numIntersections[i]%2 == 1) return true;
    }
    std::cout << "Point " << p[0] << " " << p[1] << " " << p[2] << " has no valid rays!" << std::endl;
    for( size_t i = 0; i<randomNormals.size(); i++ )
    {
      std::cout << p[0] + randomNormals[i][0] << " " << p[1] +randomNormals[i][1] << " " << p[2] + randomNormals[i][2] << std::endl; 
    }
    abort();
  }

  using CHIMAT = Real[FluidBlock::sizeZ][FluidBlock::sizeY][FluidBlock::sizeX];
  void operator()(const cubism::BlockInfo &info, ObstacleBlock* const o) const
  {
    FluidBlock &b = *(FluidBlock *)info.ptrBlock;
    if (!isTouching(info, b)) return;
    auto & SDFLAB = o->sdfLab;
    for (int iz = -1; iz < FluidBlock::sizeZ+1; ++iz)
    for (int iy = -1; iy < FluidBlock::sizeY+1; ++iy)
    for (int ix = -1; ix < FluidBlock::sizeX+1; ++ix)
    {
      Real p[3];
      info.pos(p, ix, iy, iz);
      const Real dist = signedDistance(p[0], p[1], p[2], info.blockID);
      SDFLAB[iz+1][iy+1][ix+1] = dist;
    }
  }
};
}

ExternalObstacle::ExternalObstacle(SimulationData& s, ArgumentParser& p): Obstacle(s, p)
{
  // reading coordinates / indices from file
  path = p("-externalObstaclePath").asString();
  if( ExternalObstacleObstacle::exists(path) )
  {
    if( sim.rank == 0 )
        std::cout << "[ExternalObstacle] Reading mesh from " << path << std::endl;

    // Construct the data object by reading from file
    happly::PLYData plyIn(path);

    // Get mesh-style data from the object
    std::vector<std::array<Real, 3>> vPos = plyIn.getVertexPositions();
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
      coordinates.push_back(pt);
    }

    for(const auto& indx: fInd)
    {
      Vector3<int> id = { indx[0], indx[1], indx[2] };
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
  for( size_t i = 0; i< 10; i++ ) {
    Real normRandomNormal = 0.0;
    Vector3<Real> randomNormal;
    while ( std::fabs(normRandomNormal) < 1e-7 ) {
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
  BlocksToTriangles.clear();

  const std::vector<cubism::BlockInfo> vInfo = sim.vInfo();

  BlocksToTriangles.resize(vInfo.size()); //each block has a set of indices(triangles) that are inside it

  //int total = 0;
  const int BS = std::max({FluidBlock::sizeX,FluidBlock::sizeY,FluidBlock::sizeZ});
  for (size_t j = 0 ; j < indices.size() ; j++) //loop over all triangles
  {
    int found = 0;
    for (size_t b = 0 ; b < vInfo.size(); b++)
    {
      const cubism::BlockInfo & info = vInfo[b];

      Real center[3];
      info.pos(center,FluidBlock::sizeX/2,FluidBlock::sizeY/2,FluidBlock::sizeZ/2);

      Vector3<Real> t0;
      t0[0] = position[0] + coordinates[indices[j][0]][0];
      t0[1] = position[1] + coordinates[indices[j][0]][1];
      t0[2] = position[2] + coordinates[indices[j][0]][2];

      Vector3<Real> t1;
      t1[0] = position[0] + coordinates[indices[j][1]][0];
      t1[1] = position[1] + coordinates[indices[j][1]][1];
      t1[2] = position[2] + coordinates[indices[j][1]][2];

      Vector3<Real> t2;
      t2[0] = position[0] + coordinates[indices[j][2]][0];
      t2[1] = position[1] + coordinates[indices[j][2]][1];
      t2[2] = position[2] + coordinates[indices[j][2]][2];

      Vector3<Real> centerV;
      centerV[0] = center[0];
      centerV[1] = center[1];
      centerV[2] = center[2];

      const Real sqrDist = pointTriangleSqrDistance(t0,t1,t2, centerV);
      if (sqrDist < BS * info.h) // info.h * BS/2 * sqrt(3), sqrt(3)/2=0.86, we use 1.0 to be on the safe side 
      {
        //total ++;
        found ++;
        BlocksToTriangles[b].push_back(j);
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&found,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (found == 0)
    {
      std::cout << "There is a triangle that has no nearest blocks, this is impossible." << std::endl;
      abort();
    }
  }
  //int blocks = vInfo.size();
  //MPI_Allreduce(MPI_IN_PLACE,&blocks,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE,&total ,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  //if (sim.rank == 0)
  //  std::cout << " Average triangles per block = " << ((Real) total)/blocks << std::endl;
  //after this loop, each block will have a set of triangles associated with it

  const Real h = sim.hmin;
  Mesh mesh(coordinates, indices);

  ExternalObstacleObstacle::FillBlocks K(h, position, quaternion, mesh, randomNormals, BlocksToTriangles, length);

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
