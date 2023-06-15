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

struct FillBlocks
{
  ExternalObstacle * obstacle;
  std::array<std::array<Real,2>,3> box;

  FillBlocks(ExternalObstacle * _obstacle)
  {
    obstacle = _obstacle;

    // Compute maximal extents
    Real MIN = std::numeric_limits<Real>::min();
    Real MAX = std::numeric_limits<Real>::max();
    Vector3<Real> min = { MAX, MAX, MAX };
    Vector3<Real> max = { MIN, MIN, MIN };
    for( const auto& point: obstacle->x_ )
    for( size_t i = 0; i<3; i++ )
    {
      if( point[i] < min[i] ) min[i] = point[i];
      if( point[i] > max[i] ) max[i] = point[i];
    }
    box = {{
     {min[0], max[0]},
     {min[1], max[1]},
     {min[2], max[2]}
    }};
  }

  inline bool isTouching(const BlockInfo& info, const ScalarBlock&b) const
  {
    Real MINP[3], MAXP[3];
    info.pos(MINP, 0, 0, 0);
    info.pos(MAXP, ScalarBlock::sizeX-1, ScalarBlock::sizeY-1, ScalarBlock::sizeZ-1);
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
    const std::vector<int> & myTriangles = obstacle->BlocksToTriangles[id];

    if (myTriangles.size() == 0) return -1; //very far

    const auto & x_ = obstacle->x_;
    const auto & tri_ = obstacle->tri_;
    const auto & length = obstacle->length;

    // Find the closest triangles and the distance to them.
    Vector3<Real> p = { px, py, pz };
    Vector3<Real> dummy;
    Real minSqrDist = 1e10;
    std::vector<int> closest;
    for (size_t index = 0; index < myTriangles.size(); ++index)
    {
      const int i = myTriangles[index];
      Vector3<Vector3<Real>> t{x_[tri_[i][0]],x_[tri_[i][1]],x_[tri_[i][2]],};
      const Real sqrDist = pointTriangleSqrDistance(t[0],t[1],t[2], p, dummy);
      if (std::fabs(sqrDist- minSqrDist)< length * 0.01 * length * 0.01)
      {
        if (sqrDist < minSqrDist) minSqrDist = sqrDist;
        closest.push_back(i);
      }
      else if (sqrDist < minSqrDist)
      {
        minSqrDist = sqrDist;
        closest.clear();
        closest.push_back(i);
      }
    }

    const Real dist = std::sqrt(minSqrDist);

    bool trust = true;    
    Real side = -1;
    for (size_t c = 0; c<closest.size(); c++)
    {
      const int i = closest[c];
      Vector3<Vector3<Real>> t{x_[tri_[i][0]],x_[tri_[i][1]],x_[tri_[i][2]],};
      Vector3<Real> n{};
      n = cross(t[1] - t[0], t[2] - t[0]);
      const Real delta0 = n[0]*t[0][0]+n[1]*t[0][1]+n[2]*t[0][2];
      const Real delta1 = n[0]*t[1][0]+n[1]*t[1][1]+n[2]*t[1][2];
      const Real delta2 = n[0]*t[2][0]+n[1]*t[2][1]+n[2]*t[2][2];
      const Real delta_max = std::max({delta0,delta1,delta2});
      const Real delta_min = std::min({delta0,delta1,delta2});
      const Real delta = std::fabs(delta_max) > std::fabs(delta_min) ? delta_max : delta_min;   
      const Real dot_prod = n[0]*p[0]+n[1]*p[1]+n[2]*p[2];  
      const Real newside = -(dot_prod-delta);
      if (c > 0 && newside*side < 0) trust = false;
      side = newside;
      if (!trust) break;
    }

    if (trust)
    {
      return std::copysign(dist, side);
    }
    else
    {
      return isInner(p) ? dist : -dist;
    }
  }

  inline bool isInner(const Vector3<Real> & p) const
  {
    const auto & x_ = obstacle->x_;
    const auto & tri_ = obstacle->tri_;
    const auto & randomNormals = obstacle->randomNormals;

    for( const auto & randomNormal : randomNormals)
    {
      size_t numIntersections = 0;
      bool validRay = true;
      #if 1
      for( const auto& tri: tri_ )
      {
        Vector3<Vector3<Real>> t{x_[tri[0]],x_[tri[1]],x_[tri[2]]};

        //Send ray. Return 0 for miss, 1 for hit, -1 for parallel triangle, -2 for line intersection
        Vector3<Real> intersectionPoint{};
        const int intersection = rayIntersectsTriangle( p, randomNormal, t, intersectionPoint );

        if( intersection > 0 )
        {
          numIntersections += intersection;
        }
        else if (intersection == -3) // check whether ray is invalid (corner or edge intersection)
        {
          validRay = false;
          break;
        }
      }
      #else
      const auto & position = obstacle->position;
      const auto & IJKToTriangles = obstacle->IJKToTriangles;
      const int n = obstacle->nIJK;
      const Real h = obstacle->hIJK;
      const Real extent = 0.5*n*h;
      const Real h2 = 0.5*h;

      Vector3<Real> Ray = {p[0]-position[0],p[1]-position[1],p[2]-position[2]};
      std::vector<bool> block_visited(n*n*n,false);
      int i = -1;
      int j = -1;
      int k = -1;
      bool done = false;
      while(done == false)
      {
        Ray = Ray + (0.10 * h) * randomNormal;
        i =  round( (extent-h2+Ray[0])/h );
        j =  round( (extent-h2+Ray[1])/h );
        k =  round( (extent-h2+Ray[2])/h );
        i = std::min(i,n-1);
        j = std::min(j,n-1);
        k = std::min(k,n-1);
        i = std::max(i,0);
        j = std::max(j,0);
        k = std::max(k,0);
        const int b = i+j*n+k*n*n;
        if (!block_visited[b])
        {
          block_visited[b] = true;
          Vector3<Real> centerV;
          centerV[0] = -extent + h2 + i*h + position[0];
          centerV[1] = -extent + h2 + j*h + position[1];
          centerV[2] = -extent + h2 + k*h + position[2];
          for (auto & tr : IJKToTriangles[b])
          {
            Vector3<Vector3<Real>> t{x_[tri_[tr][0]],x_[tri_[tr][1]],x_[tri_[tr][2]],};
            Vector3<Real> intersectionPoint{};
            const int intersection = rayIntersectsTriangle( p, randomNormal, t, intersectionPoint );
            if( intersection > 0 )
            {
              if((std::fabs(intersectionPoint[0]-centerV[0]) <= h2) &&
                 (std::fabs(intersectionPoint[1]-centerV[1]) <= h2) &&
                 (std::fabs(intersectionPoint[2]-centerV[2]) <= h2) )
              {
                numIntersections += intersection;
              }
            }
            else if (intersection == -3)
            {
              validRay = false;
              done = true;
              break;
            }
          }
        }
        done = (std::fabs(Ray[0]) > extent || std::fabs(Ray[1]) > extent || std::fabs(Ray[2]) > extent);
      }

      #endif

      if (validRay) return (numIntersections%2 == 1);
    }

    std::cout << "Point " << p[0] << " " << p[1] << " " << p[2] << " has no valid rays!" << std::endl;
    for( size_t i = 0; i<randomNormals.size(); i++ )
    {
      std::cout << p[0] + randomNormals[i][0] << " " << p[1] +randomNormals[i][1] << " " << p[2] + randomNormals[i][2] << std::endl; 
    }
    abort();
  }

  using CHIMAT = Real[ScalarBlock::sizeZ][ScalarBlock::sizeY][ScalarBlock::sizeX];
  void operator()(const cubism::BlockInfo &info, ObstacleBlock* const o) const
  {
    for (int iz = -1; iz < ScalarBlock::sizeZ+1; ++iz)
    for (int iy = -1; iy < ScalarBlock::sizeY+1; ++iy)
    for (int ix = -1; ix < ScalarBlock::sizeX+1; ++ix)
    {
      Real p[3];
      info.pos(p, ix, iy, iz);
      const Real dist = signedDistance(p[0], p[1], p[2], info.blockID);
      o->sdfLab[iz+1][iy+1][ix+1] = dist;
    }
  }
};
}

ExternalObstacle::ExternalObstacle(SimulationData& s, ArgumentParser& p): Obstacle(s, p)
{
  path = p("-externalObstaclePath").asString();
  if( ExternalObstacleObstacle::exists(path) )
  {
    if( sim.rank == 0 )
      std::cout << "[ExternalObstacle] Reading mesh from " << path << std::endl;

    //1.Construct the data object by reading from file and get mesh-style data from the object read
    happly::PLYData plyIn(path);
    std::vector<std::array<Real, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<int>> fInd = plyIn.getFaceIndices<int>();

    //2.Compute maximal extent and ExternalObstacle's center of mass
    Real MIN = std::numeric_limits<Real>::min();
    Real MAX = std::numeric_limits<Real>::max();
    Vector3<Real> min = { MAX, MAX, MAX };
    Vector3<Real> max = { MIN, MIN, MIN };
    Vector3<Real> mean = { 0,0,0 };
    for(const auto& point: vPos)
    for( size_t i = 0; i<3; i++ )
    {
      mean[i] += point[i];
      if( point[i] < min[i] ) min[i] = point[i];
      if( point[i] > max[i] ) max[i] = point[i];
    }
    mean[0] /= vPos.size();
    mean[1] /= vPos.size();
    mean[2] /= vPos.size();
    Vector3<Real> diff = max-min;
    const Real maxSize = std::max({diff[0], diff[1], diff[2]});
    const Real scalingFac = length / maxSize;

    //3.Initialize vectors of Vector3 required by triangleMeshSDF
    for(const auto& point: vPos)
    {
      Vector3<Real> pt = { scalingFac*(point[0]-mean[0]), scalingFac*(point[1]-mean[1]), scalingFac*(point[2]-mean[2]) };
      x_.push_back(pt);
    }
    for(const auto& indx: fInd)
    {
      Vector3<int> id = { indx[0], indx[1], indx[2] };
      tri_.push_back(id);
    }

    if( sim.rank == 0 )
    {
      std::cout << "[ExternalObstacle] Largest extent = " << maxSize << ", target length = " << length << ", scaling factor = " << scalingFac << std::endl;
      std::cout << "[ExternalObstacle] Read grid with nPoints = " << vPos.size() << ", nTriangles = " << fInd.size() << std::endl;      
    }
  }
  else
  {
    fprintf(stderr, "[ExternalObstacle] ERROR: Unable to find %s file\n", path.c_str());
    fflush(0); abort();
  }

  // create 10 random vectors to determine if point is inside obstacle
  gen = std::mt19937();
  normalDistribution = std::normal_distribution<Real>(0.0, 1.0);
  for( size_t i = 0; i< 10; i++ )
  {
    Real normRandomNormal = 0.0;
    Vector3<Real> randomNormal;
    while ( std::fabs(normRandomNormal) < 1e-7 )
    {
      randomNormal = { normalDistribution(gen), normalDistribution(gen), normalDistribution(gen) };
      normRandomNormal = std::sqrt(dot( randomNormal, randomNormal )); //norm of the vector
    }
    randomNormal = ( 1/normRandomNormal )*randomNormal;
    randomNormals.push_back(randomNormal);
  }

  rotate();
}

void ExternalObstacle::create()
{
  const std::vector<cubism::BlockInfo> & chiInfo = sim.chiInfo();
  BlocksToTriangles.clear();
  BlocksToTriangles.resize(chiInfo.size()); //each block has a set of indices(triangles) that are inside it
  const int BS = std::max({ScalarBlock::sizeX,ScalarBlock::sizeY,ScalarBlock::sizeZ});

  #pragma omp parallel for
  for (size_t b = 0 ; b < chiInfo.size(); b++)
  {
    Vector3<Real> dummy;
    const cubism::BlockInfo & info = chiInfo[b];
    Real center[3];
    info.pos(center,ScalarBlock::sizeX/2,ScalarBlock::sizeY/2,ScalarBlock::sizeZ/2);
    Vector3<Real> centerV;
    centerV[0] = center[0]-0.5*info.h;
    centerV[1] = center[1]-0.5*info.h;
    centerV[2] = center[2]-0.5*info.h;

    for (size_t tr = 0 ; tr < tri_.size() ; tr++) //loop over all triangles
    {
      const int v0 = tri_[tr][0];
      const int v1 = tri_[tr][1];
      const int v2 = tri_[tr][2];
      const Vector3<Real> & t0 = x_[v0];
      const Vector3<Real> & t1 = x_[v1];
      const Vector3<Real> & t2 = x_[v2];
      const Real sqrDist = pointTriangleSqrDistance(t0,t1,t2, centerV, dummy);
      if (sqrDist < BS * info.h * BS * info.h * 0.75) // = (info.h * BS/2 * sqrt(3))^2
      {
        #pragma omp critical
        {
          BlocksToTriangles[b].push_back(tr);
        }
      }
    }
  }
  /*
  nIJK = 32;
  hIJK = 1.1*length/nIJK; // about 3% of the object's length, arbitrary choice
  const Real extent = nIJK*hIJK;
  #if 0 //serial implementation
  if (IJKToTriangles.size() == 0)
  {
    IJKToTriangles.resize(nIJK*nIJK*nIJK);
    #pragma omp parallel for collapse (3)
    for (int k = 0 ; k < nIJK; k++)
    for (int j = 0 ; j < nIJK; j++)
    for (int i = 0 ; i < nIJK; i++)
    {
      const int idx = i + nIJK*(j + nIJK*k);
      Vector3<Real> centerV;
      Vector3<Real> rpt;
      rpt[0]=0;
      rpt[1]=0;
      rpt[2]=0;
      centerV[0] = -0.5*extent + 0.5*hIJK + i*hIJK + position[0];
      centerV[1] = -0.5*extent + 0.5*hIJK + j*hIJK + position[1];
      centerV[2] = -0.5*extent + 0.5*hIJK + k*hIJK + position[2];
      for (size_t tr = 0 ; tr < tri_.size() ; tr++) //loop over all triangles
      {
        const int v0 = tri_[tr][0];
        const int v1 = tri_[tr][1];
        const int v2 = tri_[tr][2];
        const Vector3<Real> & t0 = x_[v0];
        const Vector3<Real> & t1 = x_[v1];
        const Vector3<Real> & t2 = x_[v2];
        pointTriangleSqrDistance(t0,t1,t2,centerV,rpt);
        if((std::fabs(rpt[0]-centerV[0]) < 0.51*hIJK) &&
           (std::fabs(rpt[1]-centerV[1]) < 0.51*hIJK) &&
           (std::fabs(rpt[2]-centerV[2]) < 0.51*hIJK) )
        {
          #pragma omp critical
          {
            IJKToTriangles[idx].push_back(tr);
          }
        }
      }
    }
  }
  #else //MPI implementation
  if (IJKToTriangles.size() == 0)
  {
    int size;
    MPI_Comm_size(sim.comm,&size);
    int rank = sim.rank;

    const int total_load = nIJK*nIJK*nIJK;
    const int my_load = (rank < total_load % size) ? (total_load / size + 1) : (total_load / size);
    int mystart = (total_load / size) * rank;
    mystart += (rank < (total_load % size)) ? rank : (total_load % size);
    const int myend = mystart + my_load;

    //Create a local vector of vectors
    std::vector<std::vector<int>> local_data(my_load);

    #pragma omp parallel for
    for (int idx = mystart ; idx < myend ; idx++)
    {
      const int k = idx/(nIJK*nIJK);
      const int j = (idx - k*nIJK*nIJK)/nIJK;
      const int i = (idx - k*nIJK*nIJK-j*nIJK)%nIJK;
      Vector3<Real> centerV;
      Vector3<Real> rpt;
      rpt[0]=0;
      rpt[1]=0;
      rpt[2]=0;
      centerV[0] = -0.5*extent + 0.5*hIJK + i*hIJK + position[0];
      centerV[1] = -0.5*extent + 0.5*hIJK + j*hIJK + position[1];
      centerV[2] = -0.5*extent + 0.5*hIJK + k*hIJK + position[2];
      for (size_t tr = 0 ; tr < tri_.size() ; tr++) //loop over all triangles
      {
        const int v0 = tri_[tr][0];
        const int v1 = tri_[tr][1];
        const int v2 = tri_[tr][2];
        const Vector3<Real> & t0 = x_[v0];
        const Vector3<Real> & t1 = x_[v1];
        const Vector3<Real> & t2 = x_[v2];
        pointTriangleSqrDistance(t0,t1,t2,centerV,rpt);
        if((std::fabs(rpt[0]-centerV[0]) < 0.51*hIJK) &&
           (std::fabs(rpt[1]-centerV[1]) < 0.51*hIJK) &&
           (std::fabs(rpt[2]-centerV[2]) < 0.51*hIJK) )
        {
          #pragma omp critical
          {
            local_data[idx-mystart].push_back(tr);
          }
        }
      }
    }

    // Flatten the local vectors for communication
    std::vector<int> send_data;
    for (const auto& vec : local_data)
      send_data.insert(send_data.end(), vec.begin(), vec.end());

    // Communicate the local vectors among ranks
    std::vector<int> recv_counts(size);
    const int send_count = send_data.size();
    MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, sim.comm);

    std::vector<int> displacements(size);
    int total_size = 0;
    for (int i = 0; i < size; i++)
    {
      displacements[i] = total_size;
      total_size += recv_counts[i];
    }

    // Allocate memory for the received global vectors
    std::vector<int> recv_data(total_size);

    // Communicate the local vectors and receive the global vectors
    MPI_Allgatherv(send_data.data(), send_data.size(), MPI_INT, recv_data.data(), recv_counts.data(), displacements.data(), MPI_INT, sim.comm);

    std::vector<int> vector_displacements(total_load,0);
    for (int idx = mystart ; idx < myend ; idx++)
    {
      vector_displacements[idx] = local_data[idx-mystart].size();
    }
    MPI_Allreduce(MPI_IN_PLACE, vector_displacements.data(), vector_displacements.size(), MPI_INT, MPI_SUM, sim.comm);

    //Reconstruct the global vector of vectors
    size_t current_pos = 0;
    for (int idx = 0; idx < total_load; idx++)
    {
      int count = vector_displacements[idx];
      std::vector<int> vec(recv_data.begin() + current_pos, recv_data.begin() + current_pos + count);
      IJKToTriangles.push_back(vec);
      current_pos += count;
    }
  }
  #endif
  */
  ExternalObstacleObstacle::FillBlocks K(this);

  create_base<ExternalObstacleObstacle::FillBlocks>(K);
}

CubismUP_3D_NAMESPACE_END


