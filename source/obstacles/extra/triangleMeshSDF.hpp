#include <cmath>
#include <vector>

#ifndef CubismUP_3D_triangleMeshSDF_h
#define CubismUP_3D_triangleMeshSDF_h

template <typename T>
struct Vector3 {
    T &operator[](int k) { return x_[k]; }
    const T &operator[](int k) const { return x_[k]; }

    friend Vector3 operator+(const Vector3 &a, const Vector3 &b) {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    }

    friend Vector3 operator-(const Vector3 &a, const Vector3 &b) {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }

    friend Vector3 operator*(const T &a, const Vector3 &b) {
        return {a*b[0], a*b[1], a*b[2]};
    }

    friend Vector3 cross(const Vector3 &a, const Vector3 &b) {
        return Vector3{
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        };
    }

    friend auto dot(Vector3 a, Vector3 b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    friend auto norm(Vector3 &a) {
        return std::sqrt( dot(a,a) );
    }

    T x_[3];
};

// From https://stackoverflow.com/a/253874
inline bool approximatelyEqual( Real a, Real b, Real eps ) {
    auto epsilon = eps; // std::numeric_limits<Real>::epsilon();
    return std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * epsilon;
};

// Müller-Trumbore algorithm (from https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm)
inline int rayIntersectsTriangle(const Vector3<Real> &rayOrigin,
                                 const Vector3<Real> &rayVector,
                                 const Vector3<Vector3<Real>> &triangle,
                                 Vector3<Real> &intersectionPoint ) {
    // read triangle vertices
    Vector3<Real> vertex0 = triangle[0];
    Vector3<Real> vertex1 = triangle[1];
    Vector3<Real> vertex2 = triangle[2];

    // compute triangle edges
    Vector3<Real> edge1 = vertex1 - vertex0;
    Vector3<Real> edge2 = vertex2 - vertex0;

    // compute determinant
    Vector3<Real> h = cross( rayVector, edge2 );
    Real a = dot( edge1, h );

    // if determinant is close to zero, triangle is parallel
    const Real EPS = std::numeric_limits<Real>::epsilon(); //1e-8;
    if ( std::abs(a) < norm(edge1) * norm(edge2) * EPS  )
        return -1;

    // invert determinant
    Real f = 1.0/a;

    // solve for u and return if miss
    Vector3<Real> s = rayOrigin - vertex0;
    Real u = f * dot( s, h );
    if (u < 0.0 || u > 1.0)
        return 0;

    // solve for v and return if miss
    Vector3<Real> q = cross( s, edge1 );
    Real v = f * dot( rayVector, q );
    if (v < 0.0 || u + v > 1.0)
        return 0;

    // compute t and intersection point
    Real t = f * dot( edge2, q );

    // ray is in back
    if ( t < 0 )
        return -2;

    intersectionPoint = rayOrigin + t * rayVector;
    return 1;
}

inline Real pointTriangleSqrDistance(
        const Vector3<Vector3<Real>> &tri,
        const Vector3<Real> p) {
    const auto n = cross(tri[1] - tri[0], tri[2] - tri[0]);
    const auto p0 = p - tri[0];
    const auto p1 = p - tri[1];
    const auto p2 = p - tri[2];
    const auto n0 = cross(p1, p2);
    const auto n1 = cross(p2, p0);
    const auto n2 = cross(p0, p1);
    const bool inside0 = dot(n0, n) >= 0;
    const bool inside1 = dot(n1, n) >= 0;
    const bool inside2 = dot(n2, n) >= 0;
    if (inside0 && inside1 && inside2) {
        auto dot_ = dot(p - tri[0], n);
        return dot_ * dot_ / dot(n, n);
    }
    Vector3<Real> pu;
    Vector3<Real> vu;
    if (!inside0){
        vu = tri[2] - tri[1];
        pu = p1;
    } else if (!inside1) {
        vu = tri[0] - tri[2];
        pu = p2;
    } else {
        vu = tri[1] - tri[0];
        pu = p0;
    }
    const auto vuvu = dot(vu, vu);
    const auto vupu = dot(vu, pu);
    if (vupu <= 0) {
        return dot(pu, pu);
    } else if (vupu >= vuvu) {
        const auto pv = pu - vu;
        return dot(pv, pv);
    } else {
        const auto out = dot(pu, pu) - vupu * vupu / vuvu;
        return out >= 0 ? out : 0;
    }
}

class Mesh {
public:
    Mesh(const std::vector<Vector3<Real>> &x,
         const std::vector<Vector3<int>>  &tri) :
        x_{x}, tri_{tri}
    { }

    void rotate( const Real Rmatrix[3][3], const double position[3]) {
        for( auto& pt: x_ ) {
            // printf( "Before transformation [%f, %f, %f]\n", pt[0], pt[1], pt[2] );
            // rotate point
            Vector3<Real> pRot = {
              Rmatrix[0][0]*pt[0] + Rmatrix[1][0]*pt[1] + Rmatrix[2][0]*pt[2],
              Rmatrix[0][1]*pt[0] + Rmatrix[1][1]*pt[1] + Rmatrix[2][1]*pt[2],
              Rmatrix[0][2]*pt[0] + Rmatrix[1][2]*pt[1] + Rmatrix[2][2]*pt[2]
            };
            // printf( "after rotation [%f, %f, %f]\n", pt[0], pt[1], pt[2] );
            // translate point
            pt = { pRot[0]+position[0], pRot[1]+position[1], pRot[2]+position[2] };
            // printf( "after transformation [%f, %f, %f]\n", pt[0], pt[1], pt[2] );
        }
    }

    Real nonConvexSDF(Vector3<Real> p, std::vector<Vector3<Real>> randomNormals) const {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        // Find the closest triangles and the distance to them.
        std::vector<Vector3<Vector3<Real>>> closest{};
        Real minSqrDist = 1e100;
        for (int i = 0; i < (int)tri_.size(); ++i) {
          Vector3<Vector3<Real>> t{
              x_[tri_[i][0]],
              x_[tri_[i][1]],
              x_[tri_[i][2]],
          };
          const Real sqrDist = pointTriangleSqrDistance(t, p);
          if( approximatelyEqual( sqrDist, minSqrDist, 1e-1 ) )
              closest.push_back(t);
          else if (sqrDist < minSqrDist) {
              minSqrDist = sqrDist;
              closest.clear();
              closest.push_back(t);
          }
        }
        const auto dist = std::sqrt(minSqrDist);

        // Check on which side of the closest triangle we are. Compute normal and direction vector
        Vector3<Real> n{};
        Vector3<Real> dir{};
        if( closest.size() == 1 ) {
          n = cross(closest[0][1] - closest[0][0], closest[0][2] - closest[0][0]);
          dir = p - closest[0][0];
          #if 1 // NEEDS TO BE 0 if PSEUDO-NORMAL IS 1 !
          const auto side = dot(n, dir);
          return -std::copysign(dist, side);
          #endif
        }
        #if 0// PSEUDO-NORMAL TO COMPUTE SIGN (http://www2.compute.dtu.dk/pubdb/pubs/1833-full.html) - DOES NOT WORK
        else if( closest.size() == 2 ) {
            // Closest point is on an edge, average the normals
            auto n1 = cross(closest[0][1] - closest[0][0], closest[0][2] - closest[0][0]);
            auto n2 = cross(closest[1][1] - closest[1][0], closest[1][2] - closest[1][0]);
            n = n1 + n2;
            // Get closest point
            for( size_t i = 0; i < 3; i++ )
            {
                dir = p - closest[0][i];
                if( approximatelyEqual( minSqrDist, dot(dir, dir), 1e-1 ) )
                    break;
            }
        }
        else {
            // Closest point is on a vertex, angle-weighted average
            for( const auto& triangle : closest ){
                size_t i = 0;
                for( ; i < 3; i++ )
                {
                    dir = p - triangle[i];
                    if( approximatelyEqual( minSqrDist, dot(dir, dir), 1e-1 ) )
                        break;
                }
                // compute angle at edge of triangle
                auto edgeVec1 = triangle[ (i+1)%3 ] - triangle[i];
                auto edgeVec2 = triangle[ (i+2)%3 ] - triangle[i];
                Real cosalphai = dot(edgeVec1, edgeVec2) / (norm(edgeVec1)*norm(edgeVec2));
                Real alphai  = std::acos(cosalphai);
                // compute normal of triangle
                auto ni = cross(triangle[1] - triangle[0], triangle[2] - triangle[0]);
                // normalize
                ni = ( 1 / norm(ni) ) * ni;
                // compute angle averaged normal
                n = n + alphai * ni;
            }
        }
        const auto side = dot(n, dir);
        return -std::copysign(dist, side);
        #else // RAY TO COMPUTE SIGN
        else
        {
          // shoot ray to check whether gridpoint is in- or outside of surface
          bool bInvalid;
          size_t numIntersections;;
          for( size_t i = 0; i<randomNormals.size(); i++ ) {
              bInvalid = false;
              numIntersections = 0;
              size_t nInvalidTriangles = 0;
              size_t nParallelRays = 0;
              for( const auto& tri: tri_ ) {
                  // get triangle points
                  Vector3<Vector3<Real>> t{ x_[tri[0]],
                                            x_[tri[1]],
                                            x_[tri[2]] };

                  // check area of triangle
                  Vector3<Real> side1 = t[1] - t[0];
                  Vector3<Real> side2 = t[2] - t[0];
                  Vector3<Real> areaVector = cross( side1, side2 );
                  Real area = norm( areaVector );
                  if( area < norm(side1) * norm(side2) * std::numeric_limits<Real>::epsilon() ) {
                        nInvalidTriangles++;
                        if( rank == 0 ) {
                            #pragma omp master
                            fprintf(stderr, "WARNING: Invalid triangle detected. Area = |[%e,%e,%e]x[%e,%e,%e]| = |[%e,%e,%e]| = %e\n Coordinates: t[0][0]=%f, t[0][1]=%f, t[0][2]=%f, t[1][0]=%f, t[1][1]=%f, t[1][2]=%f, t[2][0]=%f, t[2][1]=%f, t[2][2]=%f, tri[0]=%d, tri[1]=%d, tri[2]=%d\n", side1[0], side1[1], side1[2], side2[0], side2[1], side2[2], areaVector[0], areaVector[1], areaVector[2], area, t[0][0], t[0][1], t[0][2], t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2], tri[0], tri[1], tri[2]);
                        }
                        continue;
                  }

                  // send ray
                  Vector3<Real> intersectionPoint{};
                  // returns 0 for miss, 1 for hit, -1 for parallel triangle, and -2 for line intersection
                  int intersection = rayIntersectsTriangle( p, randomNormals[i], t, intersectionPoint );

                  // if ray is parallel to triangle
                  if( intersection == -1 ) {
                    nParallelRays++;
                    // fprintf(stderr, "WARNING: Ray is parallel to triangle: randomNormals[%ld][0]=%f, randomNormals[%ld][1]=%f, randomNormals[%ld][2]=%f, t[0][0]=%f, t[0][1]=%f, t[0][2]=%f, t[1][0]=%f, t[1][1]=%f, t[1][2]=%f, t[2][0]=%f, t[2][1]=%f, t[2][2]=%f, tri[0]=%d, tri[1]=%d, tri[2]=%d\n", i, randomNormals[i][0], i, randomNormals[i][1], i, randomNormals[i][2], t[0][0], t[0][1], t[0][2], t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2], tri[0], tri[1], tri[2]);
                    // bInvalid = true;
                  }
                  // if( bInvalid ) break;

                  // check whether ray is invalid
                  if( intersection >= 0 ) {
                      // .. if ray intersects corner
                      for( size_t j = 0; j<3; j++ )
                      if( approximatelyEqual(t[j][0], intersectionPoint[0], std::numeric_limits<Real>::epsilon() ) &&
                          approximatelyEqual(t[j][1], intersectionPoint[1], std::numeric_limits<Real>::epsilon() ) &&
                          approximatelyEqual(t[j][2], intersectionPoint[2], std::numeric_limits<Real>::epsilon() ) ) 
                      {
                          fprintf(stderr, "WARNING: Ray interesects a corner: intersectionPoint[0]=%f, intersectionPoint[1]=%f, intersectionPoint[2]=%f, t[%ld][0]=%f, t[%ld][1]=%f, t[%ld][2]=%f\n", intersectionPoint[0], intersectionPoint[1], intersectionPoint[2], j, t[j][0], j, t[j][1], j, t[j][2]);
                          bInvalid = true;
                      }
                      if( bInvalid ) break;

                      // .. if ray intersects edge (use triangle inequality)
                      for( size_t j = 0; j<3; j++ ){
                          Vector3<double> vecA= t[(j+1)%3] - intersectionPoint;
                          Vector3<double> vecB= intersectionPoint - t[j];
                          Vector3<double> vecC= t[(j+1)%3] - t[j];
                          Real normA = norm( vecA );
                          Real normB = norm( vecB );
                          Real normC = norm( vecC );
                          if( approximatelyEqual( normA+normB, normC, std::numeric_limits<Real>::epsilon() ) ) {
                              fprintf(stderr, "WARNING: Ray [%f,%f,%f] interesects an edge between t[%ld] = [%f,%f,%f] and t[%ld]= [%f,%f,%f] at [%f,%f,%f] (a=%f, b=%f, c=%f)\n", randomNormals[i][0], randomNormals[i][1], randomNormals[i][2], j, t[j][0], t[j][1], t[j][2], (j+1)%3, t[(j+1)%3][0], t[(j+1)%3][1], t[(j+1)%3][2], intersectionPoint[0], intersectionPoint[1], intersectionPoint[2], normA, normB, normC);
                                bInvalid = true;
                          }
                      }
                      if( bInvalid ) break;

                      // count intersection
                      numIntersections += intersection;
                  }
                }
                if( rank == 0 ) {
                    #pragma omp master
                    fprintf(stderr, "[ExternalObstacle] Ignored %ld/%ld triangles and %ld parallel Rays.\n", nInvalidTriangles, tri_.size(), nParallelRays);
                }
                if( not bInvalid ) break;
            }
            if( bInvalid ) {
              fprintf(stderr, "ERROR: Unable to find a valid ray. Aborting..\n");
              fflush(0); abort();
            }
            return numIntersections%2 == 1 ? dist : -dist;
        }
        #endif
    }

    std::vector<Vector3<Real>> x_;
    std::vector<Vector3<int>> tri_;
};

#endif // CubismUP_3D_triangleMeshSDF_h
