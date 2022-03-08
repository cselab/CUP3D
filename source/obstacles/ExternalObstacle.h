//
//  Cubism3D
//  Copyright (c) 2021 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#ifndef CubismUP_3D_ExternalObstacle_h
#define CubismUP_3D_ExternalObstacle_h

#include "Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

template <typename T>
struct Vector3
{
    T &operator[](int k) { return x_[k]; }
    const T &operator[](int k) const { return x_[k]; }

    friend Vector3 operator+(const Vector3 &a, const Vector3 &b)
    {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    }

    friend Vector3 operator-(const Vector3 &a, const Vector3 &b)
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }

    friend Vector3 operator*(const T &a, const Vector3 &b)
    {
        return {a*b[0], a*b[1], a*b[2]};
    }

    friend Vector3 cross(const Vector3 &a, const Vector3 &b)
    {
        return Vector3{
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        };
    }

    friend auto dot(Vector3 a, Vector3 b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    friend auto norm(Vector3 &a)
    {
        return std::sqrt( dot(a,a) );
    }

    T x_[3];
};

// Müller-Trumbore algorithm (from https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm)
inline int rayIntersectsTriangle(const Vector3<Real> &rayOrigin, const Vector3<Real> &rayVector, const Vector3<Vector3<Real>> &triangle, Vector3<Real> &intersectionPoint )
{
    // compute triangle edges
    Vector3<Real> edge1 = triangle[1] - triangle[0];
    Vector3<Real> edge2 = triangle[2] - triangle[0];

    // compute determinant
    Vector3<Real> h = cross( rayVector, edge2 );
    Real a = dot( edge1, h );

    // if cos(theta) is close to zero (theta is the angle between edge1 and h, triangle is parallel
    if (std::abs(a) <= norm(edge1) * norm(h) * 1e-10) return -1;

    // invert determinant
    Real f = 1.0/a;

    // solve for u and return if miss
    Vector3<Real> s = rayOrigin - triangle[0];
    Real u = f * dot( s, h );
    if (u < 0.0 || u > 1.0) return 0;

    // solve for v and return if miss
    Vector3<Real> q = cross( s, edge1 );
    Real v = f * dot( rayVector, q );
    if (v < 0.0 || u + v > 1.0) return 0;

    // compute t and intersection point
    Real t = f * dot( edge2, q );

    // ray is in back
    if ( t < 0 ) return -2;

    intersectionPoint = rayOrigin + t * rayVector;

    return 1;
}

inline Vector3<Real> ProjectToLine(const Vector3<Real> a, const Vector3<Real> b, const Vector3<Real> p)
{
    const Real norm_ab = std::sqrt((b[0]-a[0])*(b[0]-a[0])+(b[1]-a[1])*(b[1]-a[1])+(b[2]-a[2])*(b[2]-a[2]));
    const Real proj_a = std::fabs(dot(p-a,b-a))/norm_ab;
    const Real proj_b = std::fabs(dot(p-b,b-a))/norm_ab;
    if (proj_a <= norm_ab && proj_b <= norm_ab) return a + (proj_a/norm_ab) * (b-a); //point is between a and b
    else if (proj_a < proj_b) return a;
    else return b;
}
inline Real pointTriangleSqrDistance( const Vector3<Real> a, const Vector3<Real> b, const Vector3<Real> c, const Vector3<Real> r)
{
    // 1. project in 2D plane containing the triangle
    const auto ac = c - a;
    const auto ab = b - a;
    const auto ar = r - a;
    const auto n = cross(ab, ac);
    const Real alpha = dot(ar, n) / dot(n, n);
    const auto rp = r - alpha * n;

    // 2. compute the barycentric coordinates
    Real u = 0; 
    Real v = 0; 
    Real w = 0; 
    const Vector3<Real> temp = rp - a;
    const Real d00 = dot(ab  , ab);
    const Real d01 = dot(ab  , ac);
    const Real d11 = dot(ac  , ac);
    const Real d20 = dot(temp, ab);
    const Real d21 = dot(temp, ac);
    const Real denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;

    Vector3<Real> rpt;

    if (std::fabs(denom) < 1e-23) //then triangle is a line 
    {
       const Real ab1 = std::fabs(dot(ab,ab));
       const Real ac1 = std::fabs(dot(ac,ac));
       const Real bc1 = std::fabs(dot(b-c,b-c));
       if      (ab1 <= ac1 && ab1 <= bc1) rpt = ProjectToLine(a,c,rp);
       else if (ac1 <= ab1 && ac1 <= bc1) rpt = ProjectToLine(b,c,rp);
       else if (bc1 <= ab1 && bc1 <= ac1) rpt = ProjectToLine(a,b,rp);
       else MPI_Abort(MPI_COMM_WORLD,666);
    }
    else if (u>=0 && v>=0 && w>=0) //point is inside the triangle
    {
       rpt = u * a + v * b + w * c; //projected point in real space
    }
    else if (u < 0 && v < 0)rpt = c;
    else if (u < 0 && w < 0)rpt = b;
    else if (v < 0 && w < 0)rpt = a;
    else if (u < 0) rpt = ProjectToLine(b,c,rp);
    else if (v < 0) rpt = ProjectToLine(a,c,rp);
    else if (w < 0) rpt = ProjectToLine(a,b,rp);

    Real retval = std::sqrt((r[0]-rpt[0])*(r[0]-rpt[0])+(r[1]-rpt[1])*(r[1]-rpt[1])+(r[2]-rpt[2])*(r[2]-rpt[2]));
    return retval;
}

class Mesh
{
  public:
    Mesh(const std::vector<Vector3<Real>> &x, const std::vector<Vector3<int>>  &tri) : x_{x}, tri_{tri}{}

    void rotate( const Real Rmatrix[3][3], const Real position[3])
    {
        for( auto& pt: x_ )
        {
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
    std::vector<Vector3<Real>> x_;
    std::vector<Vector3<int>> tri_;
};

class ExternalObstacle : public Obstacle
{
  std::string path;
  std::vector<Vector3<Real>> coordinates;
  std::vector<Vector3<int>> indices;

  std::mt19937 gen;
  std::normal_distribution<Real> normalDistribution;
  std::vector<Vector3<Real>> randomNormals;

public:
  ExternalObstacle(SimulationData&s,cubism::ArgumentParser&p);

  void create() override;
  void finalize() override;
  void computeVelocities() override;

  std::vector<std::vector<int>> BlocksToTriangles;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ExternalObstacle_h
