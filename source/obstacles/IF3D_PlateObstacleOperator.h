//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Ivica Kicic (ikicic@ethz.ch).
//

#pragma once

#include "IF3D_ObstacleOperator.h"

/*
 * Plate operator.
 *
 * Defined by center location, side lengths and inclination angle
 *      OR by center location, side lengths, normal vector and the vector of
 *            direction of one of the edge.
 *
 * Factory example:
 *     IF3D_PlateObstacle L=0.1 xpos=0.25 ypos=0.25 zpos=0.25 a=0.1 b=0.2 thickness=0.02 alpha=30
 *     IF3D_PlateObstacle L=0.1 xpos=0.25 ypos=0.25 zpos=0.25 a=0.1 b=0.2 thickness=0.02 nx=-1 ny=0 nz=0 ax=0 ay=1 az=0
 *
 * Factory arguments:
 *     xpos, ypos, zpos - Position of the center.
 *     a                - Half-size in the direction of A.
 *     b                - Half-size in the direction of B.
 *     thickness        - MUST be at least few cell sizes.
 *     alpha            - Inclination angle. (*)
 *     nx, ny, nz       - Normal vector (automatically normalized).
 *     ax, ay, az       - Direction vector of one of the edges. (**)
 *
 * (*) If specified, vectors (nx, ny, nz) and (ax, ay, az) are automatically computed.
 *     For alpha=0,  (nx, ny, nz) = (1, 0, 0) and (ax, ay, az) = (0, 1, 0).
 *     For alpha=45, (nx, ny, nz) = (0.71, 0.71, 0) and (ax, ay, az) = (-0.71, 0.71, 0).
 *
 * (**) Doesn't have to be normalized or perpendicular to n. It's automatically
 *      adjusted to n. Namely, the exact algorithm of fixing n, A and B is:
 *          N = normalized(N)
 *          B = cross(N, A)
 *          A = cross(B, N)
 */
class IF3D_PlateObstacleOperator : public IF3D_ObstacleOperator
{
    // Vectors n, a and b are unit vectors and mutually orthogonal.
    double nx, ny, nz;      // Normal.
    double ax, ay, az;      // A-side vector.
    double bx, by, bz;      // B-side vector.
    double half_a;          // Half-size in A direction.
    double half_b;          // Half-size in B direction.
    double half_thickness;

public:
    IF3D_PlateObstacleOperator(SimulationData & s, ArgumentParser &p);
    IF3D_PlateObstacleOperator(
            SimulationData & s, ObstacleArguments &args,
            double a, double b, double thickness,
            double alpha);
    IF3D_PlateObstacleOperator(
            SimulationData & s, ObstacleArguments &args,
            double a, double b, double thickness,
            double nx, double ny, double nz, double ax, double ay, double az);

    void _from_alpha(double alpha);
    void _init(void);
    void create(int step_id, double time, double dt, const Real *Uinf) override;
    void finalize(int step_id, double time, double dt, const Real *Uinf) override;
};
