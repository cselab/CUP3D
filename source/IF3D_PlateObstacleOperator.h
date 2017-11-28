//
//  CubismUP_3D
//
//  Written by Ivica Kicic (ikicic@ethz.ch).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#pragma once

#include "IF3D_ObstacleOperator.h"

/*
 * Plate operator.
 *
 * Defined by center location, normal, vector A in the plane, half-size in
 * the A direction, half-size in the remaining direction B (perpendicular to
 * normal and A).
 *
 * The plate has rounded edges.
 *
 * Factory example:
 *     IF3D_PlateObstacle L=0.1 xpos=0.25 ypos=0.25 zpos=0.25 nx=-1 ny=0 nz=0 ax=0 ay=1 az=0 half-a=0.05 half-b=0.1 half-thickness=0.01
 *
 * Factory arguments:
 *     xpos, ypos, zpos  - Position of the center.
 *     nx, ny, nz        - Normal vector (*)
 *     ax, ay, az        - A vector (**)
 *     half-a            - Half-size in the direction of A.
 *     half-b            - Half-size in the direction of B.
 *     half-thickness    - MUST be at least few cell sizes.
 *
 * (*)  Doesn't have to be normalized.
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
    double half_thickness;  // MUST be few times larger than cell size.

public:
    IF3D_PlateObstacleOperator(FluidGridMPI *g,
                               ArgumentParser &p,
                               const Real * const u);

    void _parseArguments(ArgumentParser &parser) override;

    void create(int step_id, double time, double dt, const Real *Uinf) override;
    void finalize(int step_id, double time, double dt, const Real *Uinf) override;
};

