//
//  CubismUP_3D
//
//  Written by Ivica Kicic (kicici@ethz.ch).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_PlateObstacleOperator.h"
#include "IF3D_ObstacleLibrary.h"


static constexpr Real EPSILON = std::numeric_limits<Real>::epsilon();
static inline Real sqr(const Real x) { return x * x; }

static inline void _normalize(double * const x,
                              double * const y,
                              double * const z) {
    const double norm = std::sqrt(*x * *x + *y * *y + *z * *z);
    assert(norm > 1e-9);
    const double inv = 1.0 / norm;
    *x = inv * *x;
    *y = inv * *y;
    *z = inv * *z;
}

static inline void _normalized_cross(
        const double ax,
        const double ay,
        const double az,
        const double bx,
        const double by,
        const double bz,
        double * const cx,
        double * const cy,
        double * const cz) {
    const double x = ay * bz - az * by;
    const double y = az * bx - ax * bz;
    const double z = ax * by - ay * bx;
    const double norm = std::sqrt(x * x + y * y + z * z);
    assert(norm > 1e-9);
    const double inv = 1.0 / norm;
    *cx = inv * x;
    *cy = inv * y;
    *cz = inv * z;
}


////////////////////////////////////////////////////////////
// PLATE FILL BLOCKS
////////////////////////////////////////////////////////////

namespace PlateObstacle
{
struct FillBlocks : FillBlocksBase<FillBlocks>
{
    const Real cx, cy, cz;  // Center.
    const Real nx, ny, nz;  // Normal. NORMALIZED.
    const Real ax, ay, az;  // A-side vector. NORMALIZED.
    const Real bx, by, bz;  // A-side vector. NORMALIZED.
    const Real half_a;      // Half-size in A direction.
    const Real half_b;      // Half-size in B direction.
    const Real half_thickness;  // Half-thickess. Edges are rounded.

    Real box[3][2];   // Axis-aligned bounding box.

    FillBlocks(const Real cx, const Real cy, const Real cz,
               const Real nx, const Real ny, const Real nz,
               const Real ax, const Real ay, const Real az,
               const Real bx, const Real by, const Real bz,
               const Real half_a,
               const Real half_b,
               const Real half_thickness);

    // Required by FillBlocksBase.
    bool isTouching(const BlockInfo& info, const int buffer_dx = 0) const;
    Real signedDistance(const Real x, const Real y, const Real z) const;

    // Private.
    bool _isTouching(const Real min_pos[3], const Real max_pos[3]) const;
};

FillBlocks::FillBlocks(
        const Real _cx, const Real _cy, const Real _cz,
        const Real _nx, const Real _ny, const Real _nz,
        const Real _ax, const Real _ay, const Real _az,
        const Real _bx, const Real _by, const Real _bz,
        const Real _half_a,
        const Real _half_b,
        const Real _half_thickness)
    : cx(_cx), cy(_cy), cz(_cz),
      nx(_nx), ny(_ny), nz(_nz),
      ax(_ax), ay(_ay), az(_az),
      bx(_bx), by(_by), bz(_bz),
      half_a(_half_a),
      half_b(_half_b),
      half_thickness(_half_thickness)
{
    // Assert normalized.
    assert(std::fabs(nx * nx + ny * ny + nz * nz - 1) < (Real)1e-9);
    assert(std::fabs(ax * ax + ay * ay + az * az - 1) < (Real)1e-9);
    assert(std::fabs(bx * bx + by * by + bz * bz - 1) < (Real)1e-9);

    // Assert n, a and b are mutually orthogonal.
    assert(std::fabs(nx * ax + ny * ay + nz * az) < (Real)1e-9);
    assert(std::fabs(nx * bx + ny * by + nz * bz) < (Real)1e-9);
    assert(std::fabs(ax * bx + ay * by + az * bz) < (Real)1e-9);

    const Real corners[3][4] = {{
        cx - ax * half_a - bx * half_b,
        cx - ax * half_a + bx * half_b,
        cx + ax * half_a - bx * half_b,
        cx + ax * half_a + bx * half_b,
    }, {
        cy - ay * half_a - by * half_b,
        cy - ay * half_a + by * half_b,
        cy + ay * half_a - by * half_b,
        cy + ay * half_a + by * half_b,
    }, {
        cz - az * half_a - bz * half_b,
        cz - az * half_a + bz * half_b,
        cz + az * half_a - bz * half_b,
        cz + az * half_a + bz * half_b,
    }};
    auto min = [](const Real coord[4]) {
        return std::min(std::min(coord[0], coord[1]),
                        std::min(coord[2], coord[3]));
    };
    auto max = [](const Real coord[4]) {
        return std::max(std::max(coord[0], coord[1]),
                        std::max(coord[2], coord[3]));
    };
    box[0][0] = min(corners[0]) - half_thickness;
    box[0][1] = max(corners[0]) + half_thickness;
    box[1][0] = min(corners[1]) - half_thickness;
    box[1][1] = max(corners[1]) + half_thickness;
    box[2][0] = min(corners[2]) - half_thickness;
    box[2][1] = max(corners[2]) + half_thickness;
}

bool FillBlocks::_isTouching(const Real min_pos[3],
                             const Real max_pos[3]) const
{
    return box[0][0] <= max_pos[0] && box[0][1] >= min_pos[0]
        && box[1][0] <= max_pos[1] && box[1][1] >= min_pos[1]
        && box[2][0] <= max_pos[2] && box[2][1] >= min_pos[2];
}

bool FillBlocks::isTouching(const BlockInfo& info,
                            const int buffer_dx) const
{
    Real min_pos[3];
    Real max_pos[3];

    info.pos(min_pos, 0, 0, 0);
    info.pos(max_pos,
             FluidBlock::sizeX,
             FluidBlock::sizeY,
             FluidBlock::sizeZ);
    for (int i = 0; i < 3; ++i) {
      min_pos[i] -= buffer_dx * info.h_gridpoint;
      max_pos[i] += buffer_dx * info.h_gridpoint;
    }
    return _isTouching(min_pos, max_pos);
}

Real FillBlocks::signedDistance(const Real x,
                                const Real y,
                                const Real z) const {
    // Move plane to the center.
    const Real dx = x - cx;
    const Real dy = y - cy;
    const Real dz = z - cz;
    const Real dotn = dx * nx + dy * ny + dz * nz;

    // Project (x, y, z) to the centered plane.
    const Real px = dx - dotn * nx;
    const Real py = dy - dotn * ny;
    const Real pz = dz - dotn * nz;

    // Project into directions a and b.
    const Real dota = px * ax + py * ay + pz * az;
    const Real dotb = px * bx + py * by + pz * bz;

    // Distance to the rectangle edges in the plane coordinate system.
    const Real a = std::fabs(dota) - half_a;
    const Real b = std::fabs(dotb) - half_b;
    const Real n = std::fabs(dotn) - half_thickness;

    if (a <= 0 && b <= 0 && n <= 0) {
        // Inside, return a positive number.
        return -std::min(n, std::min(a, b));
    } else {
        // Outside, return a negative number.
        const Real a0 = std::max((Real)0, a);
        const Real b0 = std::max((Real)0, b);
        const Real n0 = std::max((Real)0, n);
        return -std::sqrt(a0 * a0 + b0 * b0 + n0 * n0);
    }

    // ROUNDED EDGES.
    // return half_thickness - std::sqrt(dotn * dotn + a0 * a0 + b0 * b0);
}
}  // Namespace PlateObstacle.

////////////////////////////////////////////////////////////
// PLATE OBSTACLE OPERATOR
////////////////////////////////////////////////////////////

IF3D_PlateObstacleOperator::IF3D_PlateObstacleOperator(
        FluidGridMPI *g,
        ArgumentParser &p,
        const Real * const u)
    : IF3D_ObstacleOperator(g, p, u)
{
    p.set_strict_mode();
    half_a = (Real)0.5 * p("-a").asDouble();
    half_b = (Real)0.5 * p("-b").asDouble();
    half_thickness = (Real)0.5 * p("-thickness").asDouble();
    p.unset_strict_mode();

    bool has_alpha = p.check("-alpha");
    if (has_alpha) {
        const double alpha = M_PI / 180. * p("-alpha").asDouble();
        nx = std::cos(alpha);
        ny = std::sin(alpha);
        nz = 0;
        ax = -std::sin(alpha);
        ay = std::cos(alpha);
        az = 0;
    } else {
        p.set_strict_mode();
        nx = p("-nx").asDouble();
        ny = p("-ny").asDouble();
        nz = p("-nz").asDouble();
        ax = p("-ax").asDouble();
        ay = p("-ay").asDouble();
        az = p("-az").asDouble();
        p.unset_strict_mode();
    }

    _normalize(&nx, &ny, &nz);
    _normalized_cross(nx, ny, nz, ax, ay, az, &bx, &by, &bz);
    _normalized_cross(bx, by, bz, nx, ny, nz, &ax, &ay, &az);
}

void IF3D_PlateObstacleOperator::create(const int step_id,
                                        const double time,
                                        const double dt,
                                        const Real * const Uinf)
{
    for(auto & entry : obstacleBlocks) delete entry.second;
    obstacleBlocks.clear();

    PlateObstacle::FillBlocks kernel(position[0], position[1], position[2],
                                     nx, ny, nz,
                                     ax, ay, az,
                                     bx, by, bz,
                                     half_a, half_b, half_thickness);

    for (const BlockInfo &info : vInfo) {
        if (kernel.isTouching(info)) {
            assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
            ObstacleBlock * const block = new ObstacleBlock();
            block->clear();
            obstacleBlocks[info.blockID] = block;
        }
    }
}

void IF3D_PlateObstacleOperator::finalize(const int step_id,
                                          const double time,
                                          const double dt,
                                          const Real * const Uinf)
{
#pragma omp parallel
    {
        PlateObstacle::FillBlocks kernel(position[0], position[1], position[2],
                                         nx, ny, nz,
                                         ax, ay, az,
                                         bx, by, bz,
                                         half_a, half_b, half_thickness);


#pragma omp for schedule(static)
        for (size_t i = 0; i < vInfo.size(); ++i) {
            BlockInfo info = vInfo[i];
            auto pos = obstacleBlocks.find(info.blockID);
            if (pos == obstacleBlocks.end()) continue;
            kernel(info, pos->second);
        }
    }
    for (auto &o : obstacleBlocks) o.second->allocate_surface();
}
