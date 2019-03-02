//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

#include "obstacles/IF3D_ExternalObstacleOperator.h"
#include "obstacles/extra/IF3D_ObstacleLibrary.h"

using namespace cubismup3d;

namespace {

// TODO: The position shift should be done here, not in the external code.
struct FillBlocksExternal : FillBlocksBase<FillBlocksExternal>
{
    const ExternalObstacleSettings &S;

    FillBlocksExternal(const ExternalObstacleSettings &_S) : S(_S) { }

    bool isTouching(const BlockInfo &info, const int buffer_dx = 0) const
    {
        // Get block bounding box.
        Real min[3], max[3];
        info.pos(min, -buffer_dx, -buffer_dx, -buffer_dx);
        info.pos(max,
                 FluidBlock::sizeX - 1 + buffer_dx,
                 FluidBlock::sizeY - 1 + buffer_dx,
                 FluidBlock::sizeZ - 1 + buffer_dx);

        // Convert to double.
        const double lo[3] = {(double)min[0], (double)min[1], (double)min[2] };
        const double hi[3] = {(double)max[0], (double)max[1], (double)max[2] };

        // Ask the external code if it the block is overlapping the box.
        return S.is_touching_fn(S.obj, lo, hi);
    }

    Real signedDistance(const Real x, const Real y, const Real z) const
    {
        // Ask the external code what's the signed distance.
        return S.signed_distance_fn(S.obj, (double)x, (double)y, (double)z);
    }

    /*
     * Fill out `ObstacleBlock::udef` by calling the external velocity function.
     */
    void setVelocity(const BlockInfo &info, ObstacleBlock * const o) const {
        for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
        for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
        for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
            Real p[3];
            info.pos(p, ix, iy, iz);

            double udef[3];
            S.velocity_fn(S.obj,
                          (double)p[0],
                          (double)p[1],
                          (double)p[2],
                          udef);

            o->udef[iz][iy][ix][0] = (Real)udef[0];
            o->udef[iz][iy][ix][1] = (Real)udef[1];
            o->udef[iz][iy][ix][2] = (Real)udef[2];
        }
    }
};

}  // namespace


IF3D_ExternalObstacleOperator::IF3D_ExternalObstacleOperator(
        SimulationData&s, const ObstacleArguments &args)
        : IF3D_ObstacleOperator(s, args)
{
    bForcedInSimFrame = {true, true, true};
    bFixFrameOfRef = {true, true, true};
    bBlockRotation = {true, true, true};
}

void IF3D_ExternalObstacleOperator::computeVelocities() {
    IF3D_ObstacleOperator::computeVelocities();

    if (settings.com_velocity_fn != nullptr) {
        settings.com_velocity_fn(settings.obj, transVel_imposed);
        transVel[0] = transVel_imposed[0];
        transVel[1] = transVel_imposed[1];
        transVel[2] = transVel_imposed[2];
    }
}

/*
 * Finds out which blocks are affected (overlap with the object).
 */
void IF3D_ExternalObstacleOperator::create()
{
    if (!settings.is_touching_fn
            || !settings.signed_distance_fn
            || !settings.velocity_fn) {
        // `Simulation` automatically invokes `create` during initalization.
        // Instead of changing the code there, we ignore that `create` request here.
        return;
    }
    assert(settings.obj != nullptr);

    const FillBlocksExternal kernel(settings);
    create_base<FillBlocksExternal>(kernel);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)vInfo.size(); ++i) {
        const BlockInfo &info = vInfo[i];
        if (obstacleBlocks[info.blockID] == nullptr) continue;
        kernel.setVelocity(info, obstacleBlocks[info.blockID]);
    }
}

void IF3D_ExternalObstacleOperator::finalize()
{
  // this method allows any computation that requires the char function
  // to be computed. E.g. compute the effective center of mass or removing
  // momenta from udef
}
