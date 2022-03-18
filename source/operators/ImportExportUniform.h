//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch)
//

#pragma once

#include "../Definitions.h"

CubismUP_3D_NAMESPACE_BEGIN

/// Export `component`th component of the fluid elements to the contiguous
/// fully refined matrix.
void exportGridFieldToUniformMatrix(
    FluidGridMPI *grid, int component, Real * __restrict__ out);

/// Import `component`th component of the fluid elements from the contiguous
/// fully refined matrix.
void importGridFieldFromUniformMatrix(
    FluidGridMPI *grid, int component, const Real * __restrict__ out);

CubismUP_3D_NAMESPACE_END
