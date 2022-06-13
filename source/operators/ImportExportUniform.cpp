//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch)
//

#include "ImportExportUniform.h"
#include <Cubism/ImportExport.hh>

CubismUP_3D_NAMESPACE_BEGIN

void exportGridFieldToUniformMatrix(FluidGridMPI *grid, int component, Real *out)
{
  exportGridToUniformMatrix<LabMPI>(
      grid,
      [component](const LabMPI::ElementType &element) {
        return element.member(component);
      },            // Getter.
      {component},  // Components vector.
      out);
}

void importGridFieldFromUniformMatrix(
    FluidGridMPI *grid, int component, const Real * __restrict__ in)
{
  importGridFromUniformMatrix(
      grid,
      [component](LabMPI::ElementType &element, Real in) {
        element.member(component) = in;
      },  // Setter.
      in);
}

CubismUP_3D_NAMESPACE_END
