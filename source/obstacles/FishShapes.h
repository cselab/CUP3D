//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef CubismUP_3D_FishShapes_h
#define CubismUP_3D_FishShapes_h

#include "../Definitions.h"

#include <cmath>

CubismUP_3D_NAMESPACE_BEGIN

namespace MidlineShapes
{
  /*
  function inputs: xc, yc are n sized arrays which contain the control points of the cubic b spline
  function outputs onto res: assumed to be either the width or the height
  */
  void integrateBSpline(const Real*const xc, const Real*const yc,
  const int n, const Real length, Real*const rS,Real*const res,const int Nm);

  void naca_width(const Real t_ratio, const Real L, Real*const rS,
    Real*const res, const int Nm);
  void stefan_width(const Real L, Real*const rS, Real*const res, const int Nm);
  void stefan_height(const Real L, Real*const rS, Real*const res, const int Nm);
  void danio_width(const Real L, Real*const rS, Real*const res, const int Nm);
  void danio_height(const Real L, Real*const rS, Real*const res, const int Nm);

  void computeWidthsHeights(const std::string &heightName, const std::string &widthName,
                            Real L, Real *rS, Real *height, Real *width, int nM, int mpirank);
}

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_FishShapes_h
