//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Wim van Rees on 17/10/14.
//

#ifndef CubismUP_3D_Frenet_h
#define CubismUP_3D_Frenet_h

#include "../../Definitions.h"

#include <limits>

CubismUP_3D_NAMESPACE_BEGIN

struct Frenet2D
{
  static void solve(const int Nm, const Real * const rS,
    const Real * const curv, const Real * const curv_dt,
    Real * const rX, Real * const rY, Real * const vX, Real * const vY,
    Real * const norX, Real * const norY, Real * const vNorX, Real * const vNorY)
    {
        // initial conditions
        rX[0] = 0.0;
        rY[0] = 0.0;
        norX[0] = 0.0;
        norY[0] = 1.0;
        Real ksiX = 1.0;
        Real ksiY = 0.0;

        // velocity variables
        vX[0] = 0.0;
        vY[0] = 0.0;
        vNorX[0] = 0.0;
        vNorY[0] = 0.0;
        Real vKsiX = 0.0;
        Real vKsiY = 0.0;

        for(int i=1; i<Nm; i++)
        {
            // compute derivatives positions
            const Real dksiX = curv[i-1]*norX[i-1];
            const Real dksiY = curv[i-1]*norY[i-1];
            const Real dnuX = -curv[i-1]*ksiX;
            const Real dnuY = -curv[i-1]*ksiY;

            // compute derivatives velocity
            const Real dvKsiX = curv_dt[i-1]*norX[i-1] + curv[i-1]*vNorX[i-1];
            const Real dvKsiY = curv_dt[i-1]*norY[i-1] + curv[i-1]*vNorY[i-1];
            const Real dvNuX = -curv_dt[i-1]*ksiX - curv[i-1]*vKsiX;
            const Real dvNuY = -curv_dt[i-1]*ksiY - curv[i-1]*vKsiY;

            // compute current ds
            const Real ds = rS[i] - rS[i-1];

            // update
            rX[i] = rX[i-1] + ds*ksiX;
            rY[i] = rY[i-1] + ds*ksiY;
            norX[i] = norX[i-1] + ds*dnuX;
            norY[i] = norY[i-1] + ds*dnuY;
            ksiX += ds * dksiX;
            ksiY += ds * dksiY;

            // update velocities
            vX[i] = vX[i-1] + ds*vKsiX;
            vY[i] = vY[i-1] + ds*vKsiY;
            vNorX[i] = vNorX[i-1] + ds*dvNuX;
            vNorY[i] = vNorY[i-1] + ds*dvNuY;
            vKsiX += ds * dvKsiX;
            vKsiY += ds * dvKsiY;

            // normalize unit vectors
            const Real d1 = ksiX*ksiX + ksiY*ksiY;
            const Real d2 = norX[i]*norX[i] + norY[i]*norY[i];
            if(d1>std::numeric_limits<Real>::epsilon())
            {
                const Real normfac = 1.0/std::sqrt(d1);
                ksiX*=normfac;
                ksiY*=normfac;
            }
            if(d2>std::numeric_limits<Real>::epsilon())
            {
                const Real normfac = 1.0/std::sqrt(d2);
                norX[i]*=normfac;
                norY[i]*=normfac;
            }
        }
    }
};

struct Frenet3D
{
  static void solve(const int Nm, const Real * const rS,
    const Real * const curv, const Real * const curv_dt,
    const Real * const tors, const Real * const tors_dt,
    Real * const rX, Real * const rY, Real * const rZ,
    Real * const vX, Real * const vY, Real * const vZ,
    Real * const norX, Real * const norY, Real * const norZ,
    Real * const vNorX, Real * const vNorY, Real * const vNorZ,
    Real * const binX, Real * const binY, Real * const binZ,
    Real * const vBinX, Real * const vBinY, Real * const vBinZ)
    {
        // initial conditions
        rX[0]     = 0.0;// x(s)
        rY[0]     = 0.0;// y(s)
        rZ[0]     = 0.0;// z(s)
        Real ksiX = 1.0;//tangent  vector X-component
        Real ksiY = 0.0;//tangent  vector Y-component
        Real ksiZ = 0.0;//tangent  vector Z-component
        norX[0]   = 0.0;//normal   vector X-component
        norY[0]   = 1.0;//normal   vector Y-component
        norZ[0]   = 0.0;//normal   vector Z-component
        binX[0]   = 0.0;//binormal vector X-component
        binY[0]   = 0.0;//binormal vector Y-component
        binZ[0]   = 1.0;//binormal vector Z-component

        // velocity variables (time derivatives of variables above)
        vX[0]      = 0.0;
        vY[0]      = 0.0;
        vZ[0]      = 0.0;
        Real vKsiX = 0.0;
        Real vKsiY = 0.0;
        Real vKsiZ = 0.0;
        vNorX[0]   = 0.0;
        vNorY[0]   = 0.0;
        vNorZ[0]   = 0.0;
        vBinX[0]   = 0.0;
        vBinY[0]   = 0.0;
        vBinZ[0]   = 0.0;

        for(int i=1; i<Nm; i++)
        {
            // compute RHS for positions
            const Real dksiX = curv[i-1]*norX[i-1];
            const Real dksiY = curv[i-1]*norY[i-1];
            const Real dksiZ = curv[i-1]*norZ[i-1];

            const Real dnuX = -curv[i-1]*ksiX + tors[i-1]*binX[i-1];
            const Real dnuY = -curv[i-1]*ksiY + tors[i-1]*binY[i-1];
            const Real dnuZ = -curv[i-1]*ksiZ + tors[i-1]*binZ[i-1];

            const Real dbinX = -tors[i-1]*norX[i-1];
            const Real dbinY = -tors[i-1]*norY[i-1];
            const Real dbinZ = -tors[i-1]*norZ[i-1];

            // compute RHS for velocities
            const Real dvKsiX = curv_dt[i-1]*norX[i-1] + curv[i-1]*vNorX[i-1];
            const Real dvKsiY = curv_dt[i-1]*norY[i-1] + curv[i-1]*vNorY[i-1];
            const Real dvKsiZ = curv_dt[i-1]*norZ[i-1] + curv[i-1]*vNorZ[i-1];

            const Real dvNuX = -curv_dt[i-1]*ksiX - curv[i-1]*vKsiX + tors_dt[i-1]*binX[i-1] + tors[i-1]*vBinX[i-1];
            const Real dvNuY = -curv_dt[i-1]*ksiY - curv[i-1]*vKsiY + tors_dt[i-1]*binY[i-1] + tors[i-1]*vBinY[i-1];
            const Real dvNuZ = -curv_dt[i-1]*ksiZ - curv[i-1]*vKsiZ + tors_dt[i-1]*binZ[i-1] + tors[i-1]*vBinZ[i-1];

            const Real dvBinX = -tors_dt[i-1]*norX[i-1] - tors[i-1]*vNorX[i-1];
            const Real dvBinY = -tors_dt[i-1]*norY[i-1] - tors[i-1]*vNorY[i-1];
            const Real dvBinZ = -tors_dt[i-1]*norZ[i-1] - tors[i-1]*vNorZ[i-1];

            // compute current ds
            const Real ds = rS[i] - rS[i-1];

            // update positions
            rX[i] = rX[i-1] + ds*ksiX;
            rY[i] = rY[i-1] + ds*ksiY;
            rZ[i] = rZ[i-1] + ds*ksiZ;
            norX[i] = norX[i-1] + ds*dnuX;
            norY[i] = norY[i-1] + ds*dnuY;
            norZ[i] = norZ[i-1] + ds*dnuZ;
            ksiX += ds * dksiX;
            ksiY += ds * dksiY;
            ksiZ += ds * dksiZ;
            binX[i] = binX[i-1] + ds * dbinX;
            binY[i] = binY[i-1] + ds * dbinY;
            binZ[i] = binZ[i-1] + ds * dbinZ;

            // update velocities
            vX[i] = vX[i-1] + ds*vKsiX;
            vY[i] = vY[i-1] + ds*vKsiY;
            vZ[i] = vZ[i-1] + ds*vKsiZ;
            vNorX[i] = vNorX[i-1] + ds*dvNuX;
            vNorY[i] = vNorY[i-1] + ds*dvNuY;
            vNorZ[i] = vNorZ[i-1] + ds*dvNuZ;
            vKsiX += ds * dvKsiX;
            vKsiY += ds * dvKsiY;
            vKsiZ += ds * dvKsiZ;
            vBinX[i] = vBinX[i-1] + ds * dvBinX;
            vBinY[i] = vBinY[i-1] + ds * dvBinY;
            vBinZ[i] = vBinZ[i-1] + ds * dvBinZ;

            // normalize unit vectors
            const Real d1 = ksiX*ksiX + ksiY*ksiY + ksiZ*ksiZ;
            const Real d2 = norX[i]*norX[i] + norY[i]*norY[i] + norZ[i]*norZ[i];
            const Real d3 = binX[i]*binX[i] + binY[i]*binY[i] + binZ[i]*binZ[i];
            if(d1>std::numeric_limits<Real>::epsilon())
            {
                const Real normfac = 1.0/std::sqrt(d1);
                ksiX*=normfac;
                ksiY*=normfac;
                ksiZ*=normfac;
            }
            if(d2>std::numeric_limits<Real>::epsilon())
            {
                const Real normfac = 1.0/std::sqrt(d2);
                norX[i]*=normfac;
                norY[i]*=normfac;
                norZ[i]*=normfac;
            }
            if(d3>std::numeric_limits<Real>::epsilon())
            {
                const Real normfac = 1.0/std::sqrt(d3);
                binX[i]*=normfac;
                binY[i]*=normfac;
                binZ[i]*=normfac;
            }
            #if 0
            std::cout << "Frenet3D: "<< i
            << " rX:"   << rX[i]   << " rY:"   << rY[i]   << " rZ:"   << rZ[i]
            << " ksiX:" << ksiX    << " ksiY:" << ksiY    << " ksiZ:" << ksiZ
            << " norX:" << norX[i] << " norY:" << norY[i] << " norZ:" << norZ[i]
            << " binX:" << binX[i] << " binY:" << binY[i] << " binZ:" << binZ[i]
            << " curv:" << curv[i] << " tors:" << tors[i]
            << " curv_dt:" << curv_dt[i] << " tors_dt:" << tors_dt[i] << std::endl;
            #endif
        }
    }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Frenet_h
