//
//  IF2D_Frenet.h
//  IF2D_ROCKS
//
//  Created by Wim van Rees on 17/10/14.
//
//

#ifndef IF2D_ROCKS_IF2D_Frenet_h
#define IF2D_ROCKS_IF2D_Frenet_h
#include <limits>

struct IF2D_Frenet2D
{
    static void solve(const int Nm, const double * const rS, const double * const curv, const double * const curv_dt, double * const rX, double * const rY, double * const vX, double * const vY, double * const norX, double * const norY, double * const vNorX, double * const vNorY)
    {
        // initial conditions
        rX[0] = 0.0;
        rY[0] = 0.0;
        norX[0] = 0.0;
        norY[0] = 1.0;
        double ksiX = 1.0;
        double ksiY = 0.0;
        
        // velocity variables
        vX[0] = 0.0;
        vY[0] = 0.0;
        vNorX[0] = 0.0;
        vNorY[0] = 0.0;
        double vKsiX = 0.0;
        double vKsiY = 0.0;
        
        for(unsigned int i=1; i<Nm; i++)
        {
            // compute derivatives positions
            const double dksiX = curv[i-1]*norX[i-1];
            const double dksiY = curv[i-1]*norY[i-1];
            const double dnuX = -curv[i-1]*ksiX;
            const double dnuY = -curv[i-1]*ksiY;
            
            // compute derivatives velocity
            const double dvKsiX = curv_dt[i-1]*norX[i-1] + curv[i-1]*vNorX[i-1];
            const double dvKsiY = curv_dt[i-1]*norY[i-1] + curv[i-1]*vNorY[i-1];
            const double dvNuX = -curv_dt[i-1]*ksiX - curv[i-1]*vKsiX;
            const double dvNuY = -curv_dt[i-1]*ksiY - curv[i-1]*vKsiY;
            
            // compute current ds
            const double ds = rS[i] - rS[i-1];

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
            const double d1 = ksiX*ksiX + ksiY*ksiY;
            const double d2 = norX[i]*norX[i] + norY[i]*norY[i];
            if(d1>std::numeric_limits<double>::epsilon())
            {
                const double normfac = 1.0/std::sqrt(d1);
                ksiX*=normfac;
                ksiY*=normfac;
            }
            if(d2>std::numeric_limits<double>::epsilon())
            {
                const double normfac = 1.0/std::sqrt(d2);
                norX[i]*=normfac;
                norY[i]*=normfac;
            }
        }
    }
};

#endif
