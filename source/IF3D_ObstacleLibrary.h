//
//  IF3D_ObstacleLibrary.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 8/6/13.
//
//

#ifndef IncompressibleFluids3D_IF3D_ObstacleLibrary_h
#define IncompressibleFluids3D_IF3D_ObstacleLibrary_h
#include <array>
#include "IF2D_Interpolation1D.h"
#define _4thORDER_Towers
namespace SphereObstacle
{
    struct FillBlocks
    {
        const Real radius,safe_radius;
        const Real sphere_position[3];
        Real sphere_box[3][2];
        
        void _find_sphere_box()
        {
            sphere_box[0][0] = sphere_position[0] - safe_radius;
            sphere_box[0][1] = sphere_position[0] + safe_radius;
            sphere_box[1][0] = sphere_position[1] - safe_radius;
            sphere_box[1][1] = sphere_position[1] + safe_radius;
            sphere_box[2][0] = sphere_position[2] - safe_radius;
            sphere_box[2][1] = sphere_position[2] + safe_radius;
        }
        
        FillBlocks(const Real radius, const Real max_dx, const Real pos[3]):
        radius(radius),safe_radius(radius+6*max_dx), sphere_position{pos[0],pos[1],pos[2]} //4h should be safe
        {            
            _find_sphere_box();
        }
        
        bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
        {
            Real intersection[3][2] = {
                max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
                max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
                max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
            };
            
            return
            intersection[0][1]-intersection[0][0]>0 &&
            intersection[1][1]-intersection[1][0]>0 &&
            intersection[2][1]-intersection[2][0]>0;
        }
        
        bool _is_touching(const BlockInfo& info, const int buffer_dx = 0) const
        {
            Real min_pos[3], max_pos[3];
            
            info.pos(min_pos, 0,0,0);
            info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            for(int i=0;i<3;++i) {
                min_pos[i]-=buffer_dx*info.h_gridpoint;
                max_pos[i]+=buffer_dx*info.h_gridpoint;
            }
            return _is_touching(min_pos,max_pos);
        }
        
        inline Real distanceToSphere(const Real x, const Real y, const Real z) const
        {
            return radius - std::sqrt(x*x+y*y+z*z); // pos inside, neg outside
        }
        
        inline void operator()(const BlockInfo& info,FluidBlock& b,ObstacleBlock* const defblock,surfaceBlocks* const surf) const
        {
            if(_is_touching(info)) {
            	const Real h = info.h_gridpoint;
#ifdef _4thORDER_Towers
        		const Real inv2h = 1./(12.*h);
#else
        		const Real inv2h = 0.5/h;
#endif
        		const Real eps = std::numeric_limits<Real>::epsilon();

                for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                    for(int iy=0; iy<FluidBlock::sizeY; iy++)
                        for(int ix=0; ix<FluidBlock::sizeX; ix++)
                        {
                            Real p[3];
                            info.pos(p, ix, iy, iz);
                            const Real x = p[0]-sphere_position[0];
                            const Real y = p[1]-sphere_position[1];
                            const Real z = p[2]-sphere_position[2];
                            const Real dist = distanceToSphere(x,y,z);
                            /*
                            const bool in_band= lab(ix,iy,iz)*lab(ix+1,iy,iz) <=0 ||
												lab(ix,iy,iz)*lab(ix-1,iy,iz) <=0 ||
												lab(ix,iy,iz)*lab(ix,iy+1,iz) <=0 ||
												lab(ix,iy,iz)*lab(ix,iy-1,iz) <=0 ||
												lab(ix,iy,iz)*lab(ix,iy,iz+1) <=0 ||
												lab(ix,iy,iz)*lab(ix,iy,iz-1) <=0;
							if(not in_band)
                             */
                            if(dist > +3*h || dist < -3*h) { //2 should be safe
                            	const Real H = dist > 0 ? 1.0 : 0.0;
                            	defblock->chi[iz][iy][ix] = H;
                            	b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
                            	continue;
                            }

                            const Real distPx = distanceToSphere(x+h,y,z);
                            const Real distMx = distanceToSphere(x-h,y,z);
                            const Real distPy = distanceToSphere(x,y+h,z);
                            const Real distMy = distanceToSphere(x,y-h,z);
                            const Real distPz = distanceToSphere(x,y,z+h);
                            const Real distMz = distanceToSphere(x,y,z-h);

                            const Real IplusX = distPx < 0 ? 0 : distPx;
                            const Real IminuX = distMx < 0 ? 0 : distMx;
                            const Real IplusY = distPy < 0 ? 0 : distPy;
                            const Real IminuY = distMy < 0 ? 0 : distMy;
                            const Real IplusZ = distPz < 0 ? 0 : distPz;
                            const Real IminuZ = distMz < 0 ? 0 : distMz;
                            const Real HplusX = distPx == 0 ? 0.5 : (distPx < 0 ? 0 : 1);
                            const Real HminuX = distMx == 0 ? 0.5 : (distMx < 0 ? 0 : 1);
                            const Real HplusY = distPy == 0 ? 0.5 : (distPy < 0 ? 0 : 1);
                            const Real HminuY = distMy == 0 ? 0.5 : (distMy < 0 ? 0 : 1);
                            const Real HplusZ = distPz == 0 ? 0.5 : (distPz < 0 ? 0 : 1);
                            const Real HminuZ = distMz == 0 ? 0.5 : (distMz < 0 ? 0 : 1);
#ifdef _4thORDER_Towers //the advantage of having 4th order derivatives is to spread the transition from chi 0 to 1 over 4 points
                            const Real distP2x= distanceToSphere(x+2*h,y,z);
                            const Real distM2x= distanceToSphere(x-2*h,y,z);
                            const Real distP2y= distanceToSphere(x,y+2*h,z);
                            const Real distM2y= distanceToSphere(x,y-2*h,z);
                            const Real distP2z= distanceToSphere(x,y,z+2*h);
                            const Real distM2z= distanceToSphere(x,y,z-2*h);
                            const Real Iplus2X= distP2x < 0 ? 0 : distP2x;
                            const Real Iminu2X= distM2x < 0 ? 0 : distM2x;
                            const Real Iplus2Y= distP2y < 0 ? 0 : distP2y;
                            const Real Iminu2Y= distM2y < 0 ? 0 : distM2y;
                            const Real Iplus2Z= distP2z < 0 ? 0 : distP2z;
                            const Real Iminu2Z= distM2z < 0 ? 0 : distM2z;
                            const Real Hplus2X= distP2x == 0 ? 0.5 : (distP2x < 0 ? 0 : 1);
                            const Real Hminu2X= distM2x == 0 ? 0.5 : (distM2x < 0 ? 0 : 1);
                            const Real Hplus2Y= distP2y == 0 ? 0.5 : (distP2y < 0 ? 0 : 1);
                            const Real Hminu2Y= distM2y == 0 ? 0.5 : (distM2y < 0 ? 0 : 1);
                            const Real Hplus2Z= distP2z == 0 ? 0.5 : (distP2z < 0 ? 0 : 1);
                            const Real Hminu2Z= distM2z == 0 ? 0.5 : (distM2z < 0 ? 0 : 1);
                            const Real gradUX = inv2h * (-distP2x + 8*distPx - 8*distMx + distM2x);
                            const Real gradUY = inv2h * (-distP2y + 8*distPy - 8*distMy + distM2y);
                            const Real gradUZ = inv2h * (-distP2z + 8*distPz - 8*distMz + distM2z);
                            const Real gradIX = inv2h * (-Iplus2X + 8*IplusX - 8*IminuX + Iminu2X);
                            const Real gradIY = inv2h * (-Iplus2Y + 8*IplusY - 8*IminuY + Iminu2Y);
                            const Real gradIZ = inv2h * (-Iplus2Z + 8*IplusZ - 8*IminuZ + Iminu2Z);
                            const Real gradHX = inv2h * (-Hplus2X + 8*HplusX - 8*HminuX + Hminu2X);
                            const Real gradHY = inv2h * (-Hplus2Y + 8*HplusY - 8*HminuY + Hminu2Y);
                            const Real gradHZ = inv2h * (-Hplus2Z + 8*HplusZ - 8*HminuZ + Hminu2Z);
#else if defined(_UPWIND_Towers)
                            Real gradUX, gradUY, gradUZ, gradIX, gradIY, gradIZ, gradHX, gradHY, gradHZ;
                            if (std::fabs(distPx) < std::fabs(distMx)) {
                            } else {
                            } 
                            const Real distP2x= distanceToSphere(x+2*h,y,z);
                            const Real distM2x= distanceToSphere(x-2*h,y,z);
                            const Real distP2y= distanceToSphere(x,y+2*h,z);
                            const Real distM2y= distanceToSphere(x,y-2*h,z);
                            const Real distP2z= distanceToSphere(x,y,z+2*h);
                            const Real distM2z= distanceToSphere(x,y,z-2*h);
                            const Real Iplus2X= distP2x < 0 ? 0 : distP2x;
                            const Real Iminu2X= distM2x < 0 ? 0 : distM2x;
                            const Real Iplus2Y= distP2y < 0 ? 0 : distP2y;
                            const Real Iminu2Y= distM2y < 0 ? 0 : distM2y;
                            const Real Iplus2Z= distP2z < 0 ? 0 : distP2z;
                            const Real Iminu2Z= distM2z < 0 ? 0 : distM2z;
                            const Real Hplus2X= distP2x == 0 ? 0.5 : (distP2x < 0 ? 0 : 1);
                            const Real Hminu2X= distM2x == 0 ? 0.5 : (distM2x < 0 ? 0 : 1);
                            const Real Hplus2Y= distP2y == 0 ? 0.5 : (distP2y < 0 ? 0 : 1);
                            const Real Hminu2Y= distM2y == 0 ? 0.5 : (distM2y < 0 ? 0 : 1);
                            const Real Hplus2Z= distP2z == 0 ? 0.5 : (distP2z < 0 ? 0 : 1);
                            const Real Hminu2Z= distM2z == 0 ? 0.5 : (distM2z < 0 ? 0 : 1);
                            const Real gradUX = (distPx<distMx) ? inv2h * (-distP2x + 8*distPx - 8*distMx + distM2x);
                            const Real gradUY = inv2h * (-distP2y + 8*distPy - 8*distMy + distM2y);
                            const Real gradUZ = inv2h * (-distP2z + 8*distPz - 8*distMz + distM2z);
                            const Real gradIX = inv2h * (-Iplus2X + 8*IplusX - 8*IminuX + Iminu2X);
                            const Real gradIY = inv2h * (-Iplus2Y + 8*IplusY - 8*IminuY + Iminu2Y);
                            const Real gradIZ = inv2h * (-Iplus2Z + 8*IplusZ - 8*IminuZ + Iminu2Z);
                            const Real gradHX = inv2h * (-Hplus2X + 8*HplusX - 8*HminuX + Hminu2X);
                            const Real gradHY = inv2h * (-Hplus2Y + 8*HplusY - 8*HminuY + Hminu2Y);
                            const Real gradHZ = inv2h * (-Hplus2Z + 8*HplusZ - 8*HminuZ + Hminu2Z);
#else
                            const Real gradUX = inv2h * (distPx - distMx);
                            const Real gradUY = inv2h * (distPy - distMy);
                            const Real gradUZ = inv2h * (distPz - distMz);
							       const Real gradIX = inv2h * (IplusX - IminuX);
                            const Real gradIY = inv2h * (IplusY - IminuY);
                            const Real gradIZ = inv2h * (IplusZ - IminuZ);
                            const Real gradHX = inv2h * (HplusX - HminuX);
                            const Real gradHY = inv2h * (HplusY - HminuY);
                            const Real gradHZ = inv2h * (HplusZ - HminuZ);
#endif
                            const Real gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ;
                            const Real numH    = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;
                            const Real numD    = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;
                            const Real Delta = std::abs(gradUSq) < eps ? numD : numD / gradUSq;
                            const Real H     = std::abs(gradUSq) < eps ? numH : numH / gradUSq;

                            if (Delta>1e-6) { //will always be multiplied by h^2
                            	const Real dchidx = -Delta*gradUX;
                            	const Real dchidy = -Delta*gradUY;
                            	const Real dchidz = -Delta*gradUZ;
                            	surf->add(info.blockID,ix,iy,iz,dchidx,dchidy,dchidz,Delta);
                            }

                            assert(H>=0-eps && H<=1+eps);
#ifndef NDEBUG
                            if(H<0 || H>1)
							printf("invalid H?: %10.10e %10.10e %10.10e: %10.10e\n",x,y,z,H);
#endif
                			defblock->chi[iz][iy][ix] = H;
                			b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
                        b(ix,iy,iz).tmpU = dist;
                        b(ix,iy,iz).tmpV = Delta;
                        }
            }
        }
    };
}

namespace TorusObstacle
{
    struct FillBlocks
    {
        const Real big_r, small_r, smoothing_length;
        Real position[3];
        Real sphere_box[3][2];
        
        void _find_sphere_box()
        {
            sphere_box[0][0] = position[0] - 2*(small_r + smoothing_length);
            sphere_box[0][1] = position[0] + 2*(small_r + smoothing_length);
            sphere_box[1][0] = position[1] - 2*(big_r + small_r + smoothing_length);
            sphere_box[1][1] = position[1] + 2*(big_r + small_r + smoothing_length);
            sphere_box[2][0] = position[2] - 2*(big_r + small_r + smoothing_length);
            sphere_box[2][1] = position[2] + 2*(big_r + small_r + smoothing_length);
        }
        
        FillBlocks(Real big_r, Real small_r, Real smoothing_length,  Real position[3]):
        big_r(big_r), small_r(small_r), smoothing_length(smoothing_length)
        {
            this->position[0] = position[0];
            this->position[1] = position[1];
            this->position[2] = position[2];
            
            _find_sphere_box();
        }
        
        FillBlocks(const FillBlocks& c):
        small_r(c.small_r), big_r(c.big_r), smoothing_length(c.smoothing_length)
        {
            position[0] = c.position[0];
            position[1] = c.position[1];
            position[2] = c.position[2];
            
            _find_sphere_box();
        }
        
        bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
        {
            Real intersection[3][2] = {
                max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
                max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
                max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
            };
            
            return
            intersection[0][1]-intersection[0][0]>0 &&
            intersection[1][1]-intersection[1][0]>0 &&
            intersection[2][1]-intersection[2][0]>0;
        }
        
        bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
        {
            Real min_pos[3], max_pos[3];
            
            info.pos(min_pos, 0,0,0);
            info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            for(int i=0;i<3;++i)
            {
                min_pos[i]-=buffer_dx*info.h_gridpoint;
                max_pos[i]+=buffer_dx*info.h_gridpoint;
            }
            return _is_touching(min_pos,max_pos);
        }
        
        Real mollified_heaviside(const Real x, const Real eps) const
        {
            const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
            
            return 0.5+0.5*cos(alpha);
        }

        
        inline void operator()(const BlockInfo& info, FluidBlock& b) const
        {
            if(_is_touching(info))
                //	if(true)
            {
                for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                    for(int iy=0; iy<FluidBlock::sizeY; iy++)
                        for(int ix=0; ix<FluidBlock::sizeX; ix++)
                        {
                            Real p[3];
                            info.pos(p, ix, iy, iz);
                            
                            const Real t[3] = {
                                (Real)(p[0] - position[0]),
                                (Real)(p[1] - position[1]),
                                (Real)(p[2] - position[2])
                            };
                            
                            //mappele //sqrt(t[1] * t[1] + t[2] * t[2]);//
                            const Real r = sqrt(t[1]*t[1] + t[2]*t[2]);
                            
                            if (r > 0)
                            {
                                const Real c[3] = {
                                    0,
                                    big_r*t[1]/r,
                                    big_r*t[2]/r
                                };
                                
                                const Real d = sqrt(pow(t[0]-c[0], 2) +  pow(t[1]-c[1], 2) + pow(t[2]-c[2], 2));
                                const Real chi = mollified_heaviside(d-small_r, smoothing_length);
                                b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
                            }
                        }
            }
        }
    };
}

namespace EllipsoidObstacle
{
    // code from http://www.geometrictools.com/
    
    static Real DistancePointEllipseSpecial (const Real e[2], const Real y[2], Real x[2])
    {
        Real distance = (Real)0;
        if (y[1] > (Real)0)
        {
            if (y[0] > (Real)0)
            {
                // Bisect to compute the root of F(t) for t >= -e1*e1.
                Real esqr[2] = { e[0]*e[0], e[1]*e[1] };
                Real ey[2] = { e[0]*y[0], e[1]*y[1] };
                Real t0 = -esqr[1] + ey[1];
                Real t1 = -esqr[1] + sqrt(ey[0]*ey[0] + ey[1]*ey[1]);
                Real t = t0;
                const int imax = 2*std::numeric_limits<Real>::max_exponent;
                for (int i = 0; i < imax; ++i)
                {
                    t = ((Real)0.5)*(t0 + t1);
                    if (t == t0 || t == t1)
                    {
                        break;
                    }
                    
                    Real r[2] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]) };
                    Real f = r[0]*r[0] + r[1]*r[1] - (Real)1;
                    if (f > (Real)0)
                    {
                        t0 = t;
                    }
                    else if (f < (Real)0)
                    {
                        t1 = t;
                    }
                    else
                    {
                        break;
                    }
                }
                
                x[0] = esqr[0]*y[0]/(t + esqr[0]);
                x[1] = esqr[1]*y[1]/(t + esqr[1]);
                Real d[2] = { x[0] - y[0], x[1] - y[1] };
                distance = sqrt(d[0]*d[0] + d[1]*d[1]);
            }
            else  // y0 == 0
            {
                x[0] = (Real)0;
                x[1] = e[1];
                distance = fabs(y[1] - e[1]);
            }
        }
        else  // y1 == 0
        {
            Real denom0 = e[0]*e[0] - e[1]*e[1];
            Real e0y0 = e[0]*y[0];
            if (e0y0 < denom0)
            {
                // y0 is inside the subinterval.
                Real x0de0 = e0y0/denom0;
                Real x0de0sqr = x0de0*x0de0;
                x[0] = e[0]*x0de0;
                x[1] = e[1]*sqrt(fabs((Real)1 - x0de0sqr));
                Real d0 = x[0] - y[0];
                distance = sqrt(d0*d0 + x[1]*x[1]);
            }
            else
            {
                // y0 is outside the subinterval.  The closest ellipse point has
                // x1 == 0 and is on the domain-boundary interval (x0/e0)^2 = 1.
                x[0] = e[0];
                x[1] = (Real)0;
                distance = fabs(y[0] - e[0]);
            }
        }
        return distance;
    }
    //----------------------------------------------------------------------------
    // The ellipse is (x0/e0)^2 + (x1/e1)^2 = 1.  The query point is (y0,y1).
    // The function returns the distance from the query point to the ellipse.
    // It also computes the ellipse point (x0,x1) that is closest to (y0,y1).
    //----------------------------------------------------------------------------
    
    static Real DistancePointEllipse (const Real e[2], const Real y[2], Real x[2])
    {
        // Determine reflections for y to the first quadrant.
        bool reflect[2];
        int i, j;
        for (i = 0; i < 2; ++i)
        {
            reflect[i] = (y[i] < (Real)0);
        }
        
        // Determine the axis order for decreasing extents.
        int permute[2];
        if (e[0] < e[1])
        {
            permute[0] = 1;  permute[1] = 0;
        }
        else
        {
            permute[0] = 0;  permute[1] = 1;
        }
        
        int invpermute[2];
        for (i = 0; i < 2; ++i)
        {
            invpermute[permute[i]] = i;
        }
        
        Real locE[2], locY[2];
        for (i = 0; i < 2; ++i)
        {
            j = permute[i];
            locE[i] = e[j];
            locY[i] = y[j];
            if (reflect[j])
            {
                locY[i] = -locY[i];
            }
        }
        
        Real locX[2];
        Real distance = DistancePointEllipseSpecial(locE, locY, locX);
        
        // Restore the axis order and reflections.
        for (i = 0; i < 2; ++i)
        {
            j = invpermute[i];
            if (reflect[j])
            {
                locX[j] = -locX[j];
            }
            x[i] = locX[j];
        }
        
        return distance;
    }
    //----------------------------------------------------------------------------
    // The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1 with e0 >= e1 >= e2.
    // The query point is (y0,y1,y2) with y0 >= 0, y1 >= 0, and y2 >= 0.  The
    // function returns the distance from the query point to the ellipsoid.  It
    // also computes the ellipsoid point (x0,x1,x2) in the first octant that is
    // closest to (y0,y1,y2).
    //----------------------------------------------------------------------------
    
    static Real DistancePointEllipsoidSpecial (const Real e[3], const Real y[3],
                                               Real x[3])
    {
        Real distance;
        if (y[2] > (Real)0)
        {
            if (y[1] > (Real)0)
            {
                if (y[0] > (Real)0)
                {
                    // Bisect to compute the root of F(t) for t >= -e2*e2.
                    Real esqr[3] = { e[0]*e[0], e[1]*e[1], e[2]*e[2] };
                    Real ey[3] = { e[0]*y[0], e[1]*y[1], e[2]*y[2] };
                    Real t0 = -esqr[2] + ey[2];
                    Real t1 = -esqr[2] + sqrt(ey[0]*ey[0] + ey[1]*ey[1] +
                                              ey[2]*ey[2]);
                    Real t = t0;
                    const int imax = 2*std::numeric_limits<Real>::max_exponent;
                    for (int i = 0; i < imax; ++i)
                    {
                        t = ((Real)0.5)*(t0 + t1);
                        if (t == t0 || t == t1)
                        {
                            break;
                        }
                        
                        Real r[3] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]),
                            ey[2]/(t + esqr[2]) };
                        Real f = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - (Real)1;
                        if (f > (Real)0)
                        {
                            t0 = t;
                        }
                        else if (f < (Real)0)
                        {
                            t1 = t;
                        }
                        else
                        {
                            break;
                        }
                    }
                    
                    x[0] = esqr[0]*y[0]/(t + esqr[0]);
                    x[1] = esqr[1]*y[1]/(t + esqr[1]);
                    x[2] = esqr[2]*y[2]/(t + esqr[2]);
                    Real d[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
                    distance = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                }
                else  // y0 == 0
                {
                    x[0] = (Real)0;
                    Real etmp[2] = { e[1], e[2] };
                    Real ytmp[2] = { y[1], y[2] };
                    Real xtmp[2];
                    distance = DistancePointEllipseSpecial(etmp, ytmp, xtmp);
                    x[1] = xtmp[0];
                    x[2] = xtmp[1];
                }
            }
            else  // y1 == 0
            {
                x[1] = (Real)0;
                if (y[0] > (Real)0)
                {
                    Real etmp[2] = { e[0], e[2] };
                    Real ytmp[2] = { y[0], y[2] };
                    Real xtmp[2];
                    distance = DistancePointEllipseSpecial(etmp, ytmp, xtmp);
                    x[0] = xtmp[0];
                    x[2] = xtmp[1];
                }
                else  // y0 == 0
                {
                    x[0] = (Real)0;
                    x[2] = e[2];
                    distance = fabs(y[2] - e[2]);
                }
            }
        }
        else  // y2 == 0
        {
            Real denom[2] = { e[0]*e[0] - e[2]*e[2], e[1]*e[1] - e[2]*e[2] };
            Real ey[2] = { e[0]*y[0], e[1]*y[1] };
            if (ey[0] < denom[0] && ey[1] < denom[1])
            {
                // (y0,y1) is inside the axis-aligned bounding rectangle of the
                // subellipse.  This intermediate test is designed to guard
                // against the division by zero when e0 == e2 or e1 == e2.
                Real xde[2] = { ey[0]/denom[0], ey[1]/denom[1] };
                Real xdesqr[2] = { xde[0]*xde[0], xde[1]*xde[1] };
                Real discr = (Real)1 - xdesqr[0] - xdesqr[1];
                if (discr > (Real)0)
                {
                    // (y0,y1) is inside the subellipse.  The closest ellipsoid
                    // point has x2 > 0.
                    x[0] = e[0]*xde[0];
                    x[1] = e[1]*xde[1];
                    x[2] = e[2]*sqrt(discr);
                    Real d[2] = { x[0] - y[0], x[1] - y[1] };
                    distance = sqrt(d[0]*d[0] + d[1]*d[1] + x[2]*x[2]);
                }
                else
                {
                    // (y0,y1) is outside the subellipse.  The closest ellipsoid
                    // point has x2 == 0 and is on the domain-boundary ellipse
                    // (x0/e0)^2 + (x1/e1)^2 = 1.
                    x[2] = (Real)0;
                    distance = DistancePointEllipseSpecial(e, y, x);
                }
            }
            else
            {
                // (y0,y1) is outside the subellipse.  The closest ellipsoid
                // point has x2 == 0 and is on the domain-boundary ellipse
                // (x0/e0)^2 + (x1/e1)^2 = 1.
                x[2] = (Real)0;
                distance = DistancePointEllipseSpecial(e, y, x);
            }
        }
        return distance;
    }
    //----------------------------------------------------------------------------
    // The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1.  The query point is
    // (y0,y1,y2).  The function returns the distance from the query point to the
    // ellipsoid.   It also computes the ellipsoid point (x0,x1,x2) that is
    // closest to (y0,y1,y2).
    //----------------------------------------------------------------------------
    template <typename Real>
    static Real DistancePointEllipsoid (const Real e[3], const Real y[3], Real x[3])
    {
        // Determine reflections for y to the first octant.
        bool reflect[3];
        int i, j;
        for (i = 0; i < 3; ++i)
        {
            reflect[i] = (y[i] < (Real)0);
        }
        
        // Determine the axis order for decreasing extents.
        int permute[3];
        if (e[0] < e[1])
        {
            if (e[2] < e[0])
            {
                permute[0] = 1;  permute[1] = 0;  permute[2] = 2;
            }
            else if (e[2] < e[1])
            {
                permute[0] = 1;  permute[1] = 2;  permute[2] = 0;
            }
            else
            {
                permute[0] = 2;  permute[1] = 1;  permute[2] = 0;
            }
        }
        else
        {
            if (e[2] < e[1])
            {
                permute[0] = 0;  permute[1] = 1;  permute[2] = 2;
            }
            else if (e[2] < e[0])
            {
                permute[0] = 0;  permute[1] = 2;  permute[2] = 1;
            }
            else
            {
                permute[0] = 2;  permute[1] = 0;  permute[2] = 1;
            }
        }
        
        int invpermute[3];
        for (i = 0; i < 3; ++i)
        {
            invpermute[permute[i]] = i;
        }
        
        Real locE[3], locY[3];
        for (i = 0; i < 3; ++i)
        {
            j = permute[i];
            locE[i] = e[j];
            locY[i] = y[j];
            if (reflect[j])
            {
                locY[i] = -locY[i];
            }
        }
        
        Real locX[3];
        Real distance = DistancePointEllipsoidSpecial(locE, locY, locX);
        
        // Restore the axis order and reflections.
        for (i = 0; i < 3; ++i)
        {
            j = invpermute[i];
            if (reflect[j])
            {
                locX[j] = -locX[j];
            }
            x[i] = locX[j];
        }
        
        return distance;
    }

    struct FillBlocksEllipsoid
    {
        const Real e0,e1,e2,smoothing_length;
        Real position[3];
        Real quaternion[4];
        Real sphere_box[3][2];
        
        void _find_sphere_box()
        {
            const Real maxAxis = std::max(std::max(e0,e1),e2);
            sphere_box[0][0] = position[0] - 2*(maxAxis + smoothing_length);
            sphere_box[0][1] = position[0] + 2*(maxAxis + smoothing_length);
            sphere_box[1][0] = position[1] - 2*(maxAxis + smoothing_length);
            sphere_box[1][1] = position[1] + 2*(maxAxis + smoothing_length);
            sphere_box[2][0] = position[2] - 2*(maxAxis + smoothing_length);
            sphere_box[2][1] = position[2] + 2*(maxAxis + smoothing_length);
        }
        
        FillBlocksEllipsoid(const Real e0, const Real e1, const Real e2, const Real smoothing_length, const Real position[3], const Real quaternion[4]):
        e0(e0),e1(e1),e2(e2), smoothing_length(smoothing_length)
        {
            this->position[0] = position[0];
            this->position[1] = position[1];
            this->position[2] = position[2];
            
            this->quaternion[0] = quaternion[0];
            this->quaternion[1] = quaternion[1];
            this->quaternion[2] = quaternion[2];
            this->quaternion[3] = quaternion[3];
            
            _find_sphere_box();
        }
        
        FillBlocksEllipsoid(const FillBlocksEllipsoid& c):
        e0(c.e0),e1(c.e1),e2(c.e2), smoothing_length(c.smoothing_length)
        {
            position[0] = c.position[0];
            position[1] = c.position[1];
            position[2] = c.position[2];
            
            quaternion[0] = c.quaternion[0];
            quaternion[1] = c.quaternion[1];
            quaternion[2] = c.quaternion[2];
            quaternion[3] = c.quaternion[3];
            
            _find_sphere_box();
        }
        
        Real mollified_heaviside(const Real x, const Real eps) const
        {
            const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
            
            return 0.5+0.5*cos(alpha);
        }
        
        bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
        {
            Real intersection[3][2] = {
                max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
                max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
                max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
            };
            
            return
            intersection[0][1]-intersection[0][0]>0 &&
            intersection[1][1]-intersection[1][0]>0 &&
            intersection[2][1]-intersection[2][0]>0;
            
        }
        
        bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
        {
            Real min_pos[3], max_pos[3];
            
            info.pos(min_pos, 0,0,0);
            info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            for(int i=0;i<3;++i)
            {
                min_pos[i]-=buffer_dx*info.h_gridpoint;
                max_pos[i]+=buffer_dx*info.h_gridpoint;
            }
            return _is_touching(min_pos,max_pos);
        }
        
        inline void operator()(const BlockInfo& info, FluidBlock& b) const
        {
            
            const Real w = quaternion[0];
            const Real x = quaternion[1];
            const Real y = quaternion[2];
            const Real z = quaternion[3];
            
            const Real Rmatrix[3][3] = {
                {1-2*(y*y+z*z),  2*(x*y+z*w),    2*(x*z-y*w)},
                {2*(x*y-z*w),    1-2*(x*x+z*z),  2*(y*z+x*w)},
                {2*(x*z+y*w),    2*(y*z-x*w),    1-2*(x*x+y*y)}
            };
            
            if(_is_touching(info))
            {
                for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                    for(int iy=0; iy<FluidBlock::sizeY; iy++)
                        for(int ix=0; ix<FluidBlock::sizeX; ix++)
                        {
                            Real p[3];
                            info.pos(p, ix, iy, iz);
                            
                            // translate
                            p[0] -= position[0];
                            p[1] -= position[1];
                            p[2] -= position[2];
                            
                            // rotate
                            const Real t[3] = {
                                Rmatrix[0][0]*p[0] + Rmatrix[0][1]*p[1] + Rmatrix[0][2]*p[2],
                                Rmatrix[1][0]*p[0] + Rmatrix[1][1]*p[1] + Rmatrix[1][2]*p[2],
                                Rmatrix[2][0]*p[0] + Rmatrix[2][1]*p[1] + Rmatrix[2][2]*p[2]
                            };
                            
                            // find distance
                            const Real e[3] = {e0,e1,e2};
                            Real xs[3];
                            const Real dist=DistancePointEllipsoid (e, t, xs);
                            const int sign = ( (t[0]*t[0]+t[1]*t[1]+t[2]*t[2]) > (xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2]) ) ? 1 : -1;
                            const Real chi =  mollified_heaviside(sign*dist, smoothing_length);
                            b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
                        }
            }
        }
    };
    
    
    
    struct FillBlocks_Towers
    {
        const Real e0,e1,e2,safe_distance;
        Real position[3];
        Real quaternion[4];
        Real normalI[3], normalJ[3], normalK[3];
        Real rotationMatrix[3][3];
        
        void normalizeNormals()
        {
            const Real magI = std::sqrt(normalI[0]*normalI[0]+normalI[1]*normalI[1]+normalI[2]*normalI[2]);
            const Real magJ = std::sqrt(normalJ[0]*normalJ[0]+normalJ[1]*normalJ[1]+normalJ[2]*normalJ[2]);
            const Real magK = std::sqrt(normalK[0]*normalK[0]+normalK[1]*normalK[1]+normalK[2]*normalK[2]);
            assert(magI > std::numeric_limits<Real>::epsilon());
            assert(magJ > std::numeric_limits<Real>::epsilon());
            assert(magK > std::numeric_limits<Real>::epsilon());
            const Real invMagI = 1.0/magI;
            const Real invMagJ = 1.0/magJ;
            const Real invMagK = 1.0/magK;
            
            for(int i=0;i<3;++i)
            {
                normalI[i]=std::abs(normalI[i])*invMagI;
                normalJ[i]=std::abs(normalJ[i])*invMagJ;
                normalK[i]=std::abs(normalK[i])*invMagK;
            }
        }
        
        void _computeRotationMatrix()
        {
            const Real w = quaternion[0];
            const Real x = quaternion[1];
            const Real y = quaternion[2];
            const Real z = quaternion[3];
            
            const Real Rmatrix[3][3] = {
                {1-2*(y*y+z*z),  2*(x*y+z*w),    2*(x*z-y*w)},
                {2*(x*y-z*w),    1-2*(x*x+z*z),  2*(y*z+x*w)},
                {2*(x*z+y*w),    2*(y*z-x*w),    1-2*(x*x+y*y)}
            };
            memcpy(rotationMatrix, Rmatrix, sizeof(Rmatrix));
            
            // initial axes are grid-aligned: e0 along x, e1 along y, e2 along z
            const Real nx[3] = {1.0,0.0,0.0};
            const Real ny[3] = {0.0,1.0,0.0};
            const Real nz[3] = {0.0,0.0,1.0};
            
            // now we rotate to computational frame
            for(int i=0;i<3;++i)
            {
                normalI[i] = Rmatrix[i][0]*nx[0] + Rmatrix[i][1]*nx[1] + Rmatrix[i][2]*nx[2];
                normalJ[i] = Rmatrix[i][0]*ny[0] + Rmatrix[i][1]*ny[1] + Rmatrix[i][2]*ny[2];
                normalK[i] = Rmatrix[i][0]*nz[0] + Rmatrix[i][1]*nz[1] + Rmatrix[i][2]*nz[2];
            }
            
            normalizeNormals();
        }
        
        FillBlocks_Towers(const Real e0, const Real e1, const Real e2, const Real max_dx, const Real position[3], const Real quaternion[4]):
        e0(e0),e1(e1),e2(e2), safe_distance(max_dx)
        {
            this->position[0] = position[0];
            this->position[1] = position[1];
            this->position[2] = position[2];
            
            this->quaternion[0] = quaternion[0];
            this->quaternion[1] = quaternion[1];
            this->quaternion[2] = quaternion[2];
            this->quaternion[3] = quaternion[3];
            
            _computeRotationMatrix();
        }
        
        FillBlocks_Towers(const FillBlocks_Towers& c):
        e0(c.e0),e1(c.e1),e2(c.e2), safe_distance(c.safe_distance)
        {
            position[0] = c.position[0];
            position[1] = c.position[1];
            position[2] = c.position[2];
            
            quaternion[0] = c.quaternion[0];
            quaternion[1] = c.quaternion[1];
            quaternion[2] = c.quaternion[2];
            quaternion[3] = c.quaternion[3];
            
            _computeRotationMatrix();
        }
        
        bool _is_touching(const Real start[3], const Real end[3]) const
        {
            const Real AABB_w[3] = {
                (Real)0.5*(end[0] - start[0]),
                (Real)0.5*(end[1] - start[1]),
                (Real)0.5*(end[2] - start[2])
            }; // halfwidth
            
            const Real AABB_c[3] = {
                start[0] + AABB_w[0],
                start[1] + AABB_w[1],
                start[2] + AABB_w[2]
            }; // center
            
            assert(AABB_w[0]>0);
            assert(AABB_w[1]>0);
            assert(AABB_w[2]>0);
            
            const Real w[3] = {e0,e1,e2};
            assert(w[0]>0);
            assert(w[1]>0);
            assert(w[2]>0);
            
            bool intersects = true;
            Real r;
            {
                r = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
                intersects &= ((position[0]-r <= AABB_c[0] + AABB_w[0]) && (position[0]+r >= AABB_c[0] - AABB_w[0]));
                
                r = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
                intersects &= ((position[1]-r <= AABB_c[1] + AABB_w[1]) && (position[1]+r >= AABB_c[1] - AABB_w[1]));
                
                r = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
                intersects &= ((position[2]-r <= AABB_c[2] + AABB_w[2]) && (position[2]+r >= AABB_c[2] - AABB_w[2]));
            }
            {
                r = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
                intersects &= ((AABB_c[0]-r <= position[0] + w[0]) && (AABB_c[0]+r >= position[0] - w[0]));
                
                r = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
                intersects &= ((AABB_c[1]-r <= position[1] + w[1]) && (AABB_c[1]+r >= position[1] - w[1]));
                
                r = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
                intersects &= ((AABB_c[2]-r <= position[2] + w[2]) && (AABB_c[2]+r >= position[2] - w[2]));
            }
            return intersects;
        }
        
        bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
        {
            
            const Real safe_distance_info = 2.0*info.h_gridpoint; // two points on each side needed for towers
            
            // AABB - OBB intersection.
            // ellipsoid is inside OBB with normals normalI,normalJ,normalK and widths e0,e1,e2
            // grid block is AABB
            
            Real start[3], end[3];
            
            info.pos(start, 0,0,0);
            info.pos(end, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            
            for(int i=0;i<3;++i)
            {
                start[i]-=(safe_distance_info + buffer_dx*info.h_gridpoint);
                end[i] += (safe_distance_info + buffer_dx*info.h_gridpoint);
            }
            
            return _is_touching(start,end);
        }
        
        
        inline Real distanceToEllipsoid(const Real x, const Real y, const Real z) const
        {
            const Real e[3] = {e0,e1,e2};
            const Real t[3] = {x,y,z};
            Real xs[3];
            const Real dist=DistancePointEllipsoid (e, t, xs);
            const int invSign = ( (t[0]*t[0]+t[1]*t[1]+t[2]*t[2]) > (xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2]) ) ? -1 : 1;
            return dist*invSign;// pos inside, neg outside
        }
        
        Real getHeavisideFDMH1(const Real x, const Real y, const Real z, const Real h) const
        {
            const Real dist = distanceToEllipsoid(x,y,z);
            if(dist >= +h) return 1;
            if(dist <= -h) return 0;
            assert(std::abs(dist)<=h);
            
            const Real distPx = distanceToEllipsoid(x+h,y,z);
            const Real distMx = distanceToEllipsoid(x-h,y,z);
            const Real distPy = distanceToEllipsoid(x,y+h,z);
            const Real distMy = distanceToEllipsoid(x,y-h,z);
            const Real distPz = distanceToEllipsoid(x,y,z+h);
            const Real distMz = distanceToEllipsoid(x,y,z-h);
            
            // compute first primitive of H(x): I(x) = int_0^x H(y) dy and set it to zero outside the cylinder
            const Real IplusX = distPx < 0 ? 0 : distPx;
            const Real IminuX = distMx < 0 ? 0 : distMx;
            const Real IplusY = distPy < 0 ? 0 : distPy;
            const Real IminuY = distMy < 0 ? 0 : distMy;
            const Real IplusZ = distPz < 0 ? 0 : distPz;
            const Real IminuZ = distMz < 0 ? 0 : distMz;
            
            assert(IplusX>=0);
            assert(IminuX>=0);
            assert(IplusY>=0);
            assert(IminuY>=0);
            assert(IplusZ>=0);
            assert(IminuZ>=0);
            
            
            // gradI
            const Real gradIX = 0.5/h * (IplusX - IminuX);
            const Real gradIY = 0.5/h * (IplusY - IminuY);
            const Real gradIZ = 0.5/h * (IplusZ - IminuZ);
            
            // gradU
            const Real gradUX = 0.5/h * (distPx - distMx);
            const Real gradUY = 0.5/h * (distPy - distMy);
            const Real gradUZ = 0.5/h * (distPz - distMz);
            
            const Real H = (gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ)/(gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ);
            
            //        assert(H>=0 && H<=1);
            if(H<0.0 || H>1.0)
                printf("invalid H?: %10.10e %10.10e %10.10e: %10.10e\n",x,y,z,H);
            
            return H;
        }
        
        
        inline void operator()(const BlockInfo& info, FluidBlock& b) const
        {
            if(_is_touching(info))
            {
                for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                    for(int iy=0; iy<FluidBlock::sizeY; iy++)
                        for(int ix=0; ix<FluidBlock::sizeX; ix++)
                        {
                            Real p[3];
                            info.pos(p, ix, iy, iz);
                            
                            // translate
                            p[0] -= position[0];
                            p[1] -= position[1];
                            p[2] -= position[2];
                            
                            // rotate
                            const Real t[3] = {
                                rotationMatrix[0][0]*p[0] + rotationMatrix[0][1]*p[1] + rotationMatrix[0][2]*p[2],
                                rotationMatrix[1][0]*p[0] + rotationMatrix[1][1]*p[1] + rotationMatrix[1][2]*p[2],
                                rotationMatrix[2][0]*p[0] + rotationMatrix[2][1]*p[1] + rotationMatrix[2][2]*p[2]
                            };

                            const Real chi = getHeavisideFDMH1(t[0],t[1],t[2],info.h_gridpoint);
                            b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
                        }
            }
        }
    };

}

namespace VAWTObstacle
{
    struct FillBlocks_Towers
    {
        const Real radius,height,chord,thickness; // chord and thickness are actually haflchord and halfthickness
        const Real safe_distance;
        const int nBlades=3;
        Real position[3];
        Real rotationAngle;
        Real bounding_box[3][2];
        const Real pitchAngle = 45.0*M_PI/180;
        
        void _find_bounding_box()
        {
            const Real width = thickness + chord*std::sin(pitchAngle);
            for(int i=0;i<2;++i)
            {
                bounding_box[i][0] = position[i] - radius - width - safe_distance;
                bounding_box[i][1] = position[i] + radius + width + safe_distance;
            }
            bounding_box[2][0] = position[2] - 0.5*height - safe_distance;
            bounding_box[2][1] = position[2] + 0.5*height + safe_distance;
        }
        
        FillBlocks_Towers(const Real radius, const Real height, const Real chord, const Real thickness, const Real max_dx, const Real position[3], const Real rotationAngle):
        radius(radius),height(height),chord(chord),thickness(thickness), safe_distance(max_dx),rotationAngle(rotationAngle)
        {
            this->position[0] = position[0];
            this->position[1] = position[1];
            this->position[2] = position[2];
            
            _find_bounding_box();
        }
        
        FillBlocks_Towers(const FillBlocks_Towers& c):
        radius(c.radius),height(c.height),chord(c.chord),thickness(c.thickness), safe_distance(c.safe_distance),rotationAngle(c.rotationAngle),nBlades(c.nBlades),pitchAngle(c.pitchAngle)
        {
            position[0] = c.position[0];
            position[1] = c.position[1];
            position[2] = c.position[2];
            
            _find_bounding_box();
        }
        
        bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
        {
            Real intersection[3][2] = {
                max(min_pos[0], bounding_box[0][0]), min(max_pos[0], bounding_box[0][1]),
                max(min_pos[1], bounding_box[1][0]), min(max_pos[1], bounding_box[1][1]),
                max(min_pos[2], bounding_box[2][0]), min(max_pos[2], bounding_box[2][1])
            };
            
            return
            intersection[0][1]-intersection[0][0]>0 &&
            intersection[1][1]-intersection[1][0]>0 &&
            intersection[2][1]-intersection[2][0]>0;
        }
        
        bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
        {
            Real min_pos[3], max_pos[3];
            
            info.pos(min_pos, 0,0,0);
            info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            for(int i=0;i<3;++i)
            {
                min_pos[i]-=buffer_dx*info.h_gridpoint;
                max_pos[i]+=buffer_dx*info.h_gridpoint;
            }
            return _is_touching(min_pos,max_pos);
        }
        
        inline Real distanceToEllipse(const Real x, const Real y, const Real z, const int bladeIdx) const
        {
            assert(nBlades==3);
            const Real deltaBladeAngle = 2.0*M_PI/((Real)nBlades);
            const Real angle = rotationAngle + bladeIdx*deltaBladeAngle + pitchAngle;
            
            Real t[2] = {  x*std::cos(angle) + y*std::sin(angle),
                -x*std::sin(angle) + y*std::cos(angle)};
            // wtf
            if(std::abs(t[0])<std::numeric_limits<Real>::epsilon()) t[0]=0;
            if(std::abs(t[1])<std::numeric_limits<Real>::epsilon()) t[1]=0;
            
            const Real e[2] = {thickness,chord};
            
            Real xs[2];
            const Real dist=EllipsoidObstacle::DistancePointEllipse (e, t, xs);
            const int invSign = ( (t[0]*t[0]+t[1]*t[1]) > (xs[0]*xs[0]+xs[1]*xs[1]) ) ? -1 : 1;
            
            if(std::abs(z)<=0.5*height)
            {
                return dist*invSign;// pos inside, neg outside
            }
            else
            {
                return -std::min<Real>(std::sqrt(dist*dist + z*z), std::abs(z)-0.5*height); // always negative
            }
        }
        
        Real getHeavisideFDMH1(const Real x, const Real y, const Real z, const Real h, const int bladeIdx) const
        {
            const Real dist = distanceToEllipse(x,y,z,bladeIdx);
            if(dist >= +h) return 1;
            if(dist <= -h) return 0;
            assert(std::abs(dist)<=h);
            
            const Real distPx = distanceToEllipse(x+h,y,z,bladeIdx);
            const Real distMx = distanceToEllipse(x-h,y,z,bladeIdx);
            const Real distPy = distanceToEllipse(x,y+h,z,bladeIdx);
            const Real distMy = distanceToEllipse(x,y-h,z,bladeIdx);
            const Real distPz = distanceToEllipse(x,y,z+h,bladeIdx);
            const Real distMz = distanceToEllipse(x,y,z-h,bladeIdx);
            
            // compute first primitive of H(x): I(x) = int_0^x H(y) dy and set it to zero outside the cylinder
            const Real IplusX = distPx < 0 ? 0 : distPx;
            const Real IminuX = distMx < 0 ? 0 : distMx;
            const Real IplusY = distPy < 0 ? 0 : distPy;
            const Real IminuY = distMy < 0 ? 0 : distMy;
            const Real IplusZ = distPz < 0 ? 0 : distPz;
            const Real IminuZ = distMz < 0 ? 0 : distMz;
            
            assert(IplusX>=0);
            assert(IminuX>=0);
            assert(IplusY>=0);
            assert(IminuY>=0);
            assert(IplusZ>=0);
            assert(IminuZ>=0);
            
            
            // gradI
            const Real gradIX = 0.5/h * (IplusX - IminuX);
            const Real gradIY = 0.5/h * (IplusY - IminuY);
            const Real gradIZ = 0.5/h * (IplusZ - IminuZ);
            
            // gradU
            const Real gradUX = 0.5/h * (distPx - distMx);
            const Real gradUY = 0.5/h * (distPy - distMy);
            const Real gradUZ = 0.5/h * (distPz - distMz);
            
            const Real H = (gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ)/(gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ);
            
            //        assert(H>=0 && H<=1);
            if(H<0.0 || H>1.0)
                printf("invalid H?: %10.10e %10.10e %10.10e: %10.10e\n",x,y,z,H);
            
            return H;
        }
        
        inline void operator()(const BlockInfo& info, FluidBlock& b) const
        {
            if(_is_touching(info))
            {
                std::vector<Real> bladePosX(nBlades),bladePosY(nBlades);
                const Real deltaBladeAngle = 2.0*M_PI/((Real)nBlades);
                
                for(int i=0;i<nBlades;++i)
                {
                    const Real bladeAngle = rotationAngle + i*deltaBladeAngle;
                    bladePosX[i] = radius * std::cos(bladeAngle);
                    bladePosY[i] = radius * std::sin(bladeAngle);
                }
                
                for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                    for(int iy=0; iy<FluidBlock::sizeY; iy++)
                        for(int ix=0; ix<FluidBlock::sizeX; ix++)
                        {
                            Real p[3];
                            info.pos(p, ix, iy, iz);
                            
                            // translate
                            p[0] -= position[0];
                            p[1] -= position[1];
                            p[2] -= position[2];
                            
                            // find blade Idx
                            Real minDistSq=2.0;
                            int minIdx=-1;
                            for(int i=0;i<nBlades;++i)
                            {
                                const Real distSq = (p[0]-bladePosX[i])*(p[0]-bladePosX[i]) + (p[1]-bladePosY[i])*(p[1]-bladePosY[i]);
                                minIdx = (distSq<minDistSq) ? i : minIdx;
                                minDistSq = (distSq<minDistSq) ? distSq : minDistSq;
                            }
                            assert(minIdx>=0 && minIdx < nBlades);
                            
                            p[0]-=bladePosX[minIdx];
                            p[1]-=bladePosY[minIdx];

                            const Real chi = getHeavisideFDMH1(p[0],p[1],p[2],info.h_gridpoint,minIdx);
                            b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
                        }
            }
        }
    };
}

namespace Schedulers
{
    template<int Npoints>
    struct ParameterScheduler
    {
        std::array<Real, Npoints> parameters_t0; // parameters at t0
        std::array<Real, Npoints> parameters_t1; // parameters at t1
        std::array<Real, Npoints> dparameters_t0; // derivative at t0
        Real t0, t1; // t0 and t1

        void save(std::string filename)
        {
            std::ofstream savestream;
            savestream.setf(std::ios::scientific);
            savestream.precision(std::numeric_limits<Real>::digits10 + 1);
            savestream.open(filename+".txt");

            savestream << t0 << "\t" << t1 << std::endl;
            for(int i=0;i<Npoints;++i)
                savestream << parameters_t0[i] << "\t" << parameters_t1[i] << "\t" << dparameters_t0[i] << std::endl;
            savestream.close();
        }

        void restart(std::string filename)
        {
            std::ifstream restartstream;
            restartstream.open(filename+".txt");

            restartstream >> t0 >> t1;
            for(int i=0;i<Npoints;++i)
                restartstream >> parameters_t0[i] >> parameters_t1[i] >> dparameters_t0[i];
            restartstream.close();
        }

        ParameterScheduler()
        {
            t0=-1;
            parameters_t0 = std::array<Real, Npoints>();
            dparameters_t0 = std::array<Real, Npoints>();
        }

        void transition(const Real t, const Real tstart, const Real tend,
        		const std::array<Real, Npoints> parameters_tend,
				const bool UseCurrentDerivative = false)
        {
            if(t<tstart or t>tend) return; // this transition is out of scope
            if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

            // we transition from whatever state we are in to a new state
            // the start point is where we are now: lets find out
            std::array<Real, Npoints> parameters;
            std::array<Real, Npoints> dparameters;
            gimmeValues(tstart,parameters,dparameters);

            /*
            if (Npoints == 7)
                printf("[t0 t1] were [%f %f], will be [%f %f], parameters %f dpdt %f param end %f\n",
                		t0,t1,tstart,tend,parameters[3],dparameters[3], parameters_tend[3]);
            */

            // fill my members
            t0 = tstart;
            t1 = tend;
            parameters_t0 = parameters;
            parameters_t1 = parameters_tend;
            dparameters_t0 = UseCurrentDerivative ? dparameters : std::array<Real, Npoints>();
        }

        void transition(const Real t, const Real tstart, const Real tend,
        		const std::array<Real, Npoints> parameters_tstart,
				const std::array<Real, Npoints> parameters_tend)
        {
            if(t<tstart or t>tend) return; // this transition is out of scope
            if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

            // fill my members
            t0 = tstart;
            t1 = tend;
            parameters_t0 = parameters_tstart;
            parameters_t1 = parameters_tend;
        }

        void gimmeValues(const Real t,std::array<Real, Npoints> & parameters, std::array<Real, Npoints> & dparameters)
        {
            // look at the different cases
            if(t<t0 or t0<0) { // no transition, we are in state 0
                parameters = parameters_t0;
                dparameters = std::array<Real, Npoints>();
            } else if(t>t1) { // no transition, we are in state 1
                parameters = parameters_t1;
                dparameters = std::array<Real, Npoints>();
            } else { // we are within transition: interpolate
                for(int i=0;i<Npoints;++i)
                    IF2D_Interpolation1D::cubicInterpolation(t0,t1,t,parameters_t0[i],parameters_t1[i],dparameters_t0[i],0.0,parameters[i],dparameters[i]);
            }
        }

        void gimmeValues(const Real t,std::array<Real, Npoints> & parameters)
        {
            std::array<Real, Npoints> dparameters_whocares; // no derivative info
            return gimmeValues(t,parameters,dparameters_whocares);
        }
    };

    struct ParameterSchedulerScalar : ParameterScheduler<1>
    {
        void transition(const Real t, const Real tstart, const Real tend, const Real parameter_tend, const bool UseCurrentDerivative = false)
        {
            const std::array<Real, 1> myParameter = {parameter_tend};
            return ParameterScheduler<1>::transition(t,tstart,tend,myParameter,UseCurrentDerivative);
        }

        void gimmeValues(const Real t, Real & parameter, Real & dparameter)
        {
            std::array<Real, 1> myParameter, mydParameter;
            ParameterScheduler<1>::gimmeValues(t, myParameter, mydParameter);
            parameter = myParameter[0];
            dparameter = mydParameter[0];
        }

        void gimmeValues(const Real t, Real & parameter)
        {
            std::array<Real, 1> myParameter;
            ParameterScheduler<1>::gimmeValues(t, myParameter);
            parameter = myParameter[0];
        }
    };

    template<int Npoints>
    struct ParameterSchedulerVector : ParameterScheduler<Npoints>
    {
        void gimmeValues(const Real t, const std::array<Real, Npoints> & positions, const int Nfine,
        		const Real * const positions_fine, Real * const parameters_fine, Real * const dparameters_fine)
        {
            // we interpolate in space the start and end point
            Real parameters_t0_fine[Nfine];
            Real parameters_t1_fine[Nfine];
            Real dparameters_t0_fine[Nfine];

            IF2D_Interpolation1D::naturalCubicSpline(positions.data(), this->parameters_t0.data(), Npoints, positions_fine, parameters_t0_fine,  Nfine);
            IF2D_Interpolation1D::naturalCubicSpline(positions.data(), this->parameters_t1.data(), Npoints, positions_fine, parameters_t1_fine,  Nfine);
            IF2D_Interpolation1D::naturalCubicSpline(positions.data(), this->dparameters_t0.data(),Npoints, positions_fine, dparameters_t0_fine, Nfine);

            // look at the different cases
            if(t<this->t0 or this->t0<0) {
                // no transition, we are in state 0
                for(int i=0;i<Nfine;++i) {
                    parameters_fine[i] = parameters_t0_fine[i];
                    dparameters_fine[i] = 0.0;
                }
            } else if(t>this->t1) {
                // no transition, we are in state 1
                for(int i=0;i<Nfine;++i) {
                    parameters_fine[i] = parameters_t1_fine[i];
                    dparameters_fine[i] = 0.0;
                }
            } else {
                // we are within transition: interpolate in time for each point of the fine discretization
                for(int i=0;i<Nfine;++i)
                    IF2D_Interpolation1D::cubicInterpolation(this->t0,this->t1,t,parameters_t0_fine[i],parameters_t1_fine[i],dparameters_t0_fine[i],
                    											  	  	  	  	 0.0,				   parameters_fine[i],	 dparameters_fine[i]);
            }
        }

        void gimmeValues(const Real t,std::array<Real, Npoints> & parameters)
        {
            ParameterScheduler<Npoints>::gimmeValues(t, parameters);
        }

        void gimmeValues(const Real t,std::array<Real, Npoints> & parameters, std::array<Real, Npoints> & dparameters)
        {
            ParameterScheduler<Npoints>::gimmeValues(t, parameters, dparameters);
        }
    };

    template<int Npoints>
    struct ParameterSchedulerLearnWave : ParameterScheduler<Npoints>
    {
    	void gimmeValues(const Real t, const Real Twave, const Real Length, const std::array<Real, Npoints> & positions,
    					const int Nfine, const Real* const positions_fine, Real* const parameters_fine, Real* const dparameters_fine)
        {
    		const Real _1oL = 1./Length;
    		const Real _1oT = 1./Twave;
            // the fish goes through (as function of t and s) a wave function that describes the curvature
            for(int i=0;i<Nfine;++i) {
                const Real c = positions_fine[i]*_1oL - (t - this->t0)*_1oT; //traveling wave coord
            	bool bCheck = true;

                if (c < positions[0]) { // Are you before latest wave node?
                    IF2D_Interpolation1D::cubicInterpolation(c, positions[0], c,
                    		this->parameters_t0[0], this->parameters_t0[0],
							parameters_fine[i], dparameters_fine[i]);
                    bCheck = false;
                }
                else if (c > positions[Npoints-1]) {// Are you after oldest wave node?
                    IF2D_Interpolation1D::cubicInterpolation(positions[Npoints-1], c, c,
                    		this->parameters_t0[Npoints-1], this->parameters_t0[Npoints-1],
							parameters_fine[i], dparameters_fine[i]);
                    bCheck = false;
                } else {
                    for (int j=1; j<Npoints; ++j) { // Check at which point of the travelling wave we are
                        if (( c >= positions[j-1] ) && ( c <= positions[j] )) {
                            IF2D_Interpolation1D::cubicInterpolation(positions[j-1], positions[j], c,
                            		this->parameters_t0[j-1], this->parameters_t0[j],
									parameters_fine[i], dparameters_fine[i]);
                            dparameters_fine[i] = -dparameters_fine[i]*_1oT; // df/dc * dc/dt
                            bCheck = false;
                        }
                    }
                }
                if (bCheck) { std::cout << "Ciaone2!" << std::endl; abort(); }
            }
        }

        void Turn(const Real b, const Real t_turn) // each decision adds a node at the beginning of the wave (left, right, straight) and pops last node
        {
        	this->t0 = t_turn;
            for(int i=Npoints-1; i>0; --i) this->parameters_t0[i] = this->parameters_t0[i-2];
            this->parameters_t0[1] = b;
            this->parameters_t0[0] = 0;
        }
    };
}

#endif
