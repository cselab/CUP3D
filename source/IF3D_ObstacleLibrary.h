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
#include "GenericOperator.h"

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
		radius(radius),safe_radius(radius+4*max_dx), sphere_position{pos[0],pos[1],pos[2]}
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

	inline Real sign(const Real& val) const {
		return (0. < val) - (val < 0.);
	}
#if 1
inline void operator()(const BlockInfo& info,FluidBlock& b,ObstacleBlock* const defblock,surfaceBlocks* const surf) const
{
	if(_is_touching(info)) {
		const Real eps = std::numeric_limits<Real>::epsilon();
		const Real h = info.h_gridpoint;
		const Real invh = 1./h;
		const Real fac1 = .5/h;
		const Real fac2 = .5/(std::sqrt(2.)*h);
		const Real fac[9] = {fac1,fac1,fac1,fac2,fac2,fac2,fac2,fac2,fac2};

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++) {
					Real p[3];
					info.pos(p, ix, iy, iz);
					const Real x = p[0]-sphere_position[0];
					const Real y = p[1]-sphere_position[1];
					const Real z = p[2]-sphere_position[2];
					const Real U = distanceToSphere(x,y,z);
					if(U > 2*h || U < -2*h) { //2 should be safe
						const Real H = U > 0 ? 1.0 : 0.0;
						defblock->chi[iz][iy][ix] = H;
						b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
						continue;
					}
					const Real Up[9] = {
							distanceToSphere(x+h,y,z),   distanceToSphere(x,y+h,z),   distanceToSphere(x,y,z+h),
							distanceToSphere(x,y+h,z+h), distanceToSphere(x+h,y,z+h), distanceToSphere(x+h,y+h,z),
							distanceToSphere(x,y+h,z-h), distanceToSphere(x+h,y,z-h), distanceToSphere(x+h,y-h,z)
					};
					const Real Um[9] = {
							distanceToSphere(x-h,y,z),   distanceToSphere(x,y-h,z),   distanceToSphere(x,y,z-h),
							distanceToSphere(x,y-h,z-h), distanceToSphere(x-h,y,z-h), distanceToSphere(x-h,y-h,z),
							distanceToSphere(x,y-h,z+h), distanceToSphere(x-h,y,z+h), distanceToSphere(x-h,y+h,z)
					};
					Real Ip[9],Im[9];
					for (int i=0; i<9; i++) {
						Ip[i] = std::max(Up[i],(Real)0.);
						Im[i] = std::max(Um[i],(Real)0.);
					}
					Real gradU[9], gradI[9];
					for (int i=0; i<9; i++) {
						gradU[i] = fac[i]*(Up[i]-Um[i]);
						gradI[i] = fac[i]*(Ip[i]-Im[i]);
					}
					Real Hp[3], Hm[3];
					for(int i=0; i<3; i++) {
						Hp[i] = (Up[i]> h) ? h : (
								(Up[i]<-h) ? 0 :
										.5*h+(Up[i]-fac1*sign(Up[i])*Up[i]*Up[i]));
						Hm[i] = (Um[i]> h) ? h : (
								(Um[i]<-h) ? 0 :
										.5*h+(Um[i]-fac1*sign(Um[i])*Um[i]*Um[i]));
					}
					const Real gradH[3] = {.5*(Hp[0]-Hm[0]), .5*(Hp[1]-Hm[1]), .5*(Hp[2]-Hm[2])};
					Real gradUU[3], gradUI[3], gradUH[3];
					for (int i=0; i<3; i++) {
						gradUU[i] = gradU[i]*gradU[i] + gradU[i+3]*gradU[i+3] + gradU[i+6]*gradU[i+6];
						gradUI[i] = gradU[i]*gradI[i] + gradU[i+3]*gradI[i+3] + gradU[i+6]*gradI[i+6];
						gradUH[i] = gradU[i]*gradH[i];
					}
					for (int i=0; i<3; i++)  gradUU[i] = max(gradUU[i], eps);
					const Real FDD = h*(gradUH[0] + gradUH[1] + gradUH[2])/gradUU[0];
					const Real FDH = 1/3. * (gradUI[0]/gradUU[0]+gradUI[1]/gradUU[1]+gradUI[2]/gradUU[2]);

					if (FDD>1e-6) {
						const Real dchidx = -FDD*gradU[0];
						const Real dchidy = -FDD*gradU[1];
						const Real dchidz = -FDD*gradU[2];
						surf->add(info.blockID,ix,iy,iz,dchidx,dchidy,dchidz,FDD);
					}
#ifndef NDEBUG
					if(FDH<0 || FDH>1) printf("invalid H?: %9.9e %9.9e %9.9e: %9.9e\n",x,y,z,FDH);
#endif
					defblock->chi[iz][iy][ix] = FDH;
					b(ix,iy,iz).chi = std::max(FDH, b(ix,iy,iz).chi);
					//b(ix,iy,iz).tmpU = dist;
					//b(ix,iy,iz).tmpV = Delta;
				}
	}
}

#else
inline void operator()(const BlockInfo& info,FluidBlock& b,ObstacleBlock* const defblock,surfaceBlocks* const surf) const
{
	if(_is_touching(info)) {
		const Real h = info.h_gridpoint;
		const Real fac1 = 0.5/h;
		const Real eps = std::numeric_limits<Real>::epsilon();

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++) {
					Real p[3];
					info.pos(p, ix, iy, iz);
					const Real x = p[0]-sphere_position[0];
					const Real y = p[1]-sphere_position[1];
					const Real z = p[2]-sphere_position[2];
					const Real dist = distanceToSphere(x,y,z);

					if(dist > 2*h || dist < -2*h) { //2 should be safe
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
					const Real gradUX = inv2h * (distPx - distMx);
					const Real gradUY = inv2h * (distPy - distMy);
					const Real gradUZ = inv2h * (distPz - distMz);
					const Real gradIX = inv2h * (IplusX - IminuX);
					const Real gradIY = inv2h * (IplusY - IminuY);
					const Real gradIZ = inv2h * (IplusZ - IminuZ);
					const Real gradHX = inv2h * (HplusX - HminuX);
					const Real gradHY = inv2h * (HplusY - HminuY);
					const Real gradHZ = inv2h * (HplusZ - HminuZ);

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
#ifndef NDEBUG
					if(H<0 || H>1)
						printf("invalid H?: %10.10e %10.10e %10.10e: %10.10e\n",x,y,z,H);
#endif
					defblock->chi[iz][iy][ix] = H;
					b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
					//b(ix,iy,iz).tmpU = dist;
					//b(ix,iy,iz).tmpV = Delta;
				}
	}
}
#endif
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

namespace FishObstacle
{
struct VolumeSegment_OBB
{
	std::pair<int, int> s_range;
	Real normalI[3]; // should be normalized and >=0
	Real normalJ[3];
	Real normalK[3];
	Real w[3]; // halfwidth
	Real c[3]; // center

#if 1==0
	VolumeSegment_OBB(std::pair<int, int> s_range, const Real bbox[3][2])
	: s_range(s_range)
	{
		normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
		normalI[0]=normalJ[1]=normalK[2]=1.0;
		for(int i=0;i<3;++i) {
			w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
			c[i] = bbox[i][0] + w[i];
			assert(w[i]>0);
		}

	}
#endif

	VolumeSegment_OBB() { }

	void prepare(std::pair<int, int> _s_range, const Real bbox[3][2])
	{
		normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
		normalI[0]=normalJ[1]=normalK[2]=1.0;
		s_range.first = _s_range.first;
		s_range.second = _s_range.second;
		for(int i=0;i<3;++i) {
			w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
			c[i] = bbox[i][0] + w[i];
			assert(w[i]>0);
		}
	}

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

		for(int i=0;i<3;++i) {
			// also take absolute value since thats what we need when doing intersection checks later
			normalI[i]=std::abs(normalI[i])*invMagI;
			normalJ[i]=std::abs(normalJ[i])*invMagJ;
			normalK[i]=std::abs(normalK[i])*invMagK;
		}
	}

	void changeToComputationalFrame(const Real position[3], const Real quaternion[4])
	{
		// we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
		const Real w = quaternion[0];
		const Real x = quaternion[1];
		const Real y = quaternion[2];
		const Real z = quaternion[3];
		const Real Rmatrix3D[3][3] = {
				{1.-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
				{2*(x*y+z*w),    1.-2*(x*x+z*z),  2*(y*z-x*w)},
				{2*(x*z-y*w),    2*(y*z+x*w),    1.-2*(x*x+y*y)}
		};
		const Real p[3] = {c[0],c[1],c[2]};
		const Real nx[3] = {normalI[0],normalI[1],normalI[2]};
		const Real ny[3] = {normalJ[0],normalJ[1],normalJ[2]};
		const Real nz[3] = {normalK[0],normalK[1],normalK[2]};
		for(int i=0;i<3;++i) {
			c[i] = Rmatrix3D[i][0]*p[0] + Rmatrix3D[i][1]*p[1] + Rmatrix3D[i][2]*p[2];
			normalI[i] = Rmatrix3D[i][0]*nx[0] + Rmatrix3D[i][1]*nx[1] + Rmatrix3D[i][2]*nx[2];
			normalJ[i] = Rmatrix3D[i][0]*ny[0] + Rmatrix3D[i][1]*ny[1] + Rmatrix3D[i][2]*ny[2];
			normalK[i] = Rmatrix3D[i][0]*nz[0] + Rmatrix3D[i][1]*nz[1] + Rmatrix3D[i][2]*nz[2];
		}
		c[0] +=position[0];
		c[1] +=position[1];
		c[2] +=position[2];
		normalizeNormals();
	}

	bool isIntersectingWithAABB(const Real start[3],const Real end[3], const Real safe_distance = 0.0) const
	{
		//start and end are two diagonally opposed corners of grid block
		const Real AABB_w[3] = { //half block width + safe distance
				0.5*(end[0] - start[0]) + 2.0*safe_distance,
				0.5*(end[1] - start[1]) + 2.0*safe_distance,
				0.5*(end[2] - start[2]) + 2.0*safe_distance
		};
		const Real AABB_c[3] = { //block center
				start[0] + AABB_w[0] - safe_distance,
				start[1] + AABB_w[1] - safe_distance,
				start[2] + AABB_w[2] - safe_distance
		};
		assert(AABB_w[0]>0);
		assert(AABB_w[1]>0);
		assert(AABB_w[2]>0);
		bool intersects = true;
		Real r;
		{
			r = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
			intersects &= ((c[0]-r <= AABB_c[0] + AABB_w[0]) && (c[0]+r >= AABB_c[0] - AABB_w[0]));

			r = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
			intersects &= ((c[1]-r <= AABB_c[1] + AABB_w[1]) && (c[1]+r >= AABB_c[1] - AABB_w[1]));

			r = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
			intersects &= ((c[2]-r <= AABB_c[2] + AABB_w[2]) && (c[2]+r >= AABB_c[2] - AABB_w[2]));
		}
		{
			//r = AABB_w[0]*normalI[0] + AABB_w[1]*normalI[1] + AABB_w[2]*normalI[2];
			r = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
			intersects &= ((AABB_c[0]-r <= c[0] + w[0]) && (AABB_c[0]+r >= c[0] - w[0]));

			//r = AABB_w[0]*normalJ[0] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalJ[2];
			r = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
			intersects &= ((AABB_c[1]-r <= c[1] + w[1]) && (AABB_c[1]+r >= c[1] - w[1]));

			//r = AABB_w[0]*normalK[0] + AABB_w[1]*normalK[1] + AABB_w[2]*normalK[2];
			r = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
			intersects &= ((AABB_c[2]-r <= c[2] + w[2]) && (AABB_c[2]+r >= c[2] - w[2]));
		}
		return intersects;
	}
};

struct PutFishOnBlocks
{
	const Fish::FishMidlineData * cfish;
	Real position[3];
	Real quaternion[4];
	Real Rmatrix3D[3][3];

	PutFishOnBlocks(const Fish::FishMidlineData* const cfish, const Real pos[3], const Real quat[4]):
		cfish(cfish)
	{
		position[0]=pos[0];
		position[1]=pos[1];
		position[2]=pos[2];
		quaternion[0]=quat[0];
		quaternion[1]=quat[1];
		quaternion[2]=quat[2];
		quaternion[3]=quat[3];
		computeRotationMatrix();
	}

	void computeRotationMatrix()
	{
		const Real w = quaternion[0];
		const Real x = quaternion[1];
		const Real y = quaternion[2];
		const Real z = quaternion[3];
		const Real R[3][3] = {
				{1-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
				{2*(x*y+z*w),    1-2*(x*x+z*z),  2*(y*z-x*w)},
				{2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x*x+y*y)}
		};
		memcpy(Rmatrix3D, R, sizeof(R));
	}

	inline int find_closest_dist(const int s, const int dir, const Real x[3], Real & oldDistSq) const
	{
		if((s+dir)<cfish->iFishStart or (s+dir)>cfish->iFishEnd)
			return s;

		const Real newDistSq = (x[0]-cfish->rX[s+dir])*(x[0]-cfish->rX[s+dir])
        		            																				   + (x[1]-cfish->rY[s+dir])*(x[1]-cfish->rY[s+dir])
																											   + (x[2])*(x[2]);

		if(oldDistSq<=newDistSq) {
			return s;
		} else {
			oldDistSq = newDistSq;
			return s+dir;
		}
	}

	void changeVelocityToComputationalFrame(Real x[3]) const
	{
		const Real p[3] = {x[0],x[1],x[2]};

		// rotate (around CoM)
		x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[0][1]*p[1] + Rmatrix3D[0][2]*p[2];
		x[1]=Rmatrix3D[1][0]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[1][2]*p[2];
		x[2]=Rmatrix3D[2][0]*p[0] + Rmatrix3D[2][1]*p[1] + Rmatrix3D[2][2]*p[2];
	}

	template<typename T>
	void changeToComputationalFrame(T x[3]) const
	{
		const T p[3] = {x[0],x[1],x[2]};

		// rotate (around CoM)
		x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[0][1]*p[1] + Rmatrix3D[0][2]*p[2];
		x[1]=Rmatrix3D[1][0]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[1][2]*p[2];
		x[2]=Rmatrix3D[2][0]*p[0] + Rmatrix3D[2][1]*p[1] + Rmatrix3D[2][2]*p[2];

		// translate
		x[0]+=position[0];
		x[1]+=position[1];
		x[2]+=position[2];
	}

	template<typename T>
	void changeFromComputationalFrame(T x[3]) const
	{
		// translate back to CoM
		const T p[3] = {
				x[0]-position[0],
				x[1]-position[1],
				x[2]-position[2]
		};

		// rotate back around CoM
		x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[1][0]*p[1] + Rmatrix3D[2][0]*p[2];
		x[1]=Rmatrix3D[0][1]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[2][1]*p[2];
		x[2]=Rmatrix3D[0][2]*p[0] + Rmatrix3D[1][2]*p[1] + Rmatrix3D[2][2]*p[2];
	}

	Real getSmallerDistToMidline(const int start_s, const Real x[3], int & final_s) const
	{
		Real relX[3] = {x[0],x[1],x[2]};
		changeFromComputationalFrame(relX);

		const Real curDistSq =  std::pow(relX[0]-cfish->rX[start_s],2)
		+ std::pow(relX[1]-cfish->rY[start_s],2)
		+ std::pow(relX[2],2);

		Real distSq;

		distSq = curDistSq; // check right
		const int sRight = find_closest_dist(start_s, +1, relX, distSq);

		distSq = curDistSq; // check left
		const int sLeft = find_closest_dist(start_s, -1, relX, distSq);

		if(sRight==start_s and sLeft==start_s) {
			final_s = start_s;
			return distSq;
		}

		assert(sRight==start_s or sLeft==start_s);

		int curr_s = start_s;
		int new_s = sRight == start_s ? sLeft : sRight;
		const int dir = new_s-curr_s;
		while(curr_s not_eq new_s) {
			curr_s = new_s;
			new_s = find_closest_dist(curr_s,dir,relX, distSq);
		}

		final_s = new_s;
		return distSq;
	}

	//inline void operator()(const BlockInfo& info, FluidBlock3D& b) const
	void operator()(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
	{
		Real org[3];
		info.pos(org, 0, 0, 0);
		const Real invh = 1.0/info.h_gridpoint;

		{
			Real * const chi = &(defblock->chi[0][0][0]);
			Real * const udef = &(defblock->udef[0][0][0][0]);

			static const int n = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX;
			for(int i=0; i<n; i++) {
				chi[i]=-1;
				udef[3*i+0]=0;
				udef[3*i+1]=0;
				udef[3*i+2]=0;
			}
		}

		// construct the shape (P2M with min(distance) as kernel) onto defblocks
		for(int i=0;i<vSegments.size();++i) {
			//iterate over segments contained in the vSegm intersecting this block:
			const int firstSegm = std::max(vSegments[i].s_range.first,cfish->iFishStart);
			const int lastSegm = std::min(vSegments[i].s_range.second, cfish->iFishEnd);
			for(int ss=firstSegm; ss<=lastSegm; ++ss) {
				assert(ss>=cfish->iFishStart && ss<=cfish->iFishEnd);

				// fill chi
				// assume width is major axis, else correction:
				const Real offset = cfish->height[ss] > cfish->width[ss] ? 0.5*M_PI : 0.0;
				const Real ell_a = (Real)std::max(cfish->height[ss],cfish->width[ss]);
				//                    const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : 0.25*info.h[0]/ell_a;
				// maximum distance between two points is 2*ell_a * sin(theta). set this distance to dx/2 -->
				const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : std::abs(std::asin(0.5*(0.25*info.h_gridpoint)/ell_a));

				const int Ntheta = (int)std::ceil(2.0*M_PI/dtheta_target) + 1;
				const Real dtheta = 2.0*M_PI/((Real)Ntheta);

				for(int tt=0;tt<Ntheta;++tt) {
					const Real theta = tt*dtheta + offset;
					// create a surface point
					// special treatment of tail (width = 0 --> no ellipse, just line)
					const Real hght = cfish->width[ss] == 0 ? cfish->height[ss]*(2*tt/((Real)Ntheta-1) - 1)
							: cfish->height[ss]*std::sin(theta);
					Real myP[3] = {
							(cfish->rX[ss] + cfish->width[ss]*std::cos(theta)*cfish->norX[ss]),
							(cfish->rY[ss] + cfish->width[ss]*std::cos(theta)*cfish->norY[ss]),
							hght
					};
					changeToComputationalFrame(myP);

					const int iap[3] = {
							(int)std::floor((myP[0]-org[0])*invh),
							(int)std::floor((myP[1]-org[1])*invh),
							(int)std::floor((myP[2]-org[2])*invh)
					};
					// support is two points left, two points right --> Towers Chi will be one point left, one point right, but needs SDF wider
					const int start[3] = {
							std::max(-1, 0 - iap[0] ),
							std::max(-1, 0 - iap[1] ),
							std::max(-1, 0 - iap[2] )
					};
					const int end[3] = {
							std::min(+3, FluidBlock::sizeX - iap[0]),
							std::min(+3, FluidBlock::sizeY - iap[1]),
							std::min(+3, FluidBlock::sizeZ - iap[2])
					};
					const Real myP_distToMidlineSq = cfish->width[ss] == 0 ? std::pow(hght,2) :
							(Real)(std::pow(cfish->width[ss]*std::cos(theta),2) + std::pow(cfish->height[ss]*std::sin(theta),2));

					if(myP_distToMidlineSq<std::numeric_limits<Real>::epsilon()) {
						// if ss==iFishStart or ss==iFishEnd, the surface point lies on the midline --> myP_distToMidlineSq=0.
						// here our in/out criterion fails
						for(int sz=start[2]; sz<end[2];++sz)
						for(int sy=start[1]; sy<end[1];++sy)
						for(int sx=start[0]; sx<end[0];++sx) {
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);

							Real p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const Real diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const Real distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const Real distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);

							changeFromComputationalFrame(p);
							const Real distPlanar = std::sqrt( std::pow(p[0]-cfish->rX[closest_s],2) +
									std::pow(p[1]-cfish->rY[closest_s],2) );
							const Real distHeight = std::abs(p[2]);
							const Real sign = (distPlanar > cfish->width[closest_s] or distHeight > cfish->height[closest_s]) ? -1.0 : 1.0;

							defblock->chi[idx[2]][idx[1]][idx[0]] =
									(std::abs(defblock->chi[idx[2]][idx[1]][idx[0]]) > distSq) ? sign*distSq
											: defblock->chi[idx[2]][idx[1]][idx[0]];
						}
					} else {
						for(int sz=start[2]; sz<end[2];++sz)
						for(int sy=start[1]; sy<end[1];++sy)
						for(int sx=start[0]; sx<end[0];++sx) {
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);

							Real p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const Real diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const Real distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const Real distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);
							const Real sign = distToMidlineSq >= myP_distToMidlineSq ? -1.0 : 1.0;
							defblock->chi[idx[2]][idx[1]][idx[0]] =
									(std::abs(defblock->chi[idx[2]][idx[1]][idx[0]]) > distSq) ? sign*distSq
											: defblock->chi[idx[2]][idx[1]][idx[0]];
						}
					}
				}
			}
		}

		// construct the deformation velocities (P2M with hat function as kernel)

		for(int i=0;i<vSegments.size();++i) {
			for(int ss=vSegments[i].s_range.first;ss<=vSegments[i].s_range.second;++ss) {
				assert(ss>=0 && ss<=cfish->Nm-1);
				// P2M udef of a slice at this s
				const Real myWidth =  (ss < cfish->iFishStart ? cfish->width[ cfish->iFishStart]
																			  : (ss > cfish->iFishEnd   ? cfish->width[ cfish->iFishEnd] : cfish->width[ss]));
				const Real myHeight = (ss < cfish->iFishStart ? cfish->height[cfish->iFishStart]
																			  : (ss > cfish->iFishEnd   ? cfish->height[cfish->iFishEnd] : cfish->height[ss]));
				const Real ds_defGrid = info.h_gridpoint;
				// towers needs 1dx on each side, smooth needs 2dx --> make it 3 to be nice (and so we can floor!)
				const Real extension = NPPEXT*info.h_gridpoint; //G tmp changed back to 2

				const int Nh = std::floor( (myHeight+extension)/ds_defGrid );

				for(int ih=-Nh;ih<=Nh; ++ih) {
					const Real offsetH = ih*ds_defGrid;
					// add an extra extension when width == 0 (to deal with large curvatures near head and/or tail):
					const Real currentWidth = myWidth== 0 ? extension : myWidth * std::sqrt(1 - std::pow(offsetH/(myHeight+extension),2));
					const Real actualWidth = (cfish->height[ss] == 0 or std::abs(offsetH)>=cfish->height[ss]) ? 0.0
							: cfish->width[ss] * std::sqrt(1 - std::pow(offsetH/cfish->height[ss],2));
					const int Nw = std::floor( (currentWidth+extension)/ds_defGrid); // add extension here to make sure we have it in each direction

					for(int iw=-Nw;iw<=Nw; ++iw) {
						const Real offsetW = iw*ds_defGrid;
						Real xp[3] = {
								(cfish->rX[ss] + offsetW*cfish->norX[ss]),
								(cfish->rY[ss] + offsetW*cfish->norY[ss]),
								offsetH
						};
						changeToComputationalFrame(xp);
						xp[0] = (xp[0]-org[0])*invh;
						xp[1] = (xp[1]-org[1])*invh;
						xp[2] = (xp[2]-org[2])*invh;
						Real udef[3] = {
								(cfish->vX[ss] + offsetW*cfish->vNorX[ss]),
								(cfish->vY[ss] + offsetW*cfish->vNorY[ss]),
								0.0
						};
						changeVelocityToComputationalFrame(udef);
						const Real ap[3] = {
								std::floor((Real)xp[0]),
								std::floor((Real)xp[1]),
								std::floor((Real)xp[2])
						};
						const int iap[3] = {
								(int)ap[0],
								(int)ap[1],
								(int)ap[2]
						};
						// now we P2M
						const int start[3] = {
								std::max(0, 0 - iap[0] ),
								std::max(0, 0 - iap[1] ),
								std::max(0, 0 - iap[2] )
						};
						const int end[3] = {
								std::min(+2, FluidBlock::sizeX - iap[0]),
								std::min(+2, FluidBlock::sizeY - iap[1]),
								std::min(+2, FluidBlock::sizeZ - iap[2])
						};

						Real wghts[3][2];
						for(int c=0;c<3;++c) {
							const Real t[2] = {
									std::abs((Real)xp[c] - (ap[c]+0)),
									std::abs((Real)xp[c] - (ap[c]+1))
							};
							wghts[c][0] = 1.0 - t[0];
							wghts[c][1] = 1.0 - t[1];
						}

						const bool isInside = (std::abs(offsetW) < actualWidth) && (std::abs(offsetH) < cfish->height[ss]);
						for(int sz=start[2]; sz<end[2];++sz) {
						const Real wz = wghts[2][sz];
						for(int sy=start[1];sy<end[1];++sy) {
						const Real wywz = wz*wghts[1][sy];
						for(int sx=start[0];sx<end[0];++sx) {
							const Real wxwywz = wywz*wghts[0][sx];
							assert(wxwywz>=0 && wxwywz<=1);
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);
							defblock->udef[idx[2]][idx[1]][idx[0]][0] += wxwywz*udef[0];
							defblock->udef[idx[2]][idx[1]][idx[0]][1] += wxwywz*udef[1];
							defblock->udef[idx[2]][idx[1]][idx[0]][2] += wxwywz*udef[2];
							b(idx[0],idx[1],idx[2]).tmpU += wxwywz;
							// set sign for all interior points
							if( (std::abs(defblock->chi[idx[2]][idx[1]][idx[0]] + 1) < 5*std::numeric_limits<Real>::epsilon()) && isInside)
								defblock->chi[idx[2]][idx[1]][idx[0]] = 1.0;
						}
							}
						}
					}
				}
			}
		}

		// finalize signed distance function in tmpU
		{
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++) {
				const Real normfac = b(ix,iy,iz).tmpU > numeric_limits<Real>::epsilon() ? b(ix,iy,iz).tmpU : 1.;
				defblock->udef[iz][iy][ix][0] /= normfac;
				defblock->udef[iz][iy][ix][1] /= normfac;
				defblock->udef[iz][iy][ix][2] /= normfac;
				// change from signed squared distance function to normal sdf
				b(ix,iy,iz).tmpU = defblock->chi[iz][iy][ix] > 0 ?
						sqrt( defblock->chi[iz][iy][ix]) :
						-sqrt(-defblock->chi[iz][iy][ix]);
				b(ix,iy,iz).tmpV = defblock->udef[iz][iy][ix][0];
				b(ix,iy,iz).tmpW = defblock->udef[iz][iy][ix][1];
			}
		}
	}
};

struct PutFishOnBlocks_Finalize : public GenericLabOperator
{
	Real t;
	int stencil_start[3], stencil_end[3];
	array<Real,4>* const momenta;
	surfaceBlocks* const surface;
	std::map<int,ObstacleBlock*>* const obstacleBlocks;


	//PutFishOnBlocks_Finalize finalize(obstacleBlocks,dataPerThread[tid],tmp,blockID);

	PutFishOnBlocks_Finalize(map<int,ObstacleBlock*>* const obstacleBlocks, //to write chi
			surfaceBlocks* const surface, 					//to write gradChi
			array<Real,4>* const momenta)     			//compute CM
	: t(0), momenta(momenta),surface(surface), obstacleBlocks(obstacleBlocks)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, true, 1, 5);
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
	}
#if 1
	inline Real sign(const Real& val) const {
		return (0. < val) - (val < 0.);
	}
	template <typename Lab, typename BlockType>
	void operator()(Lab& lab, const BlockInfo& info, BlockType& b)
	{
		if(obstacleBlocks->find(info.blockID) == obstacleBlocks->end()) return;
		ObstacleBlock* const defblock = obstacleBlocks->find(info.blockID)->second;

		const Real eps  = std::numeric_limits<Real>::epsilon();
		const Real h    = info.h_gridpoint;
		const Real fac1 = .5/h;
		const Real fac2 = .5/(std::sqrt(2.)*h);
		const Real fac[9] = {fac1,fac1,fac1,fac2,fac2,fac2,fac2,fac2,fac2};

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU > +2*h || lab(ix,iy,iz).tmpU < -2*h) {
				const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
				(*momenta)[0] += H;
				(*momenta)[1] += p[0]*H;
				(*momenta)[2] += p[1]*H;
				(*momenta)[3] += p[2]*H;
				b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = H;
				continue;
			}
			const Real Up[9] = {
					lab(ix+1,iy  ,iz  ).tmpU, lab(ix  ,iy+1,iz  ).tmpU, lab(ix  ,iy  ,iz+1).tmpU,
					lab(ix  ,iy+1,iz+1).tmpU, lab(ix+1,iy  ,iz+1).tmpU, lab(ix+1,iy+1,iz  ).tmpU,
					lab(ix  ,iy+1,iz-1).tmpU, lab(ix+1,iy  ,iz-1).tmpU, lab(ix+1,iy-1,iz  ).tmpU
			};
			const Real Um[9] = {
					lab(ix-1,iy  ,iz  ).tmpU, lab(ix  ,iy-1,iz  ).tmpU, lab(ix  ,iy  ,iz-1).tmpU,
					lab(ix  ,iy-1,iz-1).tmpU, lab(ix-1,iy  ,iz-1).tmpU, lab(ix-1,iy-1,iz  ).tmpU,
					lab(ix  ,iy-1,iz+1).tmpU, lab(ix-1,iy  ,iz+1).tmpU, lab(ix-1,iy+1,iz  ).tmpU
			};
			Real Ip[9],Im[9];
			for (int i=0; i<9; i++) {
				Ip[i] = std::max(Up[i],(Real)0.);
				Im[i] = std::max(Um[i],(Real)0.);
			}
			Real gradU[9], gradI[9];
			for (int i=0; i<9; i++) {
				gradU[i] = fac[i]*(Up[i]-Um[i]);
				gradI[i] = fac[i]*(Ip[i]-Im[i]);
			}
			Real Hp[3], Hm[3]; //only x y and z, the direction where i need gradChi
			for(int i=0; i<3; i++) {
				Hp[i] = (Up[i]> h) ? h : (
						(Up[i]<-h) ? 0 :
								.5*h+(Up[i]-fac1*sign(Up[i])*Up[i]*Up[i]));
				Hm[i] = (Um[i]> h) ? h : (
						(Um[i]<-h) ? 0 :
								.5*h+(Um[i]-fac1*sign(Um[i])*Um[i]*Um[i]));
			}
			const Real gradH[3] = {.5*(Hp[0]-Hm[0]), .5*(Hp[1]-Hm[1]), .5*(Hp[2]-Hm[2])};
			Real gradUU[3], gradUI[3], gradUH[3];
			for (int i=0; i<3; i++) {
				gradUU[i] = gradU[i]*gradU[i] + gradU[i+3]*gradU[i+3] + gradU[i+6]*gradU[i+6];
				gradUI[i] = gradU[i]*gradI[i] + gradU[i+3]*gradI[i+3] + gradU[i+6]*gradI[i+6];
				gradUH[i] = gradU[i]*gradH[i];
			}
			for (int i=0; i<3; i++)  gradUU[i] = max(gradUU[i], eps);
			const Real FDD = h*(gradUH[0] + gradUH[1] + gradUH[2])/gradUU[0]; // == delta * h^3
			const Real FDH = 1/3. * (gradUI[0]/gradUU[0]+gradUI[1]/gradUU[1]+gradUI[2]/gradUU[2]);


			if (FDD>eps) {
				const Real dchidx = -FDD*gradU[0];
				const Real dchidy = -FDD*gradU[1];
				const Real dchidz = -FDD*gradU[2];
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, FDD);
			}
#ifndef NDEBUG
			//		if(FDH<0 || FDH>1) printf("invalid H?: %9.9e %9.9e %9.9e: %9.9e\n",x,y,z,FDH);
#endif
			(*momenta)[0] += FDH;
			(*momenta)[1] += p[0]*FDH;
			(*momenta)[2] += p[1]*FDH;
			(*momenta)[3] += p[2]*FDH;
			defblock->chi[iz][iy][ix] = FDH;
			b(ix,iy,iz).chi = std::max(FDH, b(ix,iy,iz).chi);
		}
	}
#else
	template <typename Lab, typename BlockType>
	void operator()(Lab& lab, const BlockInfo& info, BlockType& b)
	{
		if(obstacleBlocks->find(info.blockID) == obstacleBlocks->end()) return;
		ObstacleBlock* const defblock = obstacleBlocks->find(info.blockID)->second;
		const Real eps = std::numeric_limits<Real>::epsilon();
		const Real h = info.h_gridpoint;
		const Real inv2h = .5/h;

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU >= +2*h || lab(ix,iy,iz).tmpU <= -2*h) {
				const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
				b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = H;
				(*momenta)[1] += p[0]*H;
				(*momenta)[2] += p[1]*H;
				(*momenta)[3] += p[2]*H;
				(*momenta)[0] += H;
				continue;
			}

			const Real distPx = lab(ix+1,iy,iz).tmpU;
			const Real distMx = lab(ix-1,iy,iz).tmpU;
			const Real distPy = lab(ix,iy+1,iz).tmpU;
			const Real distMy = lab(ix,iy-1,iz).tmpU;
			const Real distPz = lab(ix,iy,iz+1).tmpU;
			const Real distMz = lab(ix,iy,iz-1).tmpU;

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

			// gradU
			const Real gradUX = inv2h * (distPx - distMx);
			const Real gradUY = inv2h * (distPy - distMy);
			const Real gradUZ = inv2h * (distPz - distMz);
			const Real gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ;

			// gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
			const Real gradIX = inv2h * (IplusX - IminuX);
			const Real gradIY = inv2h * (IplusY - IminuY);
			const Real gradIZ = inv2h * (IplusZ - IminuZ);
			const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;

			const Real gradHX = inv2h * (HplusX - HminuX);
			const Real gradHY = inv2h * (HplusY - HminuY);
			const Real gradHZ = inv2h * (HplusZ - HminuZ);
			const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;

			const Real Delta = std::abs(gradUSq) < eps ? numD : numD/gradUSq;
			const Real H     = std::abs(gradUSq) < eps ? numH : numH/gradUSq;

			if (Delta>1e-6) {
				const Real dchidx = -Delta*gradUX;
				const Real dchidy = -Delta*gradUY;
				const Real dchidz = -Delta*gradUZ;
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, Delta);
			}
			(*momenta)[0] += H;
			(*momenta)[1] += p[0]*H;
			(*momenta)[2] += p[1]*H;
			(*momenta)[3] += p[2]*H;
			defblock->chi[iz][iy][ix] = H;
			b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
		}
	}
#endif
};
}

#endif
