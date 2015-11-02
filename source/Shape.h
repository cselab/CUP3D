//
//  Shape.h
//  CubismUP_3D
//
//	Virtual shape class which defines the interface
//	Default simple geometries are also provided and can be used as references
//
//	This class only contains static information (position, orientation,...), no dynamics are included (e.g. velocities,...)
//
//  Created by Christian Conti on 3/6/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_Shape_h
#define CubismUP_3D_Shape_h

class Shape
{
protected:
	Real center[2], orientation;
	const Real rhoS;
    
    const Real domainSize[2];
    const bool bPeriodic[2];
	
	const Real mollChi;
	const Real mollRho; // currently not used - need to change in rho method
	
	Real smoothHeaviside(Real rR, Real radius, Real eps)
	{
		if (rR < radius-eps*.5)
			return (Real) 1.;
		else if (rR > radius+eps*.5)
			return (Real) 0.;
		else
			return (Real) ((1.+cos(M_PI*((rR-radius)/eps+.5)))*.5);
	}
	
	// TODO: need to include boundary conditions in case of periodicity
	
public:
    Shape(Real center[2], Real orientation, const Real rhoS, const Real mollChi, const Real mollRho, bool bPeriodic[2], Real domainSize[2]) : center{center[0],center[1]}, domainSize{domainSize[0],domainSize[1]}, bPeriodic{bPeriodic[0],bPeriodic[1]}, orientation(orientation), rhoS(rhoS), mollChi(mollChi), mollRho(mollRho)
	{
		if (bPeriodic[0] || bPeriodic[1])
		{
			cout << "Periodic shapes are currently unsupported\n";
			abort();
		}
	}
	
	virtual ~Shape() {}
	
	virtual Real chi(Real p[2], Real h) const = 0;
	virtual Real getCharLength() const = 0;
	
	
	void updatePosition(Real u[2], Real omega, Real dt)
	{
		center[0] += dt*u[0];
#ifndef _MOVING_FRAME_
        center[1] += dt*u[1];
#endif
		orientation += dt*omega;
		orientation = orientation>2*M_PI ? orientation-2*M_PI : (orientation<0 ? orientation+2*M_PI : orientation);
	}
	
	void setPosition(Real com[2])
	{
		center[0] = com[0];
		center[1] = com[1];
	}
	
	void getPosition(Real com[2])
	{
		com[0] = center[0];
		com[1] = center[1];
	}
	
	Real getOrientation() const
	{
		return orientation;
	}
	
	inline Real getRhoS() const
	{
		return rhoS;
	}
	
	Real rho(Real p[2], Real h)
	{
		Real mask = chi(p,h);
		
		return rhoS*mask + 1.*(1.-mask);
	}
	
	virtual void outputSettings(ostream &outStream)
	{
		outStream << "centerX " << center[0] << endl;
		outStream << "centerY " << center[1] << endl;
		outStream << "orientation " << orientation << endl;
		outStream << "rhoS " << rhoS << endl;
		outStream << "mollChi " << mollChi << endl;
		outStream << "mollRho " << mollRho << endl;
	}
};

class Disk : public Shape
{
protected:
	Real radius;
	
public:
	Disk(Real center[2], Real radius, const Real rhoS, const Real mollChi, const Real mollRho, bool bPeriodic[2], Real domainSize[2]) : Shape(center, 0, rhoS, mollChi, mollRho, bPeriodic, domainSize), radius(radius)
    {
    }
	
	Real chi(Real p[2], Real h) const
    {
        const Real centerPeriodic[2] = {center[0] - floor(center[0]/domainSize[0]) * bPeriodic[0],
                                        center[1] - floor(center[1]/domainSize[1]) * bPeriodic[1]};
        
        const Real d[2] = { abs(p[0]-centerPeriodic[0]), abs(p[1]-centerPeriodic[1]) };
		const Real dist = sqrt(d[0]*d[0]+d[1]*d[1]);
		
		return smoothHeaviside(dist, radius, mollChi*sqrt(2)*h);
	}
	
	Real getCharLength() const
	{
		return 2 * radius;
	}
	
	void outputSettings(ostream &outStream)
	{
		outStream << "Disk\n";
		outStream << "radius " << radius << endl;
		
		Shape::outputSettings(outStream);
	}
};

class Ellipse : public Shape
{
protected:
	// these quantities are defined in the local coordinates of the ellipse
	Real semiAxis[2];
	
	// code from http://www.geometrictools.com/
	//----------------------------------------------------------------------------
	// The ellipse is (x0/semiAxis0)^2 + (x1/semiAxis1)^2 = 1.  The query point is (y0,y1).
	// The function returns the distance from the query point to the ellipse.
	// It also computes the ellipse point (x0,x1) that is closest to (y0,y1).
	//----------------------------------------------------------------------------
	Real DistancePointEllipseSpecial (const Real e[2], const Real y[2], Real x[2])
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
	
	Real DistancePointEllipse(const Real y[2], Real x[2])
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
		if (semiAxis[0] < semiAxis[1])
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
			locE[i] = semiAxis[j];
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
	
public:
	Ellipse(Real center[2], Real semiAxis[2], Real orientation, const Real rhoS, const Real mollChi, const Real mollRho, bool bPeriodic[2], Real domainSize[2]) : Shape(center, orientation, rhoS, mollChi, mollRho, bPeriodic, domainSize), semiAxis{semiAxis[0],semiAxis[1]} {}
	
	Real chi(Real p[2], Real h) const
    {
		/*
		// based on https://www.spaceroots.org/documents/distance/distance-to-ellipse.pdf
		const int maxIter = 100;
		const Real threshold  = .001;
		const Real thresholdC = .001;
		Real theta, g;
		Real a,b,c,d,k;
		Real q,r,e;
		Real t, tTmp, phiTmp, phi;
		
		p[0] -= center[0];
		p[1] -= center[1];
		const Real dist = sqrt(p[0]*p[0]+p[1]*p[1]);
		
		const Real ae = max(semiAxis[0],semiAxis[1]);
		const Real ap = min(semiAxis[0],semiAxis[1]);
		const Real f = 1. - ap/ae;
		const Real omf2 = (1.-f)*(1.-f);
		
		if (dist < threshold * ae)
		{
			g = -ap;
			return g;//smoothHeaviside(g, 0, mollChi*sqrt(2)*h);
		}
		
		const Real radius = omf2 * (p[0]*p[0]-ae*ae) + p[1]*p[1];
		const bool bInside = radius <= 0;
		
		const Real cosz = p[0]/dist;
		const Real sinz = p[1]/dist;
		t = p[1]/(p[0]+dist);
		
		// distance to the ellipse along the current line as the root of a 2nd degree
		// polynom closest to zero: ak2 −2bk +c = 0
		a = omf2 * cosz*cosz + sinz*sinz;
		b = omf2 * p[0]*cosz + p[1]*sinz;
		c = radius;
		k = c / (b - sqrt(b*b-a*c));
		phi = atan2(p[1] - k*sinz, omf2*(p[0]-k*cosz));
		
		if (abs(k) < threshold*dist)
		{
			g = k;
			return g;//smoothHeaviside(g, 0, mollChi*sqrt(2)*h);
		}
		
		for (int i=0; i<maxIter; i++)
		{
			a = omf2*((p[0]+k)*(p[0]+k) - ae*ae) + p[1]*p[1];
			b = -4.*k*p[1]/a;
			c = 2. * (omf2*(p[0]*p[0] - k*k - ae*ae) + 2.*k*k + p[1]*p[1]) / a;
			d = b;
			
			b += t;
			c += t;
			d += t;
			
			q = (3.*c - b*b) / 9.;
			r = (b*(9.*c-2.*b*b) - 27.*d) / 54.;
			e = q*q*q + r*r;
			
			if (e >= 0)
			{
				tTmp = cbrt(r+sqrt(e)) + cbrt(r-sqrt(e)) - b/3.;
				phiTmp = atan2(p[1] * (1. + tTmp*tTmp) - 2.*k*tTmp, omf2 * (p[0] * (1.+tTmp*tTmp) - k*(1.-tTmp*tTmp)));
			}
			else
			{
				q = -q;
				theta = acos(r) / (q*sqrt(q));
				tTmp = 2.*sqrt(q)*cos(theta/3.) - b/3.;
				phiTmp = atan2(p[1] * (1. + tTmp*tTmp) - 2.*k*tTmp, omf2 * (p[0] * (1.+tTmp*tTmp) - k*(1.-tTmp*tTmp)));
				
				if (phi*phiTmp < 0)
				{
					tTmp = 2.*sqrt(q)*cos((theta+2.*M_PI)/3.) - b/3.;
					phiTmp = atan2(p[1] * (1. + tTmp*tTmp) - 2.*k*tTmp, omf2 * (p[0] * (1.+tTmp*tTmp) - k*(1.-tTmp*tTmp)));
					
					
					if (phi*phiTmp < 0)
					{
						tTmp = 2.*sqrt(q)*cos((theta+4.*M_PI)/3.) - b/3.;
						phiTmp = atan2(p[1] * (1. + tTmp*tTmp) - 2.*k*tTmp, omf2 * (p[0] * (1.+tTmp*tTmp) - k*(1.-tTmp*tTmp)));
						
					}
				}
			}
			
			const Real deltaPhi = abs(phiTmp-phi) / 2.;
			phi = (phi+phiTmp) / 2.;
			
			if (deltaPhi < thresholdC)
			{
				g = p[0]*cos(phi) + p[1]*sin(phi) - ae * sqrt(1.-f*(2.-f)*sin(phi)*sin(phi));
				return g;//smoothHeaviside(g, 0, mollChi*sqrt(2)*h);
			}
			
			const Real deltaX = p[0] - ae        * cos(phi) / sqrt(1.-f*(2.-f)*sin(phi)*sin(phi));
			const Real deltaY = p[1] - ae * omf2 * sin(phi) / sqrt(1.-f*(2.-f)*sin(phi)*sin(phi));
			k = sqrt(deltaX*deltaX + deltaY*deltaY);
			
			if (bInside) k = -k;
			
			t = deltaY/(deltaX+k);
		}
		*/
		//*
		const Real centerPeriodic[2] = {center[0] - floor(center[0]/domainSize[0]) * bPeriodic[0],
										center[1] - floor(center[1]/domainSize[1]) * bPeriodic[1]};
		Real x[2] = {0,0};
		p[0]-=centerPeriodic[0];
		p[1]-=centerPeriodic[1];
		
		const Real rotatedP[2] = { cos(orientation)*p[1] - sin(orientation)*p[0],
								   sin(orientation)*p[1] + cos(orientation)*p[0] };
		const Real dist = DistancePointEllipse(rotatedP, x);
		const int sign = ( (rotatedP[0]*rotatedP[0]+rotatedP[1]*rotatedP[1]) > (x[0]*x[0]+x[1]*x[1]) ) ? 1 : -1;
		
		return smoothHeaviside(sign*dist,0,mollChi*sqrt(2)*h);
		
		/*/
        const Real centerPeriodic[2] = {center[0] - floor(center[0]/domainSize[0]) * bPeriodic[0],
                                        center[1] - floor(center[1]/domainSize[1]) * bPeriodic[1]};
        
		const Real angle = atan2(p[1]-centerPeriodic[1],p[0]-centerPeriodic[0]) - orientation;
		const Real x = semiAxis[0]*cos(angle);
		const Real y = semiAxis[1]*sin(angle);
		const Real radius = semiAxis[0]*semiAxis[1] / sqrt(x*x + y*y);
		const Real dist = sqrt( (p[0]-centerPeriodic[0])*(p[0]-centerPeriodic[0]) + (p[1]-centerPeriodic[1])*(p[1]-centerPeriodic[1]) );
		return smoothHeaviside(dist, radius, mollChi*sqrt(2)*h);
		/*/
		/*
		 // need to account for rotations!
		const Real focus = sqrt(semiAxis[1]*semiAxis[1]-semiAxis[0]*semiAxis[0]);
		const Real distF1 = sqrt((p[0]-(center[0]-focus))*(p[0]-(center[0]-focus))+(p[1]-center[1])*(p[1]-center[1]));
		const Real distF2 = sqrt((p[0]-(center[0]+focus))*(p[0]-(center[0]+focus))+(p[1]-center[1])*(p[1]-center[1]));
		const Real dist = distF1+distF2 - (2*semiAxis[1]);
		return smoothHeaviside(dist, 0, mollChi*sqrt(2)*h);
		return smoothHeaviside(distF1+distF2, 2*semiAxis[1], mollChi*sqrt(2)*h);
		//*/
	}
	
	Real getCharLength() const
	{
		return 2 * semiAxis[1];
	}
	
	void outputSettings(ostream &outStream)
	{
		outStream << "Ellipse\n";
		outStream << "semiAxisX " << semiAxis[0] << endl;
		outStream << "semiAxisY " << semiAxis[1] << endl;
		
		Shape::outputSettings(outStream);
	}
};

#endif
