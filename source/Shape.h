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

#include "GeometryReader.h"
#include "WavefrontReader.h"
#include "CTReader.h"

class Shape
{
protected:
	// general quantities
	Geometry::Properties properties;
	
	// scale,translate
	Real scaleFactor;
	Geometry::Point translationFactor;
	
	// smoothing
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
	
public:
	Shape(Real center[3], const Real rhoS, const Real mollChi, const Real mollRho, Real scale=1, Real tX=0, Real tY=0, Real tZ=0, Geometry::Quaternion orientation=Geometry::Quaternion()) : mollChi(mollChi), mollRho(mollRho), scaleFactor(scale), translationFactor(tX,tY,tZ)
	{
		properties.com.x = center[0];
		properties.com.y = center[1];
		properties.com.z = center[2];
		
		properties.centroid.x = center[0];
		properties.centroid.y = center[1];
		properties.centroid.z = center[2];
		
		properties.ut.x = 0;
		properties.ut.y = 0;
		properties.ut.z = 0;
		
		properties.dthetadt.x = 0;
		properties.dthetadt.y = 0;
		properties.dthetadt.z = 0;
		
		properties.density = rhoS;
		
		properties.q = orientation;
		
		properties.rotation[0][0] = 1-2*(properties.q.y*properties.q.y + properties.q.z*properties.q.z);
		properties.rotation[0][1] =   2*(properties.q.x*properties.q.y - properties.q.w*properties.q.z);
		properties.rotation[0][2] =   2*(properties.q.x*properties.q.z + properties.q.w*properties.q.y);
		properties.rotation[1][0] =   2*(properties.q.x*properties.q.y + properties.q.w*properties.q.z);
		properties.rotation[1][1] = 1-2*(properties.q.x*properties.q.x + properties.q.z*properties.q.z);
		properties.rotation[1][2] =   2*(properties.q.y*properties.q.z - properties.q.w*properties.q.x);
		properties.rotation[2][0] =   2*(properties.q.x*properties.q.z - properties.q.w*properties.q.y);
		properties.rotation[2][1] =   2*(properties.q.y*properties.q.z + properties.q.w*properties.q.x);
		properties.rotation[2][2] = 1-2*(properties.q.x*properties.q.x + properties.q.y*properties.q.y);
		
		//cout << properties.q.w << " " << properties.q.x << " " << properties.q.y << " " << properties.q.z << endl;
		//cout << properties.rotation[0][0] << " " << properties.rotation[0][1] << " " << properties.rotation[0][2] << endl;
		//cout << properties.rotation[1][0] << " " << properties.rotation[1][1] << " " << properties.rotation[1][2] << endl;
		//cout << properties.rotation[2][0] << " " << properties.rotation[2][1] << " " << properties.rotation[2][2] << endl;
		
		//cout << "WARNING - shape not completely initialized\n";
		// to initialize: mass, minb, maxb, J
	}
	
	Shape()
	{
	}
	
	virtual ~Shape() {}
	
	virtual Real chi(Real p[3], Real h) const = 0;
	virtual Real getCharLength() const = 0;
	
	virtual void updatePosition(Real u[3], Real dthetadt[3], double J[6], Real mass, Real dt)
	{
		properties.ut.x = u[0];
		properties.ut.y = u[1];
		properties.ut.z = u[2];
		
#ifdef _J0_
		int idxI[3][3] = { 0,3,4, 3,1,5, 4,5,2 };
		
		Real J0[3][3], Jtmp[3][3], Jin[3][3];
		for (int k=0; k<3; k++)
			for (int i=0; i<3; i++)
			{
				J0[i][k] = 0;
				Jtmp[i][k] = 0;
				Jin[i][k] = properties.J0[idxI[i][k]];
			}
		
		for (int k=0; k<3; k++)
			for (int j=0; j<3; j++)
				for (int i=0; i<3; i++)
					Jtmp[i][k] += Jin[i][j] * properties.rotation[k][j];
		
		for (int k=0; k<3; k++)
			for (int j=0; j<3; j++)
				for (int i=0; i<3; i++)
					J0[i][k] += properties.rotation[i][j] * Jtmp[j][k];
		
		properties.J[0] = J0[0][0];
		properties.J[1] = J0[1][1];
		properties.J[2] = J0[2][2];
		properties.J[3] = J0[0][1];
		properties.J[4] = J0[0][2];
		properties.J[5] = J0[1][2];
		
		//cout << "J comp\t" << properties.J[0] << " " << properties.J[1] << " " << properties.J[2] << " " << properties.J[3] << " " << properties.J[4] << " " << properties.J[5] << endl;
		//cout << "J input\t" << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << " " << J[4] << " " << J[5] << endl;
		
#else
		for (int i=0; i<6; i++)
			properties.J[i] = J[i];
#endif
		
		//cout << "J\t" << properties.J[0] << " " << properties.J[1] << " " << properties.J[2] << " " << properties.J[3] << " " << properties.J[4] << " " << properties.J[5] << endl;
		const double detJ = properties.J[0]*(properties.J[1]*properties.J[2] - properties.J[5]*properties.J[5]) +
							properties.J[3]*(properties.J[4]*properties.J[5] - properties.J[2]*properties.J[3]) +
							properties.J[4]*(properties.J[3]*properties.J[5] - properties.J[1]*properties.J[4]);
		assert(abs(detJ)>numeric_limits<double>::epsilon());

		const double invDetJ = 1./detJ;
		
		const double invJ[6] = {
			invDetJ * (properties.J[1]*properties.J[2] - properties.J[5]*properties.J[5]),
			invDetJ * (properties.J[0]*properties.J[2] - properties.J[4]*properties.J[4]),
			invDetJ * (properties.J[0]*properties.J[1] - properties.J[3]*properties.J[3]),
			invDetJ * (properties.J[4]*properties.J[5] - properties.J[2]*properties.J[3]),
			invDetJ * (properties.J[3]*properties.J[5] - properties.J[1]*properties.J[4]),
			invDetJ * (properties.J[3]*properties.J[4] - properties.J[0]*properties.J[5])
		};
		
		//cout << "J\t" << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << " " << J[4] << " " << J[5] << endl;
		//cout << "J0\t" << properties.J[0] << " " << properties.J[1] << " " << properties.J[2] << " " << properties.J[3] << " " << properties.J[4] << " " << properties.J[5] << endl;
		//cout << "invJ\t" << invJ[0] << " " << invJ[1] << " " << invJ[2] << " " << invJ[3] << " " << invJ[4] << " " << invJ[5] << endl;
		
		// J-1 * dthetadt
		// angular velocity from angular momentum
		properties.dthetadt.x = (invJ[0]*dthetadt[0] + invJ[3]*dthetadt[1] + invJ[4]*dthetadt[2]);
		properties.dthetadt.y = (invJ[3]*dthetadt[0] + invJ[1]*dthetadt[1] + invJ[5]*dthetadt[2]);
		properties.dthetadt.z = (invJ[4]*dthetadt[0] + invJ[5]*dthetadt[1] + invJ[2]*dthetadt[2]);
		
		//cout << dthetadt[0] << " " << dthetadt[1] << " " << dthetadt[2] << endl;
		//cout << properties.dthetadt.x << " " << properties.dthetadt.y << " " << properties.dthetadt.z << endl;
		
		properties.update(dt);
	}
	
	void getAngularVelocity(Real omega[3])
	{
		omega[0] = properties.dthetadt.x;
		omega[1] = properties.dthetadt.y;
		omega[2] = properties.dthetadt.z;
	}
	
	void getCentroid(Real centroid[3])
	{
		centroid[0] = properties.centroid.x;
		centroid[1] = properties.centroid.y;
		centroid[2] = properties.centroid.z;
	}
	
	void getCenterOfMass(Real com[3])
	{
		com[0] = properties.com.x;
		com[1] = properties.com.y;
		com[2] = properties.com.z;
	}
	
	void setCentroid(Real centroid[3])
	{
		properties.centroid.x = centroid[0];
		properties.centroid.y = centroid[1];
		properties.centroid.z = centroid[2];
	}
	
	void setCenterOfMass(Real com[3])
	{
		properties.com.x = com[0];
		properties.com.y = com[1];
		properties.com.z = com[2];
	}
	
	double getMass()
	{
		return properties.mass;
	}
	
	void setMass(double mass)
	{
		properties.mass = mass;
	}
	
	void getInertiaMatrix(Real J[6])
	{
		for (int i=0; i<6; i++)
			J[i] = properties.J[i];
	}
	
	void setInertiaMatrix(Real J[6])
	{
		for (int i=0; i<6; i++)
			properties.J[i] = J[i];
	}
	
	void setInertiaMatrix0(double J[6])
	{
		for (int i=0; i<6; i++)
			properties.J0[i] = 0;//J[i];//
		//*
		// put into system of reference of unrotated shape - using transposed matrix
		Real J0[3][3], Jtmp[3][3], Jin[3][3];
		//Real Jinvtmp[3][3], Jinv[3][3];
		int idxI[3][3] = { 0,3,4, 3,1,5, 4,5,2 };
		for (int k=0; k<3; k++)
			for (int i=0; i<3; i++)
			{
				J0[i][k] = 0;
				Jtmp[i][k] = 0;
				//Jinvtmp[i][k] = 0;
				//Jinv[i][k] = 0;
				Jin[i][k] = J[idxI[i][k]];
			}
		
		for (int k=0; k<3; k++)
			for (int j=0; j<3; j++)
				for (int i=0; i<3; i++)
					Jtmp[i][k] += Jin[i][j] * properties.rotation[j][k];
		
		for (int k=0; k<3; k++)
			for (int j=0; j<3; j++)
				for (int i=0; i<3; i++)
					J0[i][k] += properties.rotation[j][i] * Jtmp[j][k];
		
		/*
		for (int k=0; k<3; k++)
			for (int j=0; j<3; j++)
				for (int i=0; i<3; i++)
					Jinvtmp[i][k] += J0[i][j] * properties.rotation[k][j];
		
		for (int k=0; k<3; k++)
			for (int j=0; j<3; j++)
				for (int i=0; i<3; i++)
					Jinv[i][k] += properties.rotation[i][j] * Jinvtmp[j][k];
		*/
		properties.J0[0] = J0[0][0];
		properties.J0[1] = J0[1][1];
		properties.J0[2] = J0[2][2];
		properties.J0[3] = J0[0][1];
		properties.J0[4] = J0[0][2];
		properties.J0[5] = J0[1][2];
				
		
		//cout << "J0 input\t" << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << " " << J[4] << " " << J[5] << endl;
		//cout << "J0 matrix\t" << Jin[0][0] << " " << Jin[0][1] << " " << Jin[0][2] << " " << Jin[1][0] << " " << Jin[1][1] << " " << Jin[1][2] << " " << Jin[2][0] << " " << Jin[2][1] << " " << Jin[2][2] << endl;
		//cout << "J0ref mat\t" << J0[0][0] << " " << J0[0][1] << " " << J0[0][2] << " " << J0[1][0] << " " << J0[1][1] << " " << J0[1][2] << " " << J0[2][0] << " " << J0[2][1] << " " << J0[2][2] << endl;
		//cout << "J0ref mat\t" << J0[0][0] << " " << J0[1][1] << " " << J0[2][2] << " " << J0[0][1] << " " << J0[0][2] << " " << J0[1][2] << endl;
		//cout << "Jinv mat\t" << Jinv[0][0] << " " << Jinv[0][1] << " " << Jinv[0][2] << " " << Jinv[1][0] << " " << Jinv[1][1] << " " << Jinv[1][2] << " " << Jinv[2][0] << " " << Jinv[2][1] << " " << Jinv[2][2] << endl;
		//cout << "Jinv mat\t" << Jinv[0][0] << " " << Jinv[1][1] << " " << Jinv[2][2] << " " << Jinv[0][1] << " " << Jinv[0][2] << " " << Jinv[1][2] << endl;
	}
	
	void getOrientation(Real rotation[3][3]) const
	{
		rotation[0][0] = properties.rotation[0][0];
		rotation[0][1] = properties.rotation[0][1];
		rotation[0][2] = properties.rotation[0][2];
		rotation[1][0] = properties.rotation[1][0];
		rotation[1][1] = properties.rotation[1][1];
		rotation[1][2] = properties.rotation[1][2];
		rotation[2][0] = properties.rotation[2][0];
		rotation[2][1] = properties.rotation[2][1];
		rotation[2][2] = properties.rotation[2][2];
	}
	
	void getOrientation(Real q0, Real q1, Real q2, Real q3) const
	{
		q0 = properties.q.w;
		q1 = properties.q.x;
		q2 = properties.q.y;
		q3 = properties.q.z;
	}
	
	virtual inline Real getRhoS() const
	{
		return properties.density;
	}
	
	virtual Real rho(Real p[3], Real h, Real mask) const
	{
		return properties.density*mask + 1.*(1.-mask);
	}
	
	virtual Real rho(Real p[3], Real h) const
	{
		Real mask = chi(p,h);
		return rho(p,h,mask);
	}
	
	virtual void serialize(ostream &outStream)
	{
		outStream << "scaleFactor " << scaleFactor << endl;
		outStream << "translationFactor " << translationFactor.x << " " << translationFactor.y << " " << translationFactor.z << endl;
		outStream << "mollChi " << mollChi << endl;
		outStream << "mollRho " << mollRho << endl;
		
		properties.serialize(outStream);
	}
	
	void setProperties(Geometry::Properties p)
	{
		properties = p;
	}
	
	// no need for deserialization - this task is delegated to static function of descendants
};

class Sphere : public Shape
{
protected:
	Real radius;
	
public:
	Sphere(Real center[3], Real radius, const Real rhoS, const Real mollChi, const Real mollRho) : Shape(center, rhoS, mollChi, mollRho), radius(radius)
	{
	}
	
	Real chi(Real p[3], Real h) const
	{
		const Real d[3] = { abs(p[0]-properties.centroid.x), abs(p[1]-properties.centroid.y), abs(p[2]-properties.centroid.z) };
		const Real dist = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		
		return smoothHeaviside(dist, radius, mollChi*sqrt(2)*h);
	}
	
	Real getCharLength() const
	{
		return 2 * radius;
	}
	
	/*
	void outputSettings(ostream &outStream)
	{
		outStream << "Sphere\n";
		outStream << "radius " << radius << endl;
		
		Shape::serialize(outStream);
	}
	*/
};

class Ellipsoid : public Shape
{
protected:
	Real semiAxis[3];
	
public:
	Ellipsoid(Real center[3], Real semiAxis[3], const Real rhoS, const Real mollChi, const Real mollRho) : Shape(center, rhoS, mollChi, mollRho), semiAxis{semiAxis[0],semiAxis[1],semiAxis[2]}
	{
	}
	
	Real getCharLength() const
	{
		return 2 * semiAxis[1];
	}
};

class GeometryMesh : public Shape
{
protected:
	GeometryReader * geometry;
	const Real isosurface;
	const Real charSize;
	
public:
	GeometryMesh(const string filename, const int gridsize, const Real isosurface, Real center[3], const Real charSize, const Real rhoS, const Real mollChi, const Real mollRho, Real scale=1, Real tX=0, Real tY=0, Real tZ=0, Geometry::Quaternion orientation=Geometry::Quaternion()) : Shape(center, rhoS, mollChi, mollRho, scale, tX, tY, tZ, orientation), isosurface(isosurface), charSize(charSize)
	{
		Geometry::Point transFactor(tX,tY,tZ);
		geometry = new GeometryReaderOBJ(filename,properties,gridsize,scale,transFactor,isosurface);
	}
	
	Real chi(Real p[3], Real h) const
	{
		return smoothHeaviside(geometry->distance(p[0], p[1], p[2]), isosurface, mollChi*sqrt(2)*h);
	}
	
	Real getCharLength() const
	{
		return charSize;
	}
	/*
	void outputSettings(ostream &outStream)
	{
		//cout << "Not implemented yet\n";
		//abort();
	}
	*/
	virtual void serialize(ostream &outStream)
	{
		// these are information necessary for the constructor on restart!
		outStream << "filename " << geometry->filename << endl;
		outStream << "gridsize " << geometry->gridsize << endl;
		outStream << "isosurface " << isosurface << endl;
		outStream << "charSize " << charSize << endl;
		
		Shape::serialize(outStream);
	}
	
	static GeometryMesh * deserialize(istream& inStream)
	{
		string variableName;
		string filename;
		Real scale;
		Real tX,tY,tZ;
		int gridsize;
		Real isosurface;
		Real mChi, mRho;
		Real charSize;
		
		inStream >> variableName;
		assert(variableName=="filename");
		inStream >> filename;
		inStream >> variableName;
		assert(variableName=="gridsize");
		inStream >> gridsize;
		inStream >> variableName;
		assert(variableName=="isosurface");
		inStream >> isosurface;
		inStream >> variableName;
		assert(variableName=="charSize");
		inStream >> charSize;
		
		inStream >> variableName;
		assert(variableName=="scaleFactor");
		inStream >> scale;
		inStream >> variableName;
		assert(variableName=="translationFactor");
		inStream >> tX;
		inStream >> tY;
		inStream >> tZ;
		inStream >> variableName;
		assert(variableName=="mollChi");
		inStream >> mChi;
		inStream >> variableName;
		assert(variableName=="mollRho");
		inStream >> mRho;
		
		Geometry::Properties properties;
		properties.deserialize(inStream);
		
		Real center[3] = { properties.com.x, properties.com.y, properties.com.z };
		
		GeometryMesh * mesh = new GeometryMesh(filename,gridsize,isosurface,center,charSize,properties.density,mChi,mRho,scale,tX,tY,tZ,properties.q);
		mesh->setProperties(properties);
		return mesh;
	}
};

#endif
