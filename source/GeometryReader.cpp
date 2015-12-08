//
//  GeometryReader.cpp
//  CubismUP_3D
//
//  Abstract class for geometry reader
//
//  Created by Christian Conti on 06/14/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//


#include "GeometryReader.h"


// constructor loading OBJ file
GeometryReader::GeometryReader(const string filename, Geometry::Properties & properties, int gridsize, const Real scaleFactor, const Geometry::Point transFactor) : bGeometryLoaded(false), gridsize(gridsize), bVerbose(true), filename(filename), baseScale(scaleFactor), baseTranslation(transFactor), properties(properties)
{
}

// clean up
GeometryReader::~GeometryReader()
{
}

void GeometryReader::serialize(string path) const
{
	stringstream serializedFile;
	serializedFile << path << "/geometry.dat";
	
	ofstream file(serializedFile.str().c_str());
	
	file << gridsize << " "
	<< properties.com.x << " "
	<< properties.com.y << " "
	<< properties.com.z << " "
	<< properties.ut.x << " "
	<< properties.ut.y << " "
	<< properties.ut.z << " "
	<< properties.dthetadt.x << " "
	<< properties.dthetadt.y << " "
	<< properties.dthetadt.z << " "
	<< properties.density << " "
	<< properties.mass << " "
	<< properties.minb.x << " "
	<< properties.minb.y << " "
	<< properties.minb.z << " "
	<< properties.maxb.x << " "
	<< properties.maxb.y << " "
	<< properties.maxb.z << " "
	<< properties.q.w << " "
	<< properties.q.x << " "
	<< properties.q.y << " "
	<< properties.q.z << " "
	<< properties.J[0] << " "
	<< properties.J[1] << " "
	<< properties.J[2] << " "
	<< properties.J[3] << " "
	<< properties.J[4] << " "
	<< properties.J[5] << " "
	<< properties.rotation[0][0] << " "
	<< properties.rotation[0][1] << " "
	<< properties.rotation[0][2] << " "
	<< properties.rotation[1][0] << " "
	<< properties.rotation[1][1] << " "
	<< properties.rotation[1][2] << " "
	<< properties.rotation[2][0] << " "
	<< properties.rotation[2][1] << " "
	<< properties.rotation[2][2];
}

void GeometryReader::deserialize(string path)
{
	stringstream serializedFile;
	serializedFile << path << "/geometry.dat";
	
	ifstream file(serializedFile.str().c_str());
	
	file >> gridsize;
	file.ignore(1,' ');
	file >> properties.com.x;
	file.ignore(1,' ');
	file >> properties.com.y;
	file.ignore(1,' ');
	file >> properties.com.z;
	file.ignore(1,' ');
	file >> properties.ut.x;
	file.ignore(1,' ');
	file >> properties.ut.y;
	file.ignore(1,' ');
	file >> properties.ut.z;
	file.ignore(1,' ');
	file >> properties.dthetadt.x;
	file.ignore(1,' ');
	file >> properties.dthetadt.y;
	file.ignore(1,' ');
	file >> properties.dthetadt.z;
	file.ignore(1,' ');
	file >> properties.density;
	file.ignore(1,' ');
	file >> properties.mass;
	file.ignore(1,' ');
	file >> properties.minb.x;
	file.ignore(1,' ');
	file >> properties.minb.y;
	file.ignore(1,' ');
	file >> properties.minb.z;
	file.ignore(1,' ');
	file >> properties.maxb.x;
	file.ignore(1,' ');
	file >> properties.maxb.y;
	file.ignore(1,' ');
	file >> properties.maxb.z;
	file.ignore(1,' ');
	file >> properties.q.w;
	file.ignore(1,' ');
	file >> properties.q.x;
	file.ignore(1,' ');
	file >> properties.q.y;
	file.ignore(1,' ');
	file >> properties.q.z;
	file.ignore(1,' ');
	for (int i=0; i<6; i++)
	{
		file >> properties.J[i];
		file.ignore(1,' ');
	}
	for (int j=0; j<3; j++)
		for (int i=0; i<3; i++)
		{
			file >> properties.rotation[j][i];
			file.ignore(1,' ');
		}
}

// equivalent to "finalize" from Wim
// gravity is not currently used
void GeometryReader::setBodyIntegrals(double ut[3], double dthetadt[3], double J[6], double mass, double a, double dt)
{
	//cout << "Mass difference " << abs(properties.mass - volume * properties.density) << "(" << properties.mass << " " << volume * properties.density << " " << properties.mass/(volume * properties.density) << ")" << endl;
	properties.mass = mass;
	
	properties.ut.x = properties.ut.x * a + dt * (ut[0]/properties.mass);
	properties.ut.y = properties.ut.y * a + dt * (ut[1]/properties.mass);
	properties.ut.z = properties.ut.z * a + dt * (ut[2]/properties.mass);
	
	for (int i=0; i<6; i++)
		properties.J[i] = J[i];
	//	cout << "Inertia difference for component " << i << "\t" << abs(properties.J[i] - J[i]) << "(" << properties.J[i] << " " << J[i] << " " << properties.J[i]/J[i] << ")" << endl;
	
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
	
	// J-1 * dthetadt
	// angular velocity from angular momentum
	properties.dthetadt.x = properties.dthetadt.x * a + dt * (invJ[0]*dthetadt[0] + invJ[3]*dthetadt[1] + invJ[4]*dthetadt[2]);
	properties.dthetadt.y = properties.dthetadt.y * a + dt * (invJ[3]*dthetadt[0] + invJ[1]*dthetadt[1] + invJ[5]*dthetadt[2]);
	properties.dthetadt.z = properties.dthetadt.z * a + dt * (invJ[4]*dthetadt[0] + invJ[5]*dthetadt[1] + invJ[2]*dthetadt[2]);
	
	//cout << properties.dthetadt.x << " " << properties.dthetadt.y << " " << properties.dthetadt.z << endl;
}

void GeometryReader::setConstIntegrals(double J[6], double volume)
{
	properties.mass = volume * properties.density;
	for (int i=0; i<6; i++)
		properties.J[i] = J[i] * properties.density;
}

// get bounds of mesh
void GeometryReader::getBounds(double& x1, double& x2, double& y1, double& y2, double& z1, double& z2) const
{
	if (!bGeometryLoaded)
	{
		cout << "getBounds - Geometry not loaded yet!\n";
		abort();
	}
	
	x1 = properties.minb.x;
	x2 = properties.maxb.x;
	y1 = properties.minb.y;
	y2 = properties.maxb.y;
	z1 = properties.minb.z;
	z2 = properties.maxb.z;
}

// get center of geometry - unused
void GeometryReader::getCentroid(double p[3]) const
{
	if (!bGeometryLoaded)
	{
		cout << "getCenter - Geometry not loaded yet!\n";
		abort();
	}
	
	p[0] = (properties.centroid.x-properties.com.x)*properties.rotation[0][0] + (properties.centroid.y-properties.com.y)*properties.rotation[1][0] + (properties.centroid.z-properties.com.z)*properties.rotation[2][0] + properties.com.x;
	p[1] = (properties.centroid.x-properties.com.x)*properties.rotation[0][1] + (properties.centroid.y-properties.com.y)*properties.rotation[1][1] + (properties.centroid.z-properties.com.z)*properties.rotation[2][1] + properties.com.y;
	p[2] = (properties.centroid.x-properties.com.x)*properties.rotation[0][2] + (properties.centroid.y-properties.com.y)*properties.rotation[1][2] + (properties.centroid.z-properties.com.z)*properties.rotation[2][2] + properties.com.z;
}

void GeometryReader::getCenterOfMass(double p[3]) const
{
	if (!bGeometryLoaded)
	{
		cout << "getCenterOfMass - Geometry not loaded yet!\n";
		abort();
	}
	
	p[0] = properties.com.x;
	p[1] = properties.com.y;
	p[2] = properties.com.z;
}

void GeometryReader::getUT(double ut[3]) const
{
	if (!bGeometryLoaded)
	{
		cout << "getUT - Geometry not loaded yet!\n";
		abort();
	}
	
	ut[0] = properties.ut.x;
	ut[1] = properties.ut.y;
	ut[2] = properties.ut.z;
}

void GeometryReader::getAngVel(double angVel[3]) const
{
	if (!bGeometryLoaded)
	{
		cout << "getAngVel - Geometry not loaded yet!\n";
		abort();
	}
	
	angVel[0] = properties.dthetadt.x;
	angVel[1] = properties.dthetadt.y;
	angVel[2] = properties.dthetadt.z;
}

void GeometryReader::getMomentOfInertia(double J[6]) const
{
	if (!bGeometryLoaded)
	{
		cout << "getMomentOfInertia - Geometry not loaded yet!\n";
		abort();
	}
	
	J[0] = properties.J[0];
	J[1] = properties.J[1];
	J[2] = properties.J[2];
	J[3] = properties.J[3];
	J[4] = properties.J[4];
	J[5] = properties.J[5];
}

double GeometryReader::getDensity() const
{
	return properties.density;
}

double GeometryReader::getMass() const
{
	if (!bGeometryLoaded)
	{
		cout << "getMass - Geometry not loaded yet!\n";
		abort();
	}
	
	return properties.mass;
}

// check if a point lies inside the body
int GeometryReader::isPointInside(int ix, int iy, int iz)
{
	if (!bGeometryLoaded)
	{
		cout << "isPointInside - Geometry not loaded yet!\n";
		abort();
	}
	
	// assuming the surface is closed and fully contained within the grid,
	// taking the nearest gridpoint and scanning in a given direction until a nonzero value is found,
	// can give the in/out information
	// additionally, the gridpoint can be set so that for other accesses, the search is faster
	
	float signx = 0;
	for (int px=ix; px<grid.sizeX; px++)
	{
		if (grid(px,iy,iz)<10000.)
		{
			signx = grid(px,iy,iz)>0.?1.:-1.;
			break;
		}
	}
	if (signx==0.) signx = 1.;
	grid(ix,iy,iz) = signx;
	return grid(ix,iy,iz);
}

// return signed distance to closest point on surface (neg. if inside the body).
double GeometryReader::distance(double fx, double fy, double fz)
{
	// TODO:
	//	eventually remove intermediate grid and P2M directly on the output
	
	// rotate point to find position in rotated frame - should be w.r.t. center of mass
	// careful with left vs right multiplication!
	// careful with correct rotation direction
	double com[3];
	getCenterOfMass(com);
#if 0
	// inverse angle - reverse rotation - incorrect
	double rfx = (fx-com[0])*properties.rotation[0][0] + (fy-com[1])*properties.rotation[1][0] + (fz-com[2])*properties.rotation[2][0] + com[0];
	double rfy = (fx-com[0])*properties.rotation[0][1] + (fy-com[1])*properties.rotation[1][1] + (fz-com[2])*properties.rotation[2][1] + com[1];
	double rfz = (fx-com[0])*properties.rotation[0][2] + (fy-com[1])*properties.rotation[1][2] + (fz-com[2])*properties.rotation[2][2] + com[2];
#else
	// right-hand rule
	double rfx = (fx-com[0])*properties.rotation[0][0] + (fy-com[1])*properties.rotation[0][1] + (fz-com[2])*properties.rotation[0][2] + com[0];
	double rfy = (fx-com[0])*properties.rotation[1][0] + (fy-com[1])*properties.rotation[1][1] + (fz-com[2])*properties.rotation[1][2] + com[1];
	double rfz = (fx-com[0])*properties.rotation[2][0] + (fy-com[1])*properties.rotation[2][1] + (fz-com[2])*properties.rotation[2][2] + com[2];
#endif
	
	if (!bSDFComputed)
	{
		cout << "distance - SDF not computed yet!\n";
		abort();
	}
	
	if (rfx < properties.maxb.x && rfx > properties.minb.x &&
		rfy < properties.maxb.y && rfy > properties.minb.y &&
		rfz < properties.maxb.z && rfz > properties.minb.z)
	{
		// nearest gridpoint
		int ix = (int)((rfx-properties.minb.x)/(double)grid.maxsize*grid.maxpoints);
		int iy = (int)((rfy-properties.minb.y)/(double)grid.maxsize*grid.maxpoints);
		int iz = (int)((rfz-properties.minb.z)/(double)grid.maxsize*grid.maxpoints);
		
#if 1
		// trilinear interpolation - should be ok since the signed distance function is a linear function
		int dx = (rfx-properties.minb.x)/(double)grid.maxsize*grid.maxpoints - ix;
		int dy = (rfy-properties.minb.y)/(double)grid.maxsize*grid.maxpoints - iy;
		int dz = (rfz-properties.minb.z)/(double)grid.maxsize*grid.maxpoints - iz;
		
		double v000 = grid.read(ix  ,iy  ,iz  );
		double v001 = grid.read(ix  ,iy  ,iz+1);
		double v010 = grid.read(ix  ,iy+1,iz  );
		double v100 = grid.read(ix+1,iy  ,iz  );
		double v011 = grid.read(ix  ,iy+1,iz+1);
		double v101 = grid.read(ix+1,iy  ,iz+1);
		double v110 = grid.read(ix+1,iy+1,iz  );
		double v111 = grid.read(ix+1,iy+1,iz+1);
		
		if (v000>=10000. ||
			v001>=10000. ||
			v010>=10000. ||
			v100>=10000. ||
			v011>=10000. ||
			v101>=10000. ||
			v110>=10000. ||
			v111>=10000.)
		{
			v000 = grid.read(ix  ,iy  ,iz  )<10000. ? v000 : isPointInside(ix  ,iy  ,iz  );
			v001 = grid.read(ix  ,iy  ,iz+1)<10000. ? v001 : isPointInside(ix  ,iy  ,iz+1);
			v010 = grid.read(ix  ,iy+1,iz  )<10000. ? v010 : isPointInside(ix  ,iy+1,iz  );
			v100 = grid.read(ix+1,iy  ,iz  )<10000. ? v100 : isPointInside(ix+1,iy  ,iz  );
			v011 = grid.read(ix  ,iy+1,iz+1)<10000. ? v011 : isPointInside(ix  ,iy+1,iz+1);
			v101 = grid.read(ix+1,iy  ,iz+1)<10000. ? v101 : isPointInside(ix+1,iy  ,iz+1);
			v110 = grid.read(ix+1,iy+1,iz  )<10000. ? v110 : isPointInside(ix+1,iy+1,iz  );
			v111 = grid.read(ix+1,iy+1,iz+1)<10000. ? v111 : isPointInside(ix+1,iy+1,iz+1);
		}
		
		double v00z = v000*(1.-dz) + v001*dz;
		double v01z = v010*(1.-dz) + v011*dz;
		double v10z = v100*(1.-dz) + v101*dz;
		double v11z = v110*(1.-dz) + v111*dz;
		
		double v0yz = v00z*(1.-dy) + v01z*dy;
		double v1yz = v10z*(1.-dy) + v11z*dy;
		
		return v0yz*(1.-dx) + v1yz*dx;
#else
		// TODO: use a better interpolation than nearest neighbor - seem to be of the same quality of the trilinear
		if (grid.read(ix,iy,iz) < 10000.)
		{
			return grid.read(ix,iy,iz);
		}
		else
		{
			return isPointInside(ix,iy,iz);
		}
#endif
	}
	else
		return 1.;
}

bool GeometryReader::isGeometryLoaded() const
{
	return bGeometryLoaded;
}

void GeometryReader::moveBody(double dt)
{
	properties.update(dt);
}
