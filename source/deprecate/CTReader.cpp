//
//  CTReader.cpp
//  CubismUP_3D
//
//	Reads a series of CT slices and computes the signed distance function
//
//	Created by Christian Conti on 08/15/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "CTReader.h"

void GeometryReaderCT::readFromFile(const int slice, const int startslice)
{
	stringstream filename;
	filename << "/cluster/work/infk/cconti/Seed/APL_087/t_" << setfill('0') << setw(4) << slice << ".bmp";
	
	cout << "loading file " << filename.str().c_str() << endl;
	FILE* f = fopen(filename.str().c_str(), "rb");
	
	unsigned char info[54];
	fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
	
	// extract image height and width from header
	int offset = *(int*)&info[10];
	sx = *(int*)&info[18];
	sy = *(int*)&info[22];
	short nbits = *(short*)&info[28];
	binary = !(nbits/8);
	nchannels = binary ? nbits : nbits/8; // handles binary images
	/*
	 cout << "header info:\n";
	 cout << "\tSize:\t" << sx << "x" << sy << endl;
	 cout << "\tDepth:\t" << nbits << endl;
	 cout << "\tOffset:\t" << offset << endl;
	 */
	if (binary)
	{
		cout << "1-bit format unsupported\n";
		abort();
	}
	
	// account for additional padding at the end of each line of pixels
	int sxp = ceil(sx/4.)*4;
	
	int size = nchannels * sxp * sy;
	unsigned char* data = new unsigned char[size+offset]; // allocate nchannels bytes per pixel
	fread(data, sizeof(unsigned char), size+offset, f); // read the rest of the data at once
	fclose(f);
	
	int isx = sxp;
	sx /= subsampling;
	sy /= subsampling;
	
	// create grid if first slice
	if (density.maxpoints==0)
	{
		cout << "Resetting grid size " << sx << " " << sy << " " << sz << endl;
		density.reset(sx, sy, sz);
		density.maxpoints = max(sx, max(sy,sz));
		gridsize = max(sx, max(sy,sz));
	}
	
	//cout << "Writing slice " << (slice-startslice)/subsampling << " of " << sz << endl;
#pragma omp parallel for
	for (int j=0; j<sy; j++)
		for (int i=0; i<sx; i++)
		{
			if (i*subsampling*nchannels+j*subsampling*nchannels*isx+offset >= size+offset)
				abort();
			double value = (double) data[i*subsampling*nchannels+j*subsampling*nchannels*isx+offset];
			
			if (i>=density.sizeX) cout << i << " " << density.sizeX << endl;
			if (j>=density.sizeY) cout << j << " " << density.sizeY << endl;
			if ((slice-startslice)/subsampling>=density.sizeZ) cout << (slice-startslice)/subsampling << " " << density.sizeZ << endl;
			density(i,j,(slice-startslice)/subsampling) = value;
		}
	
	delete [] data;
}

void GeometryReaderCT::parse(string filename)
{
	// need to get start, end and resolution from command line
	const int startslice = 128;
	int curslice = startslice;
	const int endslice = 8830;
	density.maxpoints = 0;
	subsampling = 4;
	
	sz = 1+(endslice-startslice)/subsampling;
	
	// read data
	// TODO: parallelize over slices - can be done as slices are independent
	//	be careful with common data set - sx, sy, grid initialization
	//	transform in for loop
	while (curslice < endslice)
	{
		//cout << "Loading slice " << curslice << endl;
		readFromFile(curslice, startslice);
		curslice += subsampling;
	}
	
	cout << "Data files read\n";
}

// constructor loading OBJ file
GeometryReaderCT::GeometryReaderCT(const string filename, Geometry::Properties & properties, int gridsize, const Real scaleFactor, const Geometry::Point transFactor) : GeometryReader(filename,properties,gridsize,scaleFactor,transFactor), buffer(5), subsampling(1)
{
	load(filename);
	sdf();
}

// clean up
GeometryReaderCT::~GeometryReaderCT()
{
}

void GeometryReaderCT::load(const string filename)
{
	if (bGeometryLoaded)
		cout << "WARNING - Overwriting previously loaded geometry\n";
	
	if (bVerbose)
		cout << "Loading data\n";
	
	stringstream fn;
	fn << filename;
	
	parse(fn.str());
	
	bGeometryLoaded = true;
}

void GeometryReaderCT::sdf(const Real scaleFactor, const Geometry::Point transFactor, double comx, double comy, double comz)
{
	if (!bGeometryLoaded)
	{
		cout << "SDF - Geometry not loaded yet!\n";
		abort();
	}
	
	if (bVerbose)
		cout << "Computing SDF\n";
	
	//===================================
	//	SDF computation
	//
	//	create a grid and mark cells which contain information as -1,
	//		empty cells as 1
	//		run a gaussian filter with appropriate size to diffuse
	//		to obtain a mollified Chi
	//	with this technique the Chi function is obtained directly
	//		without the need of computing a SDF
	//	requires the use of two grids
	//
	//	TODO: can be reduced to a single grid by first computing the center of mass
	//		and then transforming the values to -1/1 and filtering
	//
	//===================================
	
	// 0. compute grid size needed from size of bounding box
	density.maxsize = 1.;
	grid = density;
	grid.reset(density.sizeX+buffer*2, density.sizeY+buffer*2, density.sizeZ+buffer*2, 0);
	grid.maxpoints = max(grid.sizeX,max(grid.sizeY,grid.sizeZ));
	properties.minb.x = 0.;
	properties.minb.y = 0.;
	properties.minb.z = 0.;
	properties.maxb.x = (double)grid.sizeX/(double)grid.maxpoints;
	properties.maxb.y = (double)grid.sizeY/(double)grid.maxpoints;
	properties.maxb.z = (double)grid.sizeZ/(double)grid.maxpoints;
	properties.minb.x = properties.minb.x*scaleFactor + transFactor.x;
	properties.minb.y = properties.minb.y*scaleFactor + transFactor.y;
	properties.minb.z = properties.minb.z*scaleFactor + transFactor.z;
	properties.maxb.x = properties.maxb.x*scaleFactor + transFactor.x;
	properties.maxb.y = properties.maxb.y*scaleFactor + transFactor.y;
	properties.maxb.z = properties.maxb.z*scaleFactor + transFactor.z;
	grid.maxsize *= scaleFactor;
	
	// 1. setup of gaussian filter
	//	float filter[5] = {};
	
	// 2. filling of the grid with gaussian filter
#pragma omp parallel for
	for (int iz=0; iz<grid.sizeZ; iz++)
		for (int iy=0; iy<grid.sizeY; iy++)
			for (int ix=0; ix<grid.sizeX; ix++)
			{
				double r = 0.;
				
				for (int k=-2; k<3; k++)
					for (int j=-2; j<3; j++)
						for (int i=-2; i<3; i++)
							r += density.read(ix+i-buffer, iy+j-buffer, iz+k-buffer, 0.);
				
				grid(ix,iy,iz) += (r>0?-1.:1.);
			}
	
	
	if (bVerbose)
		cout << "Computing center of mass\n";
	
	Geometry::Point center(0,0,0);
	Geometry::Point centerOfMass(0,0,0);
	int count = 0;
	
	double dh = 1./(double)grid.maxpoints;
	for (int iz=0; iz<grid.sizeZ; iz++)
		for (int iy=0; iy<grid.sizeY; iy++)
			for (int ix=0; ix<grid.sizeX; ix++)
			{
				float sign = 0;
				for (int px=ix; px<grid.sizeX; px++)
				{
					if (grid(px,iy,iz)<10000.)
					{
						sign = grid(px,iy,iz)>0.?1.:-1.;
						break;
					}
				}
				if (sign==0.) sign = 1.;
				
				double centroid = sign<0?1.:0.;
				center.x += ((double)ix+.5)*dh * centroid;
				center.y += ((double)iy+.5)*dh * centroid;
				center.z += ((double)iz+.5)*dh * centroid;
				count += centroid;
				
				centerOfMass.x += ((double)ix+.5)*dh * density.read(ix-buffer,iy-buffer,iz-buffer, 0);
				centerOfMass.y += ((double)iy+.5)*dh * density.read(ix-buffer,iy-buffer,iz-buffer, 0);
				centerOfMass.z += ((double)iz+.5)*dh * density.read(ix-buffer,iy-buffer,iz-buffer, 0);
			}
	
	center.x /= count;
	center.y /= count;
	center.z /= count;
	
	centerOfMass.x /= count;
	centerOfMass.y /= count;
	centerOfMass.z /= count;
	
	double scalingsize = max(max(properties.maxb.x - properties.minb.x,properties.maxb.y - properties.minb.y),properties.maxb.z - properties.minb.z);
	center.x *= scalingsize;
	center.y *= scalingsize;
	center.z *= scalingsize;
	centerOfMass.x *= scalingsize;
	centerOfMass.y *= scalingsize;
	centerOfMass.z *= scalingsize;
	
	properties.com.x = centerOfMass.x + properties.minb.x;
	properties.com.y = centerOfMass.y + properties.minb.y;
	properties.com.z = centerOfMass.z + properties.minb.z;
	
	properties.centroid.x = center.x + properties.minb.x;
	properties.centroid.y = center.y + properties.minb.y;
	properties.centroid.z = center.z + properties.minb.z;
	
	bSDFComputed = true;
}


// check if a point lies inside the body
int GeometryReaderCT::isPointInside(int ix, int iy, int iz)
{
	cout << "isPointInside unnecessary - aborting()\n";
	abort();
}

// return signed distance to closest point on surface (neg. if inside the body).
double GeometryReaderCT::distance(double fx, double fy, double fz)
{
	// TODO:
	//	eventually remove intermediate grid and P2M directly on the output
	
	// rotate point to find position in rotated frame - should be w.r.t. center of mass
	// careful with left vs right multiplication!
	// careful with correct rotation direction
	double com[3];
	getCenterOfMass(com);
	
	// inverse angle
	double rfx = (fx-com[0])*properties.rotation[0][0] + (fy-com[1])*properties.rotation[0][1] + (fz-com[2])*properties.rotation[0][2] + com[0];
	double rfy = (fx-com[0])*properties.rotation[1][0] + (fy-com[1])*properties.rotation[1][1] + (fz-com[2])*properties.rotation[1][2] + com[1];
	double rfz = (fx-com[0])*properties.rotation[2][0] + (fy-com[1])*properties.rotation[2][1] + (fz-com[2])*properties.rotation[2][2] + com[2];
	
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
		
#if 0
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
			abort();
		}
#endif
	}
	else
		return 1.;
}
