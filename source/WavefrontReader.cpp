//
//  WavefrontReader.cpp
//  CubismUP_3D
//
//	Reads a wavefront object file and computes the signed distance function
//
//	Created by Christian Conti on 06/14/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "WavefrontReader.h"

void GeometryReaderOBJ::parse(string filename)
{
	vertices.clear();
	cout << "loading file " << filename.c_str() << endl;
	ifstream f(filename.c_str());
	
	properties.minb.x = properties.minb.y = properties.minb.z =  100.;
	properties.maxb.x = properties.maxb.y = properties.maxb.z = -100.;
	properties.com.x  = properties.com.y  = properties.com.z  = 0.;
	
	string type;
	f >> type;
	while (type != "vn")
		f >> type;
	
	while (type == "vn")
	{
		Geometry::Point n;
		f >> n.x;
		f >> n.y;
		f >> n.z;
		vnormals.push_back(n);
		
		f >> type;
		
		Geometry::Point p;
		f >> p.x;
		f >> p.y;
		f >> p.z;
		vertices.push_back(p);
		
		// quantity defined in mesh coordinate system
		properties.minb.x = min(properties.minb.x,p.x);
		properties.minb.y = min(properties.minb.y,p.y);
		properties.minb.z = min(properties.minb.z,p.z);
		properties.maxb.x = max(properties.maxb.x,p.x);
		properties.maxb.y = max(properties.maxb.y,p.y);
		properties.maxb.z = max(properties.maxb.z,p.z);
		
		f >> type;
	}
	
	while (type != "f")
		f >> type;
	
	while (type == "f" && !f.eof())
	{
		// TODO: unsafe - if the mesh is not only triangles, it might cause a problem!
		Geometry::Triangle t;
		
		f >> t.a;
		f >> type;
		f >> t.b;
		f >> type;
		f >> t.c;
		f >> type;
		triangles.push_back(t);
		
		f >> type;
	}
	
	for (int v=0; v<vertices.size(); v++)
	{
		vertices[v].x = vertices[v].x*baseScale + baseTranslation.x;
		vertices[v].y = vertices[v].y*baseScale + baseTranslation.y;
		vertices[v].z = vertices[v].z*baseScale + baseTranslation.z;
	}
	
	// transformation in global system of reference
	properties.minb.x = properties.minb.x*baseScale + baseTranslation.x;
	properties.minb.y = properties.minb.y*baseScale + baseTranslation.y;
	properties.minb.z = properties.minb.z*baseScale + baseTranslation.z;
	properties.maxb.x = properties.maxb.x*baseScale + baseTranslation.x;
	properties.maxb.y = properties.maxb.y*baseScale + baseTranslation.y;
	properties.maxb.z = properties.maxb.z*baseScale + baseTranslation.z;
}

// constructor loading OBJ file
GeometryReaderOBJ::GeometryReaderOBJ(const string filename, Geometry::Properties & properties, int gridsize, const Real scaleFactor, const Geometry::Point transFactor, const Real isosurface) : GeometryReader(filename,properties,gridsize,scaleFactor,transFactor), cgrid(NULL), gridlut(NULL), isosurface(isosurface)
{
	load(filename);
	sdf();
}

// clean up
GeometryReaderOBJ::~GeometryReaderOBJ()
{
	delete cgrid;
	delete [] gridlut;
}

void GeometryReaderOBJ::load(const string filename)
{
	if (bGeometryLoaded)
		cout << "WARNING - Overwriting previously loaded geometry\n";
	
	if (bVerbose)
		cout << "Loading data\n";
	
	stringstream fn;
	fn << filename;
	
	parse(fn.str());
	
	/*
	 // clean normals
	 for (int t=0; t<triangles.size(); t++)
	 {
		float dotab = vnormals[triangles[t].a].x * vnormals[triangles[t].b].x + vnormals[triangles[t].a].y * vnormals[triangles[t].b].y + vnormals[triangles[t].a].z * vnormals[triangles[t].b].z;
		float dotbc = vnormals[triangles[t].c].x * vnormals[triangles[t].b].x + vnormals[triangles[t].c].y * vnormals[triangles[t].b].y + vnormals[triangles[t].c].z * vnormals[triangles[t].b].z;
		float dotca = vnormals[triangles[t].a].x * vnormals[triangles[t].c].x + vnormals[triangles[t].a].y * vnormals[triangles[t].c].y + vnormals[triangles[t].a].z * vnormals[triangles[t].c].z;
	 
		if (dotab < 0 && dotca < 0)
		{
	 vnormals[triangles[t].a].x = -vnormals[triangles[t].a].x;
	 vnormals[triangles[t].a].y = -vnormals[triangles[t].a].y;
	 vnormals[triangles[t].a].z = -vnormals[triangles[t].a].z;
		}
		else if (dotab < 0 && dotbc < 0)
		{
	 vnormals[triangles[t].b].x = -vnormals[triangles[t].b].x;
	 vnormals[triangles[t].b].y = -vnormals[triangles[t].b].y;
	 vnormals[triangles[t].b].z = -vnormals[triangles[t].b].z;
		}
		else if (dotbc < 0 && dotca < 0)
		{
	 vnormals[triangles[t].c].x = -vnormals[triangles[t].c].x;
	 vnormals[triangles[t].c].y = -vnormals[triangles[t].c].y;
	 vnormals[triangles[t].c].z = -vnormals[triangles[t].c].z;
		}
	 }
	 */
	
	/*
	 // cleanup
	 vector<int> vUnref, vDuplicated;
	 
	 // TODO: remove unreferenced vertices
	 for (int v=0; v<vertices.size(); v++)
	 {
	 bool bFound = false;
	 for (int t=0; t<triangles.size(); t++)
	 if (triangles[t].a == v+1 || triangles[t].b == v+1 || triangles[t].c == v+1)
	 {
	 bFound = true;
	 break;
	 }
	 
	 if (bFound==false)
	 vUnref.push_back(v);
	 }
	 cout << "Removing " << vUnref.size() << " unreferenced vertices:\n";
	 for (int i=0; i<vUnref.size(); i++)
	 cout << "\t" << vUnref[i] << endl;
	 for (int i=vUnref.size()-1; i>=0; i--)
	 vertices.erase(vertices.begin()+vUnref[i]);
	 */
	/*
	 // TODO: remove duplicate vertices
	 //	issue: need to change triangle references as well
	 for (int v=0; v<vertices.size(); v++)
	 {
	 for (int v2=0; v2<vertices.size(); v2++)
	 if (v != v2 &&
	 vertices[v].x == vertices[v2].x &&
	 vertices[v].y == vertices[v2].y &&
	 vertices[v].z == vertices[v2].z)
	 vDuplicated.push_back(v2);
	 }
	 
	 cout << "Removing " << vDuplicated.size() << " duplicated vertices\n";
	 for (int i=vDuplicated.size()-1; i>=0; i--)
	 vertices.erase(vertices.begin()+vDuplicated[i]);
	 */
	
	
	cout << "parsed " << vertices.size() << " vertices\n";
	cout << "parsed " << triangles.size() << " triangles\n";
	bGeometryLoaded = true;
}

void GeometryReaderOBJ::sdf()
{
	if (!bGeometryLoaded)
	{
		cout << "SDF - Geometry not loaded yet!\n";
		abort();
	}
	
	if (bVerbose)
		cout << "Computing SDF\n";
	
	// 0. compute grid size needed from size of bounding box
	//	increase in size to account for cells with SDF
	//	is this enough?
	// compute maximum size of bounding box to scale seed
	float sizex = (properties.maxb.x - properties.minb.x);
	float sizey = (properties.maxb.y - properties.minb.y);
	float sizez = (properties.maxb.z - properties.minb.z);
	grid_maxsize = fmax(sizex,fmax(sizey,sizez));
	
	// compute aspect ratio
	float ratiox = sizex/grid_maxsize;
	float ratioy = sizey/grid_maxsize;
	float ratioz = sizez/grid_maxsize;
	
	// 1. initialize grid - this step has to be kept as it resets the values on the grid
	if (cgrid != NULL)
	{
		delete cgrid;
		cgrid = NULL;
	}
	const int BPDX = ceil(gridsize*ratiox/ScalarBlock::sizeX);
	const int BPDY = ceil(gridsize*ratioy/ScalarBlock::sizeY);
	const int BPDZ = ceil(gridsize*ratioz/ScalarBlock::sizeZ);
	
	//cout << "Seed grid is " << BPDX*ScalarBlock::sizeX << " " << BPDY*ScalarBlock::sizeY << " " << BPDZ*ScalarBlock::sizeZ << endl;
	
	cgrid = new GScalar(BPDX, BPDY, BPDZ);
	gridlut = new int[BPDX*BPDY*BPDZ];
	for (int i=0; i<BPDX*BPDY*BPDZ; i++) gridlut[i] = 3;
	
	vector<BlockInfo> vInfo = cgrid->getBlocksInfo();
	
	// 2. P2M with abs(x)
	// TODO: need to set a smaller value!
	const int supportStart = -19;
	const int supportEnd = 20;
	
	double invGridMax = 1./grid_maxsize;
	double dh = 1./(double)gridsize;
	
	Geometry::Point translation(.125*(properties.maxb.x-properties.minb.x)*invGridMax,
								.125*(properties.maxb.y-properties.minb.y)*invGridMax,
								.125*(properties.maxb.z-properties.minb.z)*invGridMax);
	
#pragma omp parallel
	{
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
		const int mynode = omp_get_thread_num() / cores_per_node;
		numa_run_on_node(mynode);
#endif
		
		
#pragma omp for
		for (int i=0; i<vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			ScalarBlock& b = *(ScalarBlock*)info.ptrBlock;
			
			for (int iz=0; iz<ScalarBlock::sizeZ; iz++)
				for (int iy=0; iy<ScalarBlock::sizeY; iy++)
					for (int ix=0; ix<ScalarBlock::sizeX; ix++)
						b(ix,iy,iz) = 10000.;
			
			for (int v=0; v<vertices.size(); v++)
			{
				// scale seed
				// shift needed to center seed in the grid
				Geometry::Point pos(.75*(vertices[v].x - properties.minb.x)*invGridMax + translation.x,
									.75*(vertices[v].y - properties.minb.y)*invGridMax + translation.y,
									.75*(vertices[v].z - properties.minb.z)*invGridMax + translation.z);
				
				// get index of closest gridpoint
				const int ipx = floor(pos.x*gridsize);
				const int ipy = floor(pos.y*gridsize);
				const int ipz = floor(pos.z*gridsize);
				
				
#if 0
				int startx = 0;
				int starty = 0;
				int startz = 0;
				int endx   = ScalarBlock::sizeX;
				int endy   = ScalarBlock::sizeY;
				int endz   = ScalarBlock::sizeZ;
#else
				const int sX = ScalarBlock::sizeX;
				const int sY = ScalarBlock::sizeY;
				const int sZ = ScalarBlock::sizeZ;
				int startx = max(ipx+supportStart-info.index[0]*sX,0);
				int starty = max(ipy+supportStart-info.index[1]*sY,0);
				int startz = max(ipz+supportStart-info.index[2]*sZ,0);
				int endx   = min(ipx+supportEnd  -info.index[0]*sX,sX);
				int endy   = min(ipy+supportEnd  -info.index[1]*sY,sY);
				int endz   = min(ipz+supportEnd  -info.index[2]*sZ,sZ);
#endif
				
				// set lut for boundary blocks
				if (startx<endx && starty<endy && startz<endz)
				{
					gridlut[i] = 1;
					startx = starty = startz = 0;
					endx = ScalarBlock::sizeX;
					endy = ScalarBlock::sizeY;
					endz = ScalarBlock::sizeZ;
				}
				
				for (int iz=startz; iz<endz; iz++)
					for (int iy=starty; iy<endy; iy++)
						for (int ix=startx; ix<endx; ix++)
						{
							Real p[3];
							info.pos(p,ix,iy,iz);
							
							float dx = p[0] - pos.x;
							float dy = p[1] - pos.y;
							float dz = p[2] - pos.z;
							
							// distance
							float dist = sqrt(dx*dx + dy*dy + dz*dz);
							
							// sign - there is an issue with flipped normals which creates an error in the SDF!
							float dot = dx*vnormals[v].x + dy*vnormals[v].y + dz*vnormals[v].z;
							float sign = dot>0. ? 1. : -1.;
							
							bool bBest = dist < fabs(b(ix,iy,iz));
							b(ix,iy,iz)       = bBest ? sign*dist : b(ix,iy,iz);
						}
			}
		}
	}
	
	//stringstream ss;
	//ss << "/cluster/scratch_xp/public/cconti/CubismUP/sdf.vti";
	//dumper.Write(*cgrid, ss.str());
	
	if (bVerbose)
		cout << "LUT computation\n";
	
	for (int k=0; k<BPDZ; k++)
		for (int j=0; j<BPDY; j++)
			for (int i=0; i<BPDX; i++)
				if (gridlut[i+j*BPDX+k*BPDX*BPDY] == 3)
				{
					// if a boundary is encountered, assume that the block is an "outside block"
					// else it's inside (no cavities are allowed)
					bool xp=false;
					bool xm=false;
					bool yp=false;
					bool ym=false;
					bool zp=false;
					bool zm=false;
					
					
					// xp
					for (int x=i; x<BPDX; x++)
						if (gridlut[x+j*BPDX+k*BPDX*BPDY] == 1)
							xp = true;
					
					if (xp==false)
					{
						gridlut[i+j*BPDX+k*BPDX*BPDY] = 0;
						continue;
					}
					
					// xm
					for (int x=i; x>=0; x--)
						if (gridlut[x+j*BPDX+k*BPDX*BPDY] == 1)
							xm = true;
					
					if (xm==false)
					{
						gridlut[i+j*BPDX+k*BPDX*BPDY] = 0;
						continue;
					}
					
					
					// yp
					for (int y=j; y<BPDY; y++)
						if (gridlut[i+y*BPDX+k*BPDX*BPDY] == 1)
							yp = true;
					
					if (yp==false)
					{
						gridlut[i+j*BPDX+k*BPDX*BPDY] = 0;
						continue;
					}
					
					// ym
					for (int y=j; y>=0; y--)
						if (gridlut[i+y*BPDX+k*BPDX*BPDY] == 1)
							ym = true;
					
					if (ym==false)
					{
						gridlut[i+j*BPDX+k*BPDX*BPDY] = 0;
						continue;
					}
					
					
					// zp
					for (int z=k; z<BPDZ; z++)
						if (gridlut[i+j*BPDX+z*BPDX*BPDY] == 1)
							zp = true;
					
					if (zp==false)
					{
						gridlut[i+j*BPDX+k*BPDX*BPDY] = 0;
						continue;
					}
					
					// zm
					for (int z=k; z>=0; z--)
						if (gridlut[i+j*BPDX+z*BPDX*BPDY] == 1)
							zm = true;
					
					if (zm==false)
					{
						gridlut[i+j*BPDX+k*BPDX*BPDY] = 0;
						continue;
					}
					
					gridlut[i+j*BPDX+k*BPDX*BPDY] = 2;
				}
	
	if (bVerbose)
		cout << "Computing center of mass\n";
	
	double centerx = 0, centery = 0, centerz = 0;
	double count = 0;
	
#pragma omp parallel for reduction(+:centerx) reduction(+:centery) reduction(+:centerz) reduction(+:count)
	for (int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		ScalarBlock& b = *(ScalarBlock*)info.ptrBlock;
		
		for (int iz=0; iz<ScalarBlock::sizeZ; iz++)
			for (int iy=0; iy<ScalarBlock::sizeY; iy++)
				for (int ix=0; ix<ScalarBlock::sizeX; ix++)
				{
					// need to be reloaded for every point as it is updated while iterating
					const int lutvalue = gridlut[info.index[0] + cgrid->getBlocksPerDimension(0) * (info.index[1] + cgrid->getBlocksPerDimension(1) * info.index[2])];
					
					Real p[3];
					info.pos(p,ix,iy,iz);
					
					float sign = 0;
					
					if (lutvalue==2)
					{
						sign = -1;
						b(ix,iy,iz) = -1;
					}
					else if (lutvalue==0)
					{
						sign = 1;
						b(ix,iy,iz) = 1;
					}
					else if (lutvalue==1)
					{
						sign = b(ix,iy,iz);//findInOut(ix,iy,iz,info.index[0],info.index[1],info.index[2]);
						//b(ix,iy,iz) = sign;
					}
					else
					{
						cout << "if you are here, something went wrong with the lookup table!\n";
						abort();
					}
					
					// Assumes constant density
					double volume = sign<0?1.:0.;
					centerx += p[0] * volume;
					centery += p[1] * volume;
					centerz += p[2] * volume;
					count += volume;
				}
	}
	
	centerx /= count;
	centery /= count;
	centerz /= count;
	
	double scalingsize = max(max(properties.maxb.x - properties.minb.x,properties.maxb.y - properties.minb.y),properties.maxb.z - properties.minb.z);
	cout << centerx << " " << centery << " " << centerz << endl;
	centerx *= scalingsize;
	centery *= scalingsize;
	centerz *= scalingsize;
	/*
	if (bVerbose)
		cout << "Computing inertia matrix\n";
	
	double J0 = 0, J1 = 0, J2 = 0, J3 = 0, J4 = 0, J5 = 0;
	
#pragma omp parallel for reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5)
	for (int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		ScalarBlock& b = *(ScalarBlock*)info.ptrBlock;
		
		for (int iz=0; iz<ScalarBlock::sizeZ; iz++)
			for (int iy=0; iy<ScalarBlock::sizeY; iy++)
				for (int ix=0; ix<ScalarBlock::sizeX; ix++)
				{
					// need to be reloaded for every point as it is updated while iterating
					const int lutvalue = gridlut[info.index[0] + cgrid->getBlocksPerDimension(0) * (info.index[1] + cgrid->getBlocksPerDimension(1) * info.index[2])];
					
					Real p[3];
					info.pos(p,ix,iy,iz);
					p[0] = p[0]*scalingsize - centerx;
					p[1] = p[1]*scalingsize - centery;
					p[2] = p[2]*scalingsize - centerz;
					
					float sign = 0;
					
					if (lutvalue==2)
						sign = -1;
					else if (lutvalue==0)
						sign = 1;
					else if (lutvalue==1)
						sign = b(ix,iy,iz);
					else
					{
						cout << "if you are here, something went wrong with the lookup table!\n";
						abort();
					}
					
					// Assumes constant density
					double volume = sign<0?1.:0.;
					
					J0 += volume * (p[1]*p[1] + p[2]*p[2]); //       y^2 + z^2
					J1 += volume * (p[0]*p[0] + p[2]*p[2]); // x^2 +     + z^2
					J2 += volume * (p[0]*p[0] + p[1]*p[1]); // x^2 + y^2
					J3 -= volume * p[0] * p[1]; // xy
					J4 -= volume * p[0] * p[2]; // xz
					J5 -= volume * p[1] * p[2]; // yz
				}
	}
	*/
	//const string filename_sdf = "seed_sdf";
	//const string filepath_sdf = "/cluster/scratch_xl/public/cconti/";
	//DumpHDF5<ScalarGrid, StreamerMesh>(*cgrid, 0, filename_sdf, filepath_sdf);
	
	/*
	for (int k=0; k<BPDZ; k++)
	{
		for (int j=0; j<BPDY; j++)
		{
			for (int i=0; i<BPDX; i++)
				cout << gridlut[i+j*BPDX+k*BPDX*BPDY] << " ";
			cout << endl;
		}
		cout << endl;
	}
	*/
	/*
	const double factor = 1;//128.*128.*128.*128.*128.;
	properties.J0[0] = J0 * scalingsize*scalingsize*scalingsize/factor;
	properties.J0[1] = J1 * scalingsize*scalingsize*scalingsize/factor;
	properties.J0[2] = J2 * scalingsize*scalingsize*scalingsize/factor;
	properties.J0[3] = J3 * scalingsize*scalingsize*scalingsize/factor;
	properties.J0[4] = J4 * scalingsize*scalingsize*scalingsize/factor;
	properties.J0[5] = J5 * scalingsize*scalingsize*scalingsize/factor;
	
	cout << scalingsize << endl;
	cout << centerx << " " << centery << " " << centerz << endl;
	cout << properties.J0[0] << " " << J0 << endl;
	cout << properties.J0[1] << " " << J1 << endl;
	cout << properties.J0[2] << " " << J2 << endl;
	cout << properties.J0[3] << " " << J3 << endl;
	cout << properties.J0[4] << " " << J4 << endl;
	cout << properties.J0[5] << " " << J5 << endl;
	*/
	
	// center of mass, global coordinates
	properties.com.x = centerx + properties.minb.x;
	properties.com.y = centery + properties.minb.y;
	properties.com.z = centerz + properties.minb.z;
	
	// centroid, global coordinates
	properties.centroid.x = centerx + properties.minb.x;
	properties.centroid.y = centery + properties.minb.y;
	properties.centroid.z = centerz + properties.minb.z;
	
	//cout << "SDF done\n";
	bSDFComputed = true;
}

// check if a point lies inside the body
int GeometryReaderOBJ::isPointInside(int ix, int iy, int iz)
{
	cout << "This is not needed\n";
	abort();
	
	
	if (!bGeometryLoaded)
	{
		cout << "isPointInside - Geometry not loaded yet!\n";
		abort();
	}
	
	float signx = 0;
	int bx = ix/ScalarBlock::sizeX;
	int by = iy/ScalarBlock::sizeY;
	int bz = iz/ScalarBlock::sizeZ;
	ix = ix%ScalarBlock::sizeX;
	iy = iy%ScalarBlock::sizeY;
	iz = iz%ScalarBlock::sizeZ;
	
	// this should not happen
	if (bx >= cgrid->getBlocksPerDimension(0) ||
		by >= cgrid->getBlocksPerDimension(1) ||
		bz >= cgrid->getBlocksPerDimension(2))
	{
		return 1.;
	}
	
	const int lutvalue = gridlut[bx + cgrid->getBlocksPerDimension(0) * (by + cgrid->getBlocksPerDimension(1) * bz)];
	// outside
	if (lutvalue == 0)
		signx=1.;
	// inside
	else if (lutvalue == 2)
		signx=-1.;
	// boundary
	else if (lutvalue == 1)
		signx = (*cgrid)(bx,by,bz)(ix,iy,iz);//findInOut(ix, iy, iz, bx, by, bz);
	// should not happen
	else
	{
		cout << "This should not happen (" << lutvalue << " for block " << bx << " " << by << " " << bz << ")!\n";
		abort();
	}
	
	//(*cgrid)(bx,by,bz)(ix,iy,iz) = signx;
	return (*cgrid)(bx,by,bz)(ix,iy,iz);
}

// return signed distance to closest point on surface (neg. if inside the body).
double GeometryReaderOBJ::distance(double fx, double fy, double fz)
{
	// this might be problematic as it is called from an MPI process and requires communication with another
	
	// rotate point to find position in rotated frame - should be w.r.t. center of mass
	// careful with left vs right multiplication!
	// careful with correct rotation direction
	double com[3];
	getCenterOfMass(com);
	
	//cout << properties.rotation[0][0] << " " << properties.rotation[0][1] << " " << properties.rotation[0][2] << endl;
	//cout << properties.rotation[1][0] << " " << properties.rotation[1][1] << " " << properties.rotation[1][2] << endl;
	//cout << properties.rotation[2][0] << " " << properties.rotation[2][1] << " " << properties.rotation[2][2] << endl;
	//abort();
	
	//double rfx = (fx-com[0])*properties.rotation[0][0] + (fy-com[1])*properties.rotation[0][1] + (fz-com[2])*properties.rotation[0][2] + com[0];
	//double rfy = (fx-com[0])*properties.rotation[1][0] + (fy-com[1])*properties.rotation[1][1] + (fz-com[2])*properties.rotation[1][2] + com[1];
	//double rfz = (fx-com[0])*properties.rotation[2][0] + (fy-com[1])*properties.rotation[2][1] + (fz-com[2])*properties.rotation[2][2] + com[2];
	double rfx = (fx-com[0])*properties.rotation[0][0] + (fy-com[1])*properties.rotation[1][0] + (fz-com[2])*properties.rotation[2][0] + com[0];
	double rfy = (fx-com[0])*properties.rotation[0][1] + (fy-com[1])*properties.rotation[1][1] + (fz-com[2])*properties.rotation[2][1] + com[1];
	double rfz = (fx-com[0])*properties.rotation[0][2] + (fy-com[1])*properties.rotation[1][2] + (fz-com[2])*properties.rotation[2][2] + com[2];
	
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
		int ix = (int)((rfx-properties.minb.x)/(double)grid_maxsize*gridsize);
		int iy = (int)((rfy-properties.minb.y)/(double)grid_maxsize*gridsize);
		int iz = (int)((rfz-properties.minb.z)/(double)grid_maxsize*gridsize);
		
#if 0
		// trilinear interpolation - should be ok since the signed distance function is a linear function
		int dx = (rfx-properties.minb.x)/(double)grid_maxsize*gridsize - ix;
		int dy = (rfy-properties.minb.y)/(double)grid_maxsize*gridsize - iy;
		int dz = (rfz-properties.minb.z)/(double)grid_maxsize*gridsize - iz;
		
		// local coordinates and block coordinates
		int bx = ix/ScalarBlock::sizeX;
		int by = iy/ScalarBlock::sizeY;
		int bz = iz/ScalarBlock::sizeZ;
		int ixl = ix%ScalarBlock::sizeX;
		int iyl = iy%ScalarBlock::sizeY;
		int izl = iz%ScalarBlock::sizeZ;
		int bxp1 = (ix+1)==ScalarBlock::sizeX ? bx+1 : bx;
		int byp1 = (iy+1)==ScalarBlock::sizeY ? by+1 : by;
		int bzp1 = (iz+1)==ScalarBlock::sizeZ ? bz+1 : bz;
		
		double v000 = (*cgrid)(bx  ,by  ,bz  )(ixl  ,iyl  ,izl  );
		double v001 = (*cgrid)(bx  ,by  ,bzp1)(ixl  ,iyl  ,izl+1);
		double v010 = (*cgrid)(bx  ,byp1,bz  )(ixl  ,iyl+1,izl  );
		double v100 = (*cgrid)(bxp1,by  ,bz  )(ixl+1,iyl  ,izl  );
		double v011 = (*cgrid)(bx  ,byp1,bzp1)(ixl  ,iyl+1,izl+1);
		double v101 = (*cgrid)(bxp1,by  ,bzp1)(ixl+1,iyl  ,izl+1);
		double v110 = (*cgrid)(bxp1,byp1,bz  )(ixl+1,iyl+1,izl  );
		double v111 = (*cgrid)(bxp1,byp1,bzp1)(ixl+1,iyl+1,izl+1);
		
		if (v000>=10000. ||
			v001>=10000. ||
			v010>=10000. ||
			v100>=10000. ||
			v011>=10000. ||
			v101>=10000. ||
			v110>=10000. ||
			v111>=10000.)
		{
			v000 = (*cgrid)(bx  ,by  ,bz  )(ixl  ,iyl  ,izl  )<10000. ? v000 : isPointInside(ix  ,iy  ,iz  );
			v001 = (*cgrid)(bx  ,by  ,bzp1)(ixl  ,iyl  ,izl+1)<10000. ? v001 : isPointInside(ix  ,iy  ,iz+1);
			v010 = (*cgrid)(bx  ,byp1,bz  )(ixl  ,iyl+1,izl  )<10000. ? v010 : isPointInside(ix  ,iy+1,iz  );
			v100 = (*cgrid)(bxp1,by  ,bz  )(ixl+1,iyl  ,izl  )<10000. ? v100 : isPointInside(ix+1,iy  ,iz  );
			v011 = (*cgrid)(bx  ,byp1,bzp1)(ixl  ,iyl+1,izl+1)<10000. ? v011 : isPointInside(ix  ,iy+1,iz+1);
			v101 = (*cgrid)(bxp1,by  ,bzp1)(ixl+1,iyl  ,izl+1)<10000. ? v101 : isPointInside(ix+1,iy  ,iz+1);
			v110 = (*cgrid)(bxp1,byp1,bz  )(ixl+1,iyl+1,izl  )<10000. ? v110 : isPointInside(ix+1,iy+1,iz  );
			v111 = (*cgrid)(bxp1,byp1,bzp1)(ixl+1,iyl+1,izl+1)<10000. ? v111 : isPointInside(ix+1,iy+1,iz+1);
		}
		
		double v00z = v000*(1.-dz) + v001*dz;
		double v01z = v010*(1.-dz) + v011*dz;
		double v10z = v100*(1.-dz) + v101*dz;
		double v11z = v110*(1.-dz) + v111*dz;
		
		double v0yz = v00z*(1.-dy) + v01z*dy;
		double v1yz = v10z*(1.-dy) + v11z*dy;
		
		return v0yz*(1.-dx) + v1yz*dx;
#else
		int bx = ix/ScalarBlock::sizeX;
		int by = iy/ScalarBlock::sizeY;
		int bz = iz/ScalarBlock::sizeZ;
		int ixl = ix%ScalarBlock::sizeX;
		int iyl = iy%ScalarBlock::sizeY;
		int izl = iz%ScalarBlock::sizeZ;
		
		if ((*cgrid)(bx,by,bz)(ixl,iyl,izl) < 10000.)
			return (*cgrid)(bx,by,bz)(ixl,iyl,izl);
		else
			return isPointInside(ix,iy,iz);
#endif
	}
	else
		return 1.;
}

