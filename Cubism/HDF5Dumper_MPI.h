/*
 *  HDF5Dumper_MPI.h
 *  Cubism
 *
 *  Created by Babak Hejazialhosseini on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

 #ifndef Cubism_HDF5Dumper_MPI
 #define Cubism_HDF5Dumper_MPI

#include <cassert>
#include <cstdio>

#ifdef _USE_HDF_
#include <hdf5.h>
#endif

//bad:
#include "../source/GenericCoordinator.h"
#include "../source/GenericOperator.h"

//#ifdef _SP_COMP_ _FLOAT_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
//#else
//#define HDF_REAL H5T_NATIVE_DOUBLE
//#endif

using namespace std;

#include "BlockInfo.h"

class OperatorLoad : public GenericLabOperator
{
	float* const data;
	const unsigned long int NBX, NBY, NBZ, NX, NY, NZ, NCHANNELS;
 public:
	OperatorLoad(float* const dump_data, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz, const Real sliceZ)
	: data(dump_data), NBX(nx), NBY(ny), NBZ(nz), NX(nx*FluidBlock::sizeX), NY(ny*FluidBlock::sizeY), NZ(nz*FluidBlock::sizeZ), NCHANNELS(StreamerHDF5::NCHANNELS)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0,1,2);
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& b) const
	{
		const Real inv2h = .5 / info.h_gridpoint;
		const unsigned long int idx[3] = {
				info.index[0]%NBX,
				info.index[1]%NBY,
				info.index[2]%NBZ
		};
		StreamerHDF5 streamer(b);

		for(unsigned long int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(unsigned long int iy=0; iy<FluidBlock::sizeY; iy++)
		for(unsigned long int ix=0; ix<FluidBlock::sizeX; ix++){
			const unsigned long int gx = idx[0]*FluidBlock::sizeX + ix;
			const unsigned long int gy = idx[1]*FluidBlock::sizeY + iy;
			const unsigned long int gz = idx[2]*FluidBlock::sizeZ + iz;
			const unsigned long int idx = NCHANNELS * (gx + NX * (gy + NY * gz));
			assert(idx < NX * NY * NZ * NCHANNELS);
			float * const ptr = data + idx;
			//assert(NCHANNELS*(gz + NZ * (gy + NY * gx)) < NX * NY * NZ * NCHANNELS);
			//float* const ptr = data + NCHANNELS*(gz + NZ * (gy + NY * gx));

			FluidElement& phiW = lab(ix-1,iy  ,iz  );
			FluidElement& phiE = lab(ix+1,iy  ,iz  );
			FluidElement& phiS = lab(ix  ,iy-1,iz  );
			FluidElement& phiN = lab(ix  ,iy+1,iz  );
			FluidElement& phiF = lab(ix  ,iy  ,iz-1);
			FluidElement& phiB = lab(ix  ,iy  ,iz+1);
			const Real vorticX = inv2h * ((phiN.w-phiS.w) - (phiB.v-phiF.v));
			const Real vorticY = inv2h * ((phiB.u-phiF.u) - (phiE.w-phiW.w));
			const Real vorticZ = inv2h * ((phiE.v-phiW.v) - (phiN.u-phiS.u));
			//b(ix,iy,iz).tmpU = vorticX;
			//b(ix,iy,iz).tmpV = vorticY;
			//b(ix,iy,iz).tmpW = vorticZ;

			Real output[NCHANNELS];
			for(int i=0; i<NCHANNELS; ++i) output[i] = 0;
			streamer.operate(ix, iy, iz, (Real*)output);

			const Real strainXY = inv2h * ((phiE.v-phiW.v) + (phiN.u-phiS.u)) * 0.5;
			const Real strainYZ = inv2h * ((phiN.w-phiS.w) + (phiB.v-phiF.v)) * 0.5;
			const Real strainZX = inv2h * ((phiB.u-phiF.u) + (phiE.w-phiW.w)) * 0.5;
			const Real strainXX = inv2h * (phiE.u-phiW.u);
			const Real strainYY = inv2h * (phiN.v-phiS.v);
			const Real strainZZ = inv2h * (phiB.w-phiF.w);

			const Real OO = 0.5*(vorticX*vorticX+vorticY*vorticY+vorticZ*vorticZ);
			const Real SS = strainXX*strainXX+strainYY*strainYY+strainZZ*strainZZ+  //
							     2*(strainXY*strainXY+strainYZ*strainYZ+strainZX*strainZX);

			for(int i=0; i<NCHANNELS; ++i) ptr[i] = (float)output[i];
			ptr[NCHANNELS-1] = float(.5*(OO-SS));
		}
	}
};

class OperatorLoadFlat : public GenericLabOperator
{
	float* const data;
	const unsigned long int NBX, NBY, NX, NY, NCHANNELS;
	const Real sliceZ;
 public:
	OperatorLoadFlat(float* const dump_data, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz, const Real sliceZ)
	: data(dump_data), NBX(nx), NBY(ny), NX(nx*FluidBlock::sizeX), NY(ny*FluidBlock::sizeY), NCHANNELS(StreamerHDF5::NCHANNELS), sliceZ(sliceZ)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0,1,2);
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& b) const
	{
		const Real halfH = 0.5*info.h_gridpoint;
		Real ini[3], fin[3];
		info.pos(ini, 0, 0, 0);
		info.pos(fin, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
		if (ini[2]>sliceZ+halfH || fin[2]<sliceZ-halfH) return;
		int mid = (int)std::floor((sliceZ-ini[2])/info.h_gridpoint);
		if (mid<0) mid=0;
		if (mid>=FluidBlock::sizeZ) mid=FluidBlock::sizeZ-1;
		unsigned long int iz = mid; //probably useless

		const Real inv2h = .5 / info.h_gridpoint;
		const unsigned long int idx[2] = {
				info.index[0]%NBX,
				info.index[1]%NBY
		};
		StreamerHDF5 streamer(b);

		for(unsigned long int ix=0; ix<FluidBlock::sizeX; ix++) {
			const unsigned long int gx = idx[0]*FluidBlock::sizeX + ix;
			for(unsigned long int iy=0; iy<FluidBlock::sizeY; iy++) {
				const unsigned long int gy = idx[1]*FluidBlock::sizeY + iy;
				assert(NCHANNELS*(gy + NY * gx) < NX * NY * NCHANNELS);
				float* const ptr = data + NCHANNELS*(gy + NY * gx);

				FluidElement& phiW = lab(ix-1,iy  ,iz  );
				FluidElement& phiE = lab(ix+1,iy  ,iz  );
				FluidElement& phiS = lab(ix  ,iy-1,iz  );
				FluidElement& phiN = lab(ix  ,iy+1,iz  );
				FluidElement& phiF = lab(ix  ,iy  ,iz-1);
				FluidElement& phiB = lab(ix  ,iy  ,iz+1);
				const Real vorticX = inv2h * ((phiN.w-phiS.w) - (phiB.v-phiF.v));
				const Real vorticY = inv2h * ((phiB.u-phiF.u) - (phiE.w-phiW.w));
				const Real vorticZ = inv2h * ((phiE.v-phiW.v) - (phiN.u-phiS.u));
				//b(ix,iy,iz).tmpU = vorticX;
				//b(ix,iy,iz).tmpV = vorticY;
				//b(ix,iy,iz).tmpW = vorticZ;

				Real output[NCHANNELS];
				for(int i=0; i<NCHANNELS; ++i) output[i] = 0;
				streamer.operate(ix, iy, iz, (Real*)output);

				const Real strainXY = inv2h * ((phiE.v-phiW.v) + (phiN.u-phiS.u)) * 0.5;
				const Real strainYZ = inv2h * ((phiN.w-phiS.w) + (phiB.v-phiF.v)) * 0.5;
				const Real strainZX = inv2h * ((phiB.u-phiF.u) + (phiE.w-phiW.w)) * 0.5;
				const Real strainXX = inv2h * (phiE.v-phiW.u);
				const Real strainYY = inv2h * (phiN.v-phiS.v);
				const Real strainZZ = inv2h * (phiB.w-phiF.w);
				const Real OO = .25*(vorticX*vorticX+vorticY*vorticY+vorticZ*vorticZ);
				const Real SS = strainXX*strainXX+strainYY*strainYY+strainZZ*strainZZ+  //
								strainXY*strainXY+strainYZ*strainYZ+strainZX*strainZX;

				for(int i=0; i<NCHANNELS; ++i) ptr[i] = (float)output[i];
				ptr[NCHANNELS-1] = float(.5*(OO-SS));
			}
		}
	}
};

template<typename Loader>
class CoordinatorLoad : public GenericCoordinator
{
	float* const data;
 public:
	CoordinatorLoad(FluidGridMPI *grid, float* const dump_data) : GenericCoordinator(grid), data(dump_data) { }
	void operator()(const Real dt)
	{
		const Real maxExt = grid->getBlocksPerDimension(0)*FluidBlock::sizeX;
		const Real zExt = grid->getBlocksPerDimension(2)*FluidBlock::sizeZ;
		const Real sliceZ = 0.5*zExt/maxExt;
		const unsigned long int NX = grid->getResidentBlocksPerDimension(0);
		const unsigned long int NY = grid->getResidentBlocksPerDimension(1);
		const unsigned long int NZ = grid->getResidentBlocksPerDimension(2);
		Loader kernel(data,NX,NY,NZ,sliceZ);
		compute(kernel);
	}
	string getName() { return "Load"; }
};

static void DumpHDF5_MPI(FluidGridMPI &grid, const Real absTime, const string f_name, const string dump_path="./")
{
	#ifdef _USE_HDF_
	typedef typename FluidGridMPI::BlockType B;

	int rank;
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

  MPI_Comm comm = grid.getCartComm();
	MPI_Comm_rank(comm, &rank);

	int coords[3];
	grid.peindex(coords);

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const unsigned long int NCHANNELS = StreamerHDF5::NCHANNELS;
	float* array_all = new float[NX * NY * NZ * NCHANNELS];

	hsize_t count[4] = {
		grid.getResidentBlocksPerDimension(2)*B::sizeZ,
		grid.getResidentBlocksPerDimension(1)*B::sizeY,
		grid.getResidentBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t dims[4] = {
		grid.getBlocksPerDimension(2)*B::sizeZ,
		grid.getBlocksPerDimension(1)*B::sizeY,
		grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t offset[4] = {
		coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
		coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
		coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());

	H5open();
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	status = H5Pclose(fapl_id);

	CoordinatorLoad<OperatorLoad> coord(&grid, array_all);
	coord(0.);

	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

	fspace_id = H5Screate_simple(4, dims, NULL);
	#ifndef _ON_FERMI_
		dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	#else
		dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	#endif

	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	mspace_id = H5Screate_simple(4, count, NULL);
	status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

	status = H5Sclose(mspace_id);
	status = H5Sclose(fspace_id);
	status = H5Dclose(dataset_id);
	status = H5Pclose(fapl_id);
	status = H5Fclose(file_id);
	H5close();

	delete [] array_all;

	if (rank==0)
	{
		char wrapper[256];
		sprintf(wrapper, "%s/%s.xmf", dump_path.c_str(), f_name.c_str());
		FILE *xmf = 0;
		xmf = fopen(wrapper, "w");
		fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
		fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
		fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
		fprintf(xmf, " <Domain>\n");
		fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
		//fprintf(xmf, "     <Time Value=\"%05d\"/>\n", iCounter);
		fprintf(xmf, "     <Time Value=\"%e\"/>\n", absTime);
		fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n",
																			(int)dims[0], (int)dims[1], (int)dims[2]);
		fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
		fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", 0.,0.,0.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", grid.getH(), grid.getH(), grid.getH());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Geometry>\n");

		fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Node\">\n",
																							StreamerHDF5::getAttributeName());
		fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
												(int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3]);
		fprintf(xmf, "        %s:/data\n",(f_name+".h5").c_str());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Attribute>\n");

		fprintf(xmf, "   </Grid>\n");
		fprintf(xmf, " </Domain>\n");
		fprintf(xmf, "</Xdmf>\n");
		fclose(xmf);
	}
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}

static void DumpHDF5flat_MPI(FluidGridMPI &grid, const Real absTime, const string f_name, const string dump_path="./")
{
	#ifdef _USE_HDF_
	typedef typename FluidGridMPI::BlockType B;
	int rank;
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;
	MPI_Comm_rank(grid.getCartComm(), &rank);
	int coords[3];
	grid.peindex(coords);

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = 1;
	static const unsigned long int NCHANNELS = StreamerHDF5::NCHANNELS;
	float* array_all = new float[NX * NY * NCHANNELS];
	vector<BlockInfo> vInfo = grid.getBlocksInfo();

	hsize_t count[4] = {
		grid.getResidentBlocksPerDimension(0)*B::sizeX,
		grid.getResidentBlocksPerDimension(1)*B::sizeY, 1, NCHANNELS};
	hsize_t dims[4] = {
		grid.getBlocksPerDimension(0)*B::sizeX,
		grid.getBlocksPerDimension(1)*B::sizeY, 1, NCHANNELS};
	hsize_t offset[4] = {
		coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX,
		coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY, 0, 0};

	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());

	H5open();
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	status = H5Pset_fapl_mpio(fapl_id, grid.getCartComm(), MPI_INFO_NULL);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	status = H5Pclose(fapl_id);

	CoordinatorLoad<OperatorLoadFlat> coord(&grid, array_all);
	coord(0.);

	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
	fspace_id = H5Screate_simple(4, dims, NULL);
	#ifndef _ON_FERMI_
	dataset_id = H5Dcreate( file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	#else
	dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	#endif
	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	mspace_id = H5Screate_simple(4, count, NULL);
	status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);
	status = H5Sclose(mspace_id);
	status = H5Sclose(fspace_id);
	status = H5Dclose(dataset_id);
	status = H5Pclose(fapl_id);
	status = H5Fclose(file_id);
	H5close();

	delete [] array_all;

	if (rank==0) {
		char wrapper[256];
		sprintf(wrapper, "%s/%s.xmf", dump_path.c_str(), f_name.c_str());
		FILE *xmf = 0;
		xmf = fopen(wrapper, "w");
		fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
		fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
		fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
		fprintf(xmf, " <Domain>\n");
		fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
		//fprintf(xmf, "     <Time Value=\"%05d\"/>\n", iCounter);
		fprintf(xmf, "     <Time Value=\"%e\"/>\n", absTime);
		fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n", (int)dims[0], (int)dims[1], (int)dims[2]);
		fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
		fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", 0.,0.,0.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", grid.getH(), grid.getH(), grid.getH());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Geometry>\n");
		fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Node\">\n", StreamerHDF5::getAttributeName());
		fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3]);
		fprintf(xmf, "        %s:/data\n",(f_name+".h5").c_str());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Attribute>\n");
		fprintf(xmf, "   </Grid>\n");
		fprintf(xmf, " </Domain>\n");
		fprintf(xmf, "</Xdmf>\n");
		fclose(xmf);
	}
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}

template<typename TGrid, typename Streamer, int channel>
static void DumpHDF5_MPI_Channel(TGrid &grid, const Real absTime,
      const string f_name, const string dump_path="")
{
	#ifdef _USE_HDF_
	typedef typename TGrid::BlockType B;

	int rank;
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;
  MPI_Comm comm = grid.getCartComm();
	MPI_Comm_rank(comm, &rank);

	int coords[3];
	grid.peindex(coords);
	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	float* array_all = new float[NX * NY * NZ];

	hsize_t count[4] = {
		grid.getResidentBlocksPerDimension(2)*B::sizeZ,
		grid.getResidentBlocksPerDimension(1)*B::sizeY,
		grid.getResidentBlocksPerDimension(0)*B::sizeX, 1};

	hsize_t dims[4] = {
		grid.getBlocksPerDimension(2)*B::sizeZ,
		grid.getBlocksPerDimension(1)*B::sizeY,
		grid.getBlocksPerDimension(0)*B::sizeX, 1};

	hsize_t offset[4] = {
		coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
		coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
		coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

	//for (int channel=0; channel<5; channel++)
  {
		char filename[256];
    stringstream ssR;
  	ssR<<"channel"<<channel<<"_"<<f_name;
    string f_name_channel = ssR.str();
		sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name_channel.c_str());

		//printf("Saving channel %d on file %s (%s)\n", channel, f_name_channel.c_str(), filename);
    //fflush(0);

		H5open();
		fapl_id = H5Pcreate(H5P_FILE_ACCESS);
		status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
		file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
		status = H5Pclose(fapl_id);

		#pragma omp parallel for
		for(int i=0; i<vInfo_local.size(); i++) {
			BlockInfo& info = vInfo_local[i];
			const unsigned long idx[3] = {info.index[0], info.index[1], info.index[2]};
			B & b = *(B*)info.ptrBlock;
			Streamer streamer(b);

			for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++) {
				const unsigned long int gx = idx[0]*B::sizeX + ix;
				const unsigned long int gy = idx[1]*B::sizeY + iy;
				const unsigned long int gz = idx[2]*B::sizeZ + iz;
				float * const ptr_output = array_all + (gx + NX * (gy + NY * gz));
				streamer.dump(ix, iy, iz, ptr_output, channel);
			}
		}

		fapl_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

		fspace_id = H5Screate_simple(4, dims, NULL);
		#ifndef _ON_FERMI_
			dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		#else
			dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		#endif

		fspace_id = H5Dget_space(dataset_id);
		H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		mspace_id = H5Screate_simple(4, count, NULL);
		status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

		status = H5Sclose(mspace_id);
		status = H5Sclose(fspace_id);
		status = H5Dclose(dataset_id);
		status = H5Pclose(fapl_id);
		status = H5Fclose(file_id);
		H5close();

		if (rank==0)
		{
			char wrapper[256];
			sprintf(wrapper, "%s/%s.xmf", dump_path.c_str(), f_name_channel.c_str());
			FILE *xmf = 0;
			xmf = fopen(wrapper, "w");
			fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
			fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
			fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
			fprintf(xmf, " <Domain>\n");
			fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
			//fprintf(xmf, "     <Time Value=\"%05d\"/>\n", iCounter);
			fprintf(xmf, "     <Time Value=\"%e\"/>\n", absTime);
			fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n",
																			(int)dims[0], (int)dims[1], (int)dims[2]);
			fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
			fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
			fprintf(xmf, "        %e %e %e\n", 0.,0.,0.);
			fprintf(xmf, "       </DataItem>\n");
			fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
			fprintf(xmf, "        %e %e %e\n", grid.getH(), grid.getH(), grid.getH());
			fprintf(xmf, "       </DataItem>\n");
			fprintf(xmf, "     </Geometry>\n");

			fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"Scalar\" Center=\"Node\">\n");
			fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
												(int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3]);
			fprintf(xmf, "        %s:/data\n",(f_name_channel+".h5").c_str());
			fprintf(xmf, "       </DataItem>\n");
			fprintf(xmf, "     </Attribute>\n");

			fprintf(xmf, "   </Grid>\n");
			fprintf(xmf, " </Domain>\n");
			fprintf(xmf, "</Xdmf>\n");
			fclose(xmf);
		}
	}

	delete [] array_all;
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}

template<typename TGrid, typename Streamer>
static void DumpHDF5_MPI_Vector(TGrid &grid, const Real absTime, const string f_name, const string dump_path="")
{
	#ifdef _USE_HDF_
	typedef typename TGrid::BlockType B;

	int rank;
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;
  MPI_Comm comm = grid.getCartComm();
	MPI_Comm_rank(comm, &rank);

	int coords[3];
	grid.peindex(coords);
	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const unsigned long NCHANNELS = 3;
	float* array_all = new float[NX * NY * NZ * NCHANNELS];

	hsize_t count[4] = {
		grid.getResidentBlocksPerDimension(2)*B::sizeZ,
		grid.getResidentBlocksPerDimension(1)*B::sizeY,
		grid.getResidentBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t dims[4] = {
		grid.getBlocksPerDimension(2)*B::sizeZ,
		grid.getBlocksPerDimension(1)*B::sizeY,
		grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t offset[4] = {
		coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
		coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
		coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

	char filename[256];
	string f_name_channel = f_name;
	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name_channel.c_str());

	//printf("Saving channel %d on file %s (%s)\n", channel, f_name_channel.c_str(), filename);
  //fflush(0);

	H5open();
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	status = H5Pclose(fapl_id);

	#pragma omp parallel for
	for(int i=0; i<vInfo_local.size(); i++) {
		BlockInfo& info = vInfo_local[i];
		const unsigned long idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		//Streamer streamer(b);

		for(int iz=0; iz<B::sizeZ; iz++)
		for(int iy=0; iy<B::sizeY; iy++)
		for(int ix=0; ix<B::sizeX; ix++) {
			const unsigned long int gx = idx[0]*B::sizeX + ix;
			const unsigned long int gy = idx[1]*B::sizeY + iy;
			const unsigned long int gz = idx[2]*B::sizeZ + iz;
      float * const ptr_output = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));
      ptr_output[0] = b(ix, iy, iz).u;
      ptr_output[1] = b(ix, iy, iz).v;
      ptr_output[2] = b(ix, iy, iz).w;
			//streamer.dump(ix, iy, iz, ptr_output, channel);
		}
	}

	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

	fspace_id = H5Screate_simple(4, dims, NULL);
	#ifndef _ON_FERMI_
		dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	#else
		dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	#endif

	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	mspace_id = H5Screate_simple(4, count, NULL);
	status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

	status = H5Sclose(mspace_id);
	status = H5Sclose(fspace_id);
	status = H5Dclose(dataset_id);
	status = H5Pclose(fapl_id);
	status = H5Fclose(file_id);
	H5close();

	if (rank==0)
	{
		char wrapper[256];
		sprintf(wrapper, "%s/%s.xmf", dump_path.c_str(), f_name_channel.c_str());
		FILE *xmf = 0;
		xmf = fopen(wrapper, "w");
		fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
		fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
		fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
		fprintf(xmf, " <Domain>\n");
		fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
		//fprintf(xmf, "     <Time Value=\"%05d\"/>\n", iCounter);
		fprintf(xmf, "     <Time Value=\"%e\"/>\n", absTime);
		fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n",
																		(int)dims[0], (int)dims[1], (int)dims[2]);
		fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
		fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", 0.,0.,0.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", grid.getH(), grid.getH(), grid.getH());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Geometry>\n");

		fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"Vector\" Center=\"Node\">\n");
		fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
											(int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3]);
		fprintf(xmf, "        %s:/data\n",(f_name_channel+".h5").c_str());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Attribute>\n");

		fprintf(xmf, "   </Grid>\n");
		fprintf(xmf, " </Domain>\n");
		fprintf(xmf, "</Xdmf>\n");
		fclose(xmf);
	}

	delete [] array_all;
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}

template<typename TGrid, typename Streamer>
static void ReadHDF5_MPI(TGrid &grid, const string f_name, const string dump_path=".")
{
	#ifdef _USE_HDF_
	typedef typename TGrid::BlockType B;

	int rank;
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

	MPI_Comm comm = grid.getCartComm();
	MPI_Comm_rank(comm, &rank);

	int coords[3];
	grid.peindex(coords);

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const unsigned long int NCHANNELS = Streamer::NCHANNELS;

	float* array_all = new float[NX * NY * NZ * NCHANNELS];

	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	static const unsigned long int sX = 0;
	static const unsigned long int sY = 0;
	static const unsigned long int sZ = 0;

	const unsigned long int eX = B::sizeX;
	const unsigned long int eY = B::sizeY;
	const unsigned long int eZ = B::sizeZ;

	hsize_t count[4] = {
			grid.getResidentBlocksPerDimension(2)*B::sizeZ,
			grid.getResidentBlocksPerDimension(1)*B::sizeY,
			grid.getResidentBlocksPerDimension(0)*B::sizeX, Streamer::NCHANNELS};

	hsize_t dims[4] = {
			grid.getBlocksPerDimension(2)*B::sizeZ,
			grid.getBlocksPerDimension(1)*B::sizeY,
			grid.getBlocksPerDimension(0)*B::sizeX, Streamer::NCHANNELS};

	hsize_t offset[4] = {
			coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
			coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
			coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());
  if(!rank)
    printf("About to restart from %s\n",filename);
  fflush(0);
	H5open();
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
	status = H5Pclose(fapl_id);

	dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

	mspace_id = H5Screate_simple(4, count, NULL);
	status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

	#pragma omp parallel for
	for(int i=0; i<vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const unsigned long int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		Streamer streamer(b);

		for(int iz=sZ; iz<eZ; iz++)
		for(int iy=sY; iy<eY; iy++)
		for(int ix=sX; ix<eX; ix++) {
			const unsigned long int gx = idx[0]*B::sizeX + ix;
			const unsigned long int gy = idx[1]*B::sizeY + iy;
			const unsigned long int gz = idx[2]*B::sizeZ + iz;
			float * const ptr_input = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));
			Real input[NCHANNELS];
			for(int i=0; i<NCHANNELS; ++i) input[i] = (Real)ptr_input[i];
			streamer.operate(input, ix, iy, iz);
		}
	}

	status = H5Pclose(fapl_id);
	status = H5Dclose(dataset_id);
	status = H5Sclose(fspace_id);
	status = H5Sclose(mspace_id);
	status = H5Fclose(file_id);

	H5close();

	delete [] array_all;
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}

template<typename TGrid, typename Streamer, int channel>
static void ReadHDF5_MPI_Channel(TGrid &grid, const string f_name, const string dump_path=".")
{
	#ifdef _USE_HDF_
	typedef typename TGrid::BlockType B;

	int rank;
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

	MPI_Comm comm = grid.getCartComm();
	MPI_Comm_rank(comm, &rank);

	int coords[3];
	grid.peindex(coords);

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const unsigned long int NCHANNELS = 1;

	float* array_all = new float[NX * NY * NZ * NCHANNELS];

	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	hsize_t count[4] = {
			grid.getResidentBlocksPerDimension(2)*B::sizeZ,
			grid.getResidentBlocksPerDimension(1)*B::sizeY,
			grid.getResidentBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t dims[4] = {
			grid.getBlocksPerDimension(2)*B::sizeZ,
			grid.getBlocksPerDimension(1)*B::sizeY,
			grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t offset[4] = {
			coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
			coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
			coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

	//for (int channel=0; channel<5; channel++)
	{
		string f_name_channel = "channel" + to_string(channel) + "_" + f_name;
		printf("About to read from %s\n", f_name_channel.c_str());
		sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name_channel.c_str());
		fapl_id = H5Pcreate(H5P_FILE_ACCESS);
		status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
		file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
		status = H5Pclose(fapl_id);

		dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
		fapl_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

		fspace_id = H5Dget_space(dataset_id);
		H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

		mspace_id = H5Screate_simple(4, count, NULL);
		status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

		#pragma omp parallel for
		for(int i=0; i<vInfo_local.size(); i++)
		{
			BlockInfo& info = vInfo_local[i];
			const unsigned long int idx[3] = {info.index[0], info.index[1], info.index[2]};
			B & b = *(B*)info.ptrBlock;
			Streamer streamer(b);

			for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++) {
				const unsigned long int gx = idx[0]*B::sizeX + ix;
				const unsigned long int gy = idx[1]*B::sizeY + iy;
				const unsigned long int gz = idx[2]*B::sizeZ + iz;
				float * const ptr_input = array_all + (gx + NX * (gy + NY * gz));
				streamer.load<channel>(*ptr_input, ix, iy, iz);
			}
		}

		status = H5Pclose(fapl_id);
		status = H5Dclose(dataset_id);
		status = H5Sclose(fspace_id);
		status = H5Sclose(mspace_id);
		status = H5Fclose(file_id);

		H5close();
	}

	delete [] array_all;
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}

template<typename TGrid, typename Streamer>
static void ReadHDF5_MPI_Vector(TGrid &grid, const string f_name, const string dump_path=".")
{
	#ifdef _USE_HDF_
	typedef typename TGrid::BlockType B;

	int rank;
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

	MPI_Comm comm = grid.getCartComm();
	MPI_Comm_rank(comm, &rank);

	int coords[3];
	grid.peindex(coords);

	const unsigned long int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned long int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned long int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const unsigned long int NCHANNELS = 3;

	float* array_all = new float[NX * NY * NZ * NCHANNELS];

	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	hsize_t count[4] = {
			grid.getResidentBlocksPerDimension(2)*B::sizeZ,
			grid.getResidentBlocksPerDimension(1)*B::sizeY,
			grid.getResidentBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t dims[4] = {
			grid.getBlocksPerDimension(2)*B::sizeZ,
			grid.getBlocksPerDimension(1)*B::sizeY,
			grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

	hsize_t offset[4] = {
			coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
			coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
			coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

	string f_name_channel = f_name;
  if(!rank)
	  printf("About to read from %s\n", f_name_channel.c_str());
	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name_channel.c_str());

	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
	status = H5Pclose(fapl_id);

	dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

	mspace_id = H5Screate_simple(4, count, NULL);
	status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

	#pragma omp parallel for
	for(int i=0; i<vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const unsigned long int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		Streamer streamer(b);

		for(int iz=0; iz<B::sizeZ; iz++)
		for(int iy=0; iy<B::sizeY; iy++)
		for(int ix=0; ix<B::sizeX; ix++) {
			const unsigned long int gx = idx[0]*B::sizeX + ix;
			const unsigned long int gy = idx[1]*B::sizeY + iy;
			const unsigned long int gz = idx[2]*B::sizeZ + iz;
      float * const ptr_input = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));
      b(ix, iy, iz).u = ptr_input[0];
      b(ix, iy, iz).v = ptr_input[1];
      b(ix, iy, iz).w = ptr_input[2];
		}
	}

	status = H5Pclose(fapl_id);
	status = H5Dclose(dataset_id);
	status = H5Sclose(fspace_id);
	status = H5Sclose(mspace_id);
	status = H5Fclose(file_id);

	H5close();

	delete [] array_all;
	#else
	#warning USE OF HDF WAS DISABLED AT COMPILE TIME
	#endif
}
#endif
