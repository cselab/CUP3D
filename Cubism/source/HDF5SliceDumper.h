//
//  HDF5SliceDumper.h
//  Cubism
//
//  Created by Fabian Wermelinger 09/27/2016
//  Copyright 2016 ETH Zurich. All rights reserved.
//
#ifndef HDF5SLICEDUMPER_H_QI4Y9HO7
#define HDF5SLICEDUMPER_H_QI4Y9HO7

#include <iostream>
#include <vector>
#include <sstream>

#ifdef _USE_HDF_
#include <hdf5.h>
#endif

#ifndef _HDF5_DOUBLE_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
typedef  float hdf5Real;
#else
#define HDF_REAL H5T_NATIVE_DOUBLE
typedef double hdf5Real;
#endif

#include "BlockInfo.h"

// The following requirements for the data Streamer are required:
// Streamer::NCHANNELS        : Number of data elements (1=Scalar, 3=Vector, 9=Tensor)
// Streamer::operate          : Data access methods for read and write
// Streamer::getAttributeName : Attribute name of the date ("Scalar", "Vector", "Tensor")

namespace SliceExtractor
{
    template <typename TBlock, typename TStreamer>
    void YZ(const int ix, const int width, std::vector<BlockInfo>& bInfo, hdf5Real * const data)
    {
        const unsigned int NCHANNELS = TStreamer::NCHANNELS;

#pragma omp parallel for
        for(int i = 0; i < (int)bInfo.size(); ++i)
        {
            BlockInfo& info = bInfo[i];
            const int idx[3] = {info.index[0], info.index[1], info.index[2]};
            TBlock& b = *(TBlock*)info.ptrBlock;

            for(unsigned int iz=0; iz<TBlock::sizeZ; ++iz)
                for(unsigned int iy=0; iy<TBlock::sizeY; ++iy)
                {
                    hdf5Real output[NCHANNELS];
                    for(unsigned int k=0; k<NCHANNELS; ++k)
                        output[k] = 0;

                    TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                    const unsigned int gy = idx[1]*TBlock::sizeY + iy;
                    const unsigned int gz = idx[2]*TBlock::sizeZ + iz;

                    hdf5Real * const ptr = data + NCHANNELS*(gz + width * gy);

                    for(unsigned int k=0; k<NCHANNELS; ++k)
                        ptr[k] = output[k];
                }
        }
    }

    template <typename TBlock, typename TStreamer>
    void XZ(const int iy, const int width, std::vector<BlockInfo>& bInfo, hdf5Real * const data)
    {
        const unsigned int NCHANNELS = TStreamer::NCHANNELS;

#pragma omp parallel for
        for(int i = 0; i < (int)bInfo.size(); ++i)
        {
            BlockInfo& info = bInfo[i];
            const int idx[3] = {info.index[0], info.index[1], info.index[2]};
            TBlock& b = *(TBlock*)info.ptrBlock;

            for(unsigned int iz=0; iz<TBlock::sizeZ; ++iz)
                for(unsigned int ix=0; ix<TBlock::sizeX; ++ix)
                {
                    hdf5Real output[NCHANNELS];
                    for(unsigned int k=0; k<NCHANNELS; ++k)
                        output[k] = 0;

                    TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                    const unsigned int gx = idx[0]*TBlock::sizeX + ix;
                    const unsigned int gz = idx[2]*TBlock::sizeZ + iz;

                    hdf5Real * const ptr = data + NCHANNELS*(gz + width * gx);

                    for(unsigned int k=0; k<NCHANNELS; ++k)
                        ptr[k] = output[k];
                }
        }
    }

    template <typename TBlock, typename TStreamer>
    void YX(const int iz, const int width, std::vector<BlockInfo>& bInfo, hdf5Real * const data)
    {
        const unsigned int NCHANNELS = TStreamer::NCHANNELS;

#pragma omp parallel for
        for(int i = 0; i < (int)bInfo.size(); ++i)
        {
            BlockInfo& info = bInfo[i];
            const int idx[3] = {info.index[0], info.index[1], info.index[2]};
            TBlock& b = *(TBlock*)info.ptrBlock;

            for(unsigned int iy=0; iy<TBlock::sizeY; ++iy)
                for(unsigned int ix=0; ix<TBlock::sizeX; ++ix)
                {
                    hdf5Real output[NCHANNELS];
                    for(unsigned int k=0; k<NCHANNELS; ++k)
                        output[k] = 0;

                    TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                    const unsigned int gx = idx[0]*TBlock::sizeX + ix;
                    const unsigned int gy = idx[1]*TBlock::sizeY + iy;

                    hdf5Real * const ptr = data + NCHANNELS*(gx + width * gy);

                    for(unsigned int k=0; k<NCHANNELS; ++k)
                        ptr[k] = output[k];
                }
        }
    }
}

template<typename TSlice, typename TStreamer>
void DumpSliceHDF5(const TSlice& slice, const int stepID, const Real t, const std::string fname, const std::string dpath=".", const bool bXMF=true)
{
#ifdef _USE_HDF_
    typedef typename TSlice::GridType::BlockType B;
    const typename TSlice::GridType& grid = *(slice.grid);

    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;
    const unsigned int width = slice.width;
    const unsigned int height = slice.height;

    std::cout << "Allocating " << (width * height * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB of HDF5 slice data" << std::endl;;

    hdf5Real * array_all = new hdf5Real[width * height * NCHANNELS];

    std::vector<BlockInfo> bInfo_local = grid.getBlocksInfo();
    std::vector<BlockInfo> bInfo_slice;
    for (size_t i = 0; i < bInfo_local.size(); ++i)
    {
        const int start = bInfo_local[i].index[slice.dir] * _BLOCKSIZE_;
        if (start <= slice.idx && slice.idx < (start+_BLOCKSIZE_))
            bInfo_slice.push_back(bInfo_local[i]);
    }

    // fname is the base filename without file type extension
    std::ostringstream filename;
    filename << dpath << "/" << fname << "_slice" << slice.id;

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    hsize_t count[3] = {height, width, NCHANNELS};
    hsize_t dims[3] = {height, width, NCHANNELS};
    hsize_t offset[3] = {0, 0, 0};

    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate((filename.str()+".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    status = H5Pclose(fapl_id); if(status<0) H5Eprint1(stdout);

    if (0 == slice.dir)
        SliceExtractor::YZ<B,TStreamer>(slice.idx%_BLOCKSIZE_, width, bInfo_slice, array_all);
    else if (1 == slice.dir)
        SliceExtractor::XZ<B,TStreamer>(slice.idx%_BLOCKSIZE_, width, bInfo_slice, array_all);
    else if (2 == slice.dir)
        SliceExtractor::YX<B,TStreamer>(slice.idx%_BLOCKSIZE_, width, bInfo_slice, array_all);

    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    fspace_id = H5Screate_simple(3, dims, NULL);
#ifndef _ON_FERMI_
    dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    fspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    mspace_id = H5Screate_simple(3, count, NULL);
    status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all); if(status<0) H5Eprint1(stdout);

    status = H5Sclose(mspace_id); if(status<0) H5Eprint1(stdout);
    status = H5Sclose(fspace_id); if(status<0) H5Eprint1(stdout);
    status = H5Dclose(dataset_id); if(status<0) H5Eprint1(stdout);
    status = H5Pclose(fapl_id); if(status<0) H5Eprint1(stdout);
    status = H5Fclose(file_id); if(status<0) H5Eprint1(stdout);
    H5close();

    delete [] array_all;

    // writing xmf wrapper
    if (bXMF)
    {
        FILE *xmf = 0;
        xmf = fopen((filename.str()+".xmf").c_str(), "w");
        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
        fprintf(xmf, "     <Time Value=\"%e\"/>\n", t);
        fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"1 %d %d\"/>\n", height, width);
        fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
        fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "        %e %e %e\n", 0., 0., 0.);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "        %e %e %e\n", 1.,1.,1.);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");
        fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Node\">\n", TStreamer::getAttributeName());
        fprintf(xmf, "       <DataItem Dimensions=\"1 %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", height, width, NCHANNELS);
        fprintf(xmf, "        %s:/data\n",(filename.str()+".h5").c_str());
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

#endif /* HDF5SLICEDUMPER_H_QI4Y9HO7 */
