//
//  HDF5SliceDumperMPI.h
//  Cubism
//
//  Created by Fabian Wermelinger 09/28/2016
//  Copyright 2016 ETH Zurich. All rights reserved.
//
#ifndef HDF5SLICEDUMPERMPI_H_ZENQHJA6
#define HDF5SLICEDUMPERMPI_H_ZENQHJA6

#include <mpi.h>
#include "HDF5SliceDumper.h"


template<typename TSlice, typename TStreamer>
void DumpSliceHDF5MPI(const TSlice& slice, const int stepID, const Real t, const std::string fname, const std::string dpath=".", const bool bXMF=true)
{
#ifdef _USE_HDF_
    typedef typename TSlice::GridType::BlockType B;
    const typename TSlice::GridType& grid = *(slice.grid);

    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;
    const unsigned int width = slice.localWidth;
    const unsigned int height = slice.localHeight;

    int sliceRank;
    MPI_Comm comm = grid.getCartComm();
    MPI_Comm_rank(comm, &sliceRank);

// #ifndef NDEBUG
    if (0 == sliceRank)
    {
        std::cout << "Allocating " << (width * height * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB of HDF5 slice data";
    }
// #endif /* NDEBUG */

    hdf5Real * array_all = new hdf5Real[width * height * NCHANNELS];

    std::vector<BlockInfo> bInfo_local = grid.getResidentBlocksInfo();
    std::vector<BlockInfo> bInfo_slice;
    for (size_t i = 0; i < bInfo_local.size(); ++i)
    {
        const int start = bInfo_local[i].index[slice.dir] * _BLOCKSIZE_;
        if (start <= slice.idx && slice.idx < (start+_BLOCKSIZE_))
            bInfo_slice.push_back(bInfo_local[i]);
    }

    ostringstream filename;
    filename << dpath << "/" << fname << TStreamer::postfix() << "_slice" << slice.id;

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    hsize_t count[3] = {height, width, NCHANNELS}; // local
    hsize_t dims[3] = {slice.height, slice.width, NCHANNELS}; // global
    hsize_t offset[3] = {slice.offsetHeight, slice.offsetWidth, 0}; // file offset

// #ifndef NDEBUG
    if (0 == sliceRank)
    {
        std::cout << " (Total  " << (dims[0] * dims[1] * dims[2] * sizeof(hdf5Real))/(1024.*1024.) << " MB)" << std::endl;
    }
// #endif /* NDEBUG */

    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
    file_id = H5Fcreate((filename.str()+".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    status = H5Pclose(fapl_id);

    if (0 == slice.dir)
        SliceExtractor::YZ<B,TStreamer>(slice.idx%_BLOCKSIZE_, width, bInfo_slice, array_all);
    else if (1 == slice.dir)
        SliceExtractor::XZ<B,TStreamer>(slice.idx%_BLOCKSIZE_, width, bInfo_slice, array_all);
    else if (2 == slice.dir)
        SliceExtractor::YX<B,TStreamer>(slice.idx%_BLOCKSIZE_, width, bInfo_slice, array_all);

    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

    fspace_id = H5Screate_simple(3, dims, NULL);
#ifndef _ON_FERMI_
    dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    fspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    mspace_id = H5Screate_simple(3, count, NULL);
    status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

    status = H5Sclose(mspace_id);
    status = H5Sclose(fspace_id);
    status = H5Dclose(dataset_id);
    status = H5Pclose(fapl_id);
    status = H5Fclose(file_id);
    H5close();

    delete [] array_all;

    // writing xmf wrapper
    if (bXMF && 0 == sliceRank)
    {
        FILE *xmf = 0;
        xmf = fopen((filename.str()+".xmf").c_str(), "w");
        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
        fprintf(xmf, "     <Time Value=\"%e\"/>\n", t);
        fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"1 %d %d\"/>\n", slice.height, slice.width);
        fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
        fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "        %e %e %e\n", 0., 0., 0.);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "        %e %e %e\n", 1.,1.,1.);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");
        fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Node\">\n", TStreamer::getAttributeName());
        fprintf(xmf, "       <DataItem Dimensions=\"1 %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", slice.height, slice.width, NCHANNELS);
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

#endif /* HDF5SLICEDUMPERMPI_H_ZENQHJA6 */
