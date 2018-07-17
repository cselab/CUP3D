//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Fabian Wermelinger.
//

#ifndef SLICE_H_SJXK5CW4
#define SLICE_H_SJXK5CW4

#include <mpi.h>
#include <vector>
#include <sstream>
#include <string>

#include "HDF5SliceDumperMPI.h"


template <typename TGrid>
struct Slice
{
    typedef TGrid GridType;

    TGrid * grid;
    int id;
    int dir;
    int idx;
    unsigned int width, height;
    bool valid;
    Slice() : grid(NULL), id(-1), dir(-1), idx(-1), width(0), height(0), valid(false) {}

    template <typename TSlice>
    static std::vector<TSlice> getSlices(ArgumentParser& parser, TGrid& grid)
    {
        typedef typename TGrid::BlockType B;
        int Dim[3];
        Dim[0] = grid.getBlocksPerDimension(0)*B::sizeX;
        Dim[1] = grid.getBlocksPerDimension(1)*B::sizeY;
        Dim[2] = grid.getBlocksPerDimension(2)*B::sizeZ;

        std::vector<TSlice> slices(0);
        parser.unset_strict_mode();
        const size_t nSlices = parser("nslices").asInt(1);
        for (size_t i = 0; i < nSlices; ++i)
        {
            TSlice thisOne;
            thisOne.id = i+1;
            thisOne.grid = &grid;
            assert(thisOne.grid != NULL);

            std::ostringstream identifier;
            identifier << "slice" << i+1;
            // fetch direction
            const std::string sDir = identifier.str() + "_direction";
            thisOne.dir = parser(sDir).asInt(2);
            const bool bDirOK = (thisOne.dir >= 0 && thisOne.dir < 3);

            // compute index
            const std::string sFrac  = identifier.str() + "_fraction";
            const double fraction = parser(sFrac).asDouble(0.5);
            const int idx = static_cast<int>(Dim[thisOne.dir] * fraction);
            thisOne.idx = min(Dim[thisOne.dir]-1, idx);
            const bool bIdxOK = (thisOne.idx >= 0 && thisOne.idx < Dim[thisOne.dir]);

            if (bDirOK && bIdxOK) thisOne.valid = true;
            else
            {
                std::cerr << "Slice: WARNING: Ill defined slice \"" << identifier.str() << "\"... Skipping this one" << std::endl;
                thisOne.valid = false;
                slices.push_back(thisOne);
                continue;
            }

            // define slice layout
            if (thisOne.dir == 0)
            {
                thisOne.width  = Dim[2];
                thisOne.height = Dim[1];
            }
            else if (thisOne.dir == 1)
            {
                thisOne.width  = Dim[2];
                thisOne.height = Dim[0];
            }
            else if (thisOne.dir == 2)
            {
                thisOne.width  = Dim[0];
                thisOne.height = Dim[1];
            }
            slices.push_back(thisOne);
        }
        return slices;
    }
};


template <typename TGrid>
struct SliceMPI : public Slice<TGrid>
{
    typedef TGrid GridType;

    unsigned int localWidth, localHeight;
    unsigned int offsetWidth, offsetHeight;
    SliceMPI() : localWidth(-1), localHeight(-1), offsetWidth(-1), offsetHeight(-1) {}

    template <typename TSlice>
    static std::vector<TSlice> getSlices(ArgumentParser& parser, TGrid& grid)
    {
        std::vector<TSlice> slices = Slice<TGrid>::template getSlices<TSlice>(parser, grid);

        typedef typename TGrid::BlockType B;
        int Dim[3];
        Dim[0] = grid.getResidentBlocksPerDimension(0)*B::sizeX;
        Dim[1] = grid.getResidentBlocksPerDimension(1)*B::sizeY;
        Dim[2] = grid.getResidentBlocksPerDimension(2)*B::sizeZ;

        // get slice communicators
        int myRank;
        MPI_Comm_rank(grid.getCartComm(), &myRank);
        int peIdx[3];
        grid.peindex(peIdx);
        int myStart[3], myEnd[3];
        for (int i = 0; i < 3; ++i)
        {
            myStart[i] = Dim[i]*peIdx[i];
            myEnd[i]   = myStart[i] + Dim[i];
        }
        for (size_t i = 0; i < slices.size(); ++i)
        {
            TSlice& s = slices[i];
            const int sIdx = s.idx;
            const int dir  = s.dir;
            if ( !(myStart[dir] <= sIdx && sIdx < myEnd[dir]) )
                s.valid = false;

            // scale index to process local index
            s.idx = s.idx % Dim[s.dir];

            if (s.dir == 0)
            {
                s.localWidth  = Dim[2];
                s.localHeight = Dim[1];
                s.offsetWidth = peIdx[2]*Dim[2];
                s.offsetHeight= peIdx[1]*Dim[1];
            }
            else if (s.dir == 1)
            {
                s.localWidth  = Dim[2];
                s.localHeight = Dim[0];
                s.offsetWidth = peIdx[2]*Dim[2];
                s.offsetHeight= peIdx[0]*Dim[0];
            }
            else if (s.dir == 2)
            {
                s.localWidth  = Dim[0];
                s.localHeight = Dim[1];
                s.offsetWidth = peIdx[0]*Dim[0];
                s.offsetHeight= peIdx[1]*Dim[1];
            }
        }
        return slices;
    }
};

#endif /* SLICE_H_SJXK5CW4 */
