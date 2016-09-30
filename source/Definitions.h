//
//  DataStructures.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_DataStructures_h
#define CubismUP_3D_DataStructures_h

//#include <cassert>

//#define BBURST 1
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <math.h>
#include <string>
#include <vector>
#include <array>
using namespace std;
//#include <assert.h>

// utmost import to be defined before including cubism

#ifndef _SP_COMP_
typedef double Real;
#else // _SP_COMP_
typedef float Real;
#endif // _SP_COMP_

#include <mpi.h>
#include <omp.h>

//this is all cubism file we need
#include <ArgumentParser.h>
#include <Grid.h>
#include <GridMPI.h>
#include <BlockInfo.h>
#ifdef _VTK_
#include <SerializerIO_ImageVTK.h>
#endif
//#include <HDF5Dumper_MPI.h>
//#include <ZBinDumper_MPI.h>
#include <BlockLab.h>
#include <BlockLabMPI.h>
#include <Profiler.h>
#include <StencilInfo.h>
#include "Timer.h"
#include "BoundaryConditions.h"

#include "Communicator.h"
/*
#ifndef _BS_
#define _BS_ 32
#endif // _BS_

#ifndef _BSX_
#define _BSX_ 32
#endif // _BSX_

#ifndef _BSY_
#define _BSY_ 32
#endif // _BSY_

#ifndef _BSZ_
#define _BSZ_ 32
#endif // _BSZ_
*/

struct FluidElement
{
    Real u, v, w, chi, p, tmpU, tmpV, tmpW;
	
    FluidElement() : u(0), v(0), w(0), chi(0), p(0), tmpU(0), tmpV(0), tmpW(0)
	{}
    
    void clear()
    {
        u = v = w = chi = p = tmpU = tmpV = tmpW = 0;
    }
};

struct FluidVTKStreamer
{
	static const int channels = 8;
	
	void operate(FluidElement input, Real output[channels])
	{
		output[0] = input.u;
		output[1] = input.v;
		output[2] = input.w;
		output[3] = input.chi;
		output[4] = input.p;
		output[5] = input.tmpU;
		output[6] = input.tmpV;
		output[7] = input.tmpW;
	}
};

// this is used for serialization - important that ALL the quantities are streamed
struct StreamerGridPoint
{
	static const int channels = 8;
	
	void operate(const FluidElement& input, Real output[channels]) const
	{
		abort();
		output[0] = input.u;
		output[1] = input.v;
		output[2] = input.w;
		output[3] = input.chi;
		output[4] = input.p;
		output[5] = input.tmpU;
		output[6] = input.tmpV;
		output[7] = input.tmpW;
	}
	
	void operate(const Real input[channels], FluidElement& output) const
	{
		abort();
		output.u    = input[0];
		output.v    = input[1];
		output.w    = input[2];
		output.chi  = input[3];
		output.p    = input[4];
		output.tmpU = input[5];
		output.tmpV = input[6];
		output.tmpW = input[7];
	}
};

struct StreamerGridPointASCII
{ //TODO: how to terminate operate output?? endl? " "?
	void operate(const FluidElement& input, ofstream& output) const
	{
		output << input.u << " " << input.v << " " << input.w << " " << input.chi << " "
			   << input.p << " " << input.tmpU << " " << input.tmpV << " " << input.tmpW;
	}
	
	void operate(ifstream& input, FluidElement& output) const
	{
		input >> output.u;
		input >> output.v;
		input >> output.w;
		input >> output.chi;
		input >> output.p;
		input >> output.tmpU;
		input >> output.tmpV;
		input >> output.tmpW;
	}
};

struct StreamerDiv
{
	static const int channels = 1;
	static void operate(const FluidElement& input, Real output[1])
	{
		output[0] = input.p;
    }
    
    static void operate(const Real input[1], FluidElement& output)
    {
        output.p = input[0];
    }
};

struct ObstacleBlock
{
    static const int sizeX = _BSX_;
    static const int sizeY = _BSY_;
    static const int sizeZ = _BSZ_;
    Real chi[sizeX][sizeY][sizeZ];
    Real udef[sizeX][sizeY][sizeZ][3];

    void clear()
    {
        memset(chi, 0, sizeof(Real)*sizeX*sizeY*sizeZ);
        memset(udef, 0, sizeof(Real)*sizeX*sizeY*sizeZ*3);
    }
};

struct FluidBlock
{
    //these identifiers are required by cubism!
    static const int sizeX = _BSX_;
    static const int sizeY = _BSY_;
    static const int sizeZ = _BSZ_;
	
	typedef FluidElement ElementType;
	typedef FluidElement element_type;
	
	FluidElement data[sizeZ][sizeY][sizeX];
	
    //required from Grid.h
    void clear()
    {
        FluidElement * entry = &data[0][0][0];
        const int N = sizeX*sizeY*sizeZ;
        
        for(int i=0; i<N; ++i)
            entry[i].clear();
    }
    
    FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
		assert(ix>=0); assert(ix<sizeX);
		assert(iy>=0); assert(iy<sizeY);
		assert(iz>=0); assert(iz<sizeZ);
		
        return data[iz][iy][ix];
	}
	
	template <typename Streamer>
	inline void Write(ofstream& output, Streamer streamer) const
	{
		for(int iz=0; iz<sizeZ; iz++)
			for(int iy=0; iy<sizeY; iy++)
				for(int ix=0; ix<sizeX; ix++)
					streamer.operate(data[iz][iy][ix], output);
	}
	
	template <typename Streamer>
	inline void Read(ifstream& input, Streamer streamer)
	{
		for(int iz=0; iz<sizeZ; iz++)
			for(int iy=0; iy<sizeY; iy++)
				for(int ix=0; ix<sizeX; ix++)
					streamer.operate(input, data[iz][iy][ix]);
	}
};

struct surfData
{
    int blockID, ix, iy, iz;
    Real dchidx, dchidy, dchidz, delta;

    surfData(const int _blockID, const int _ix, const int _iy, const int _iz, const Real _dchidx, const Real _dchidy, const Real _dchidz, const Real _delta)
    : blockID(_blockID), ix(_ix), iy(_iy), iz(_iz), dchidx(_dchidx), dchidy(_dchidy), dchidz(_dchidz), delta(_delta)
    {}

    void set(const int _blockID, const int _ix, const int _iy, const int _iz, const Real _dchidx, const Real _dchidy, const Real _dchidz, const Real _delta)
    {
        blockID=_blockID; ix=_ix; iy=_iy; iz=_iz; dchidx=_dchidx; dchidy=_dchidy; dchidz=_dchidz; delta=_delta;
    }
};

struct surfaceBlocks
{
public:
	int Ndata;
    vector<surfData*> Set;

	surfaceBlocks() : Ndata(0) {}

	~surfaceBlocks()
	{
		for (auto & trash : Set) {
			if(trash == nullptr) continue;
			delete trash;
			trash = nullptr;
		}
	}

    void add(const int blockID, const int ix, const int iy, const int iz, const Real dchidx, const Real dchidy, const Real dchidz, const Real delta)
    {
        if(Ndata+1>Set.size()) {
        		surfData * tmp = new surfData(blockID,ix,iy,iz,dchidx,dchidy,dchidz,delta);
            Set.push_back(tmp);
        } else Set[Ndata]->set(blockID,ix,iy,iz,dchidx,dchidy,dchidz,delta);
        Ndata++;
    }
};

struct surfacePoints
{
public:
    Real *pX, *pY, *pZ, *P, *fX, *fY, *fZ, *fxP, *fyP, *fzP, *fxV, *fyV, *fzV, *vx, *vy, *vz, *vxDef, *vyDef, *vzDef;
    int Ndata, nAlloc, nMapped, *gridMap;
    vector<surfData*> Set;

    surfacePoints() :
    	Ndata(0), nAlloc(0), nMapped(0), pX(nullptr), pY(nullptr), pZ(nullptr), P(nullptr),
    fX(nullptr), fY(nullptr), fZ(nullptr), fxP(nullptr), fyP(nullptr), fzP(nullptr), fxV(nullptr), fyV(nullptr), fzV(nullptr),
    vx(nullptr), vy(nullptr), vz(nullptr), vxDef(nullptr), vyDef(nullptr), vzDef(nullptr), gridMap(nullptr)
    { }

    ~surfacePoints()
    {
        if(pX      not_eq nullptr){delete[] pX;      pX=nullptr;     }
        if(pY      not_eq nullptr){delete[] pY;      pY=nullptr;     }
        if(pZ      not_eq nullptr){delete[] pZ;      pZ=nullptr;     }
        if(P       not_eq nullptr){delete[] P;   	 P=nullptr;  	}
        if(fX      not_eq nullptr){delete[] fX;      fX=nullptr;     }
        if(fY      not_eq nullptr){delete[] fY;      fY=nullptr;     }
        if(fZ      not_eq nullptr){delete[] fZ;      fZ=nullptr;     }
        if(fxP     not_eq nullptr){delete[] fxP;     fxP=nullptr;    }
        if(fyP     not_eq nullptr){delete[] fyP;     fyP=nullptr;    }
        if(fzP     not_eq nullptr){delete[] fzP;     fzP=nullptr;    }
        if(fxV     not_eq nullptr){delete[] fxV;     fxV=nullptr;    }
        if(fyV     not_eq nullptr){delete[] fyV;     fyV=nullptr;    }
        if(fzV     not_eq nullptr){delete[] fzV;     fzV=nullptr;    }
        if(vx      not_eq nullptr){delete[] vx;      vx=nullptr;     }
        if(vy      not_eq nullptr){delete[] vy;      vy=nullptr;     }
        if(vz      not_eq nullptr){delete[] vz;      vz=nullptr;     }
        if(vxDef   not_eq nullptr){delete[] vxDef;   vxDef=nullptr;  }
        if(vyDef   not_eq nullptr){delete[] vyDef;   vyDef=nullptr;  }
        if(vzDef   not_eq nullptr){delete[] vzDef;   vzDef=nullptr;  }
        if(gridMap not_eq nullptr){delete[] gridMap; gridMap=nullptr;}
    }

    void _add(const surfData* c)
    {
        if(Ndata+1>Set.size()) {
        		surfData * tmp = new surfData(c->blockID,c->ix,c->iy,c->iz,c->dchidx,c->dchidy,c->dchidz,c->delta);
            Set.push_back(tmp);
        } else Set[Ndata]->set(c->blockID,c->ix,c->iy,c->iz,c->dchidx,c->dchidy,c->dchidz,c->delta);
        Ndata++;
    }

    void finalizeOnGrid(vector<surfaceBlocks>& blocksPerThread)
    {
		Ndata = 0;
		for(int i=0; i<blocksPerThread.size(); i++) {
			//printf("Processing chunk %d of size %d, current size is %d (%d)\n", i,blocksPerThread[i].Set.size(),Ndata,Set.size());
			for(int j=0; j<blocksPerThread[i].Set.size(); j++)
				_add(blocksPerThread[i].Set[j]);


			//printf("Processed chunk %d of size %d, current size is %d (%d)\n", i,blocksPerThread[i].Set.size(),Ndata,Set.size());
		}
			//for (auto & elem : blocksPerThread[i].Set) _add(elem);

        if (Ndata > nAlloc) {
            nAlloc = Ndata;
            if(pX      not_eq nullptr){delete[] pX;      pX=nullptr;     }
            if(pY      not_eq nullptr){delete[] pY;      pY=nullptr;     }
            if(pZ      not_eq nullptr){delete[] pZ;      pZ=nullptr;     }
            if(fX      not_eq nullptr){delete[] fX;      fX=nullptr;     }
            if(fY      not_eq nullptr){delete[] fY;      fY=nullptr;     }
            if(fZ      not_eq nullptr){delete[] fZ;      fZ=nullptr;     }
            if(fxP     not_eq nullptr){delete[] fxP;     fxP=nullptr;    }
            if(fyP     not_eq nullptr){delete[] fyP;     fyP=nullptr;    }
            if(fzP     not_eq nullptr){delete[] fzP;     fzP=nullptr;    }
            if(fxV     not_eq nullptr){delete[] fxV;     fxV=nullptr;    }
            if(fyV     not_eq nullptr){delete[] fyV;     fyV=nullptr;    }
            if(fzV     not_eq nullptr){delete[] fzV;     fzV=nullptr;    }
            if(vx      not_eq nullptr){delete[] vx;      vx=nullptr;     }
            if(vy      not_eq nullptr){delete[] vy;      vy=nullptr;     }
            if(vz      not_eq nullptr){delete[] vz;      vz=nullptr;     }
            if(vxDef   not_eq nullptr){delete[] vxDef;   vxDef=nullptr;  }
            if(vyDef   not_eq nullptr){delete[] vyDef;   vyDef=nullptr;  }
            if(vzDef   not_eq nullptr){delete[] vzDef;   vzDef=nullptr;  }
            if(P       not_eq nullptr){delete[] P;   	P =nullptr;  	}
            if(gridMap not_eq nullptr){delete[] gridMap; gridMap=nullptr;}

            pX      = new Real[nAlloc]; pY      = new Real[nAlloc]; pZ      = new Real[nAlloc];
            fX      = new Real[nAlloc]; fY      = new Real[nAlloc]; fZ      = new Real[nAlloc];
            fxP     = new Real[nAlloc]; fyP     = new Real[nAlloc]; fzP     = new Real[nAlloc];
            fxV     = new Real[nAlloc]; fyV     = new Real[nAlloc]; fzV     = new Real[nAlloc];
            vx      = new Real[nAlloc]; vy      = new Real[nAlloc]; vz      = new Real[nAlloc];
            vxDef   = new Real[nAlloc]; vyDef   = new Real[nAlloc]; vzDef   = new Real[nAlloc];
            P       = new Real[nAlloc]; gridMap = new  int[nAlloc];
        }

#ifndef NDEBUG
		int checksum = 0;
		for(int i=0; i<blocksPerThread.size(); i++) checksum += blocksPerThread[i].Ndata;
		assert(checksum==Ndata);
		printf("Random assortment of numbers in the surface blocks stuff: %d == %d <= %d\n",checksum,Ndata, nAlloc);
#endif
    }

    void sort(const Real* const skinX, const Real* const skinY, const Real* const skinZ, const int skinN)
    { //this works if N > skinN
        nMapped = skinN;
        #pragma omp parallel for
        for(int k = 0; k < Ndata; k++) {
            Real dist = 1;
            for(int j=0; j<skinN; j++) {
                const Real _d = std::pow(skinX[j] - pX[k],2)
                               +std::pow(skinY[j] - pY[k],2)
                			   +std::pow(skinZ[j] - pZ[k],2);
                if (_d < dist) {
                    gridMap[k] = j;
                    dist = _d;
                }
            }
        }
    }

    void print(const int ID, const int stepNumber, const int rank)
    {
        ofstream fileskin;
        char buf[500];
        sprintf(buf, "skinDistrib_%1d_%07d_rank%02d.txt", ID, stepNumber, rank);
        string filename(buf);
        fileskin.open(filename, ios::trunc);

        for(int i=0; i<Ndata; ++i) {
            fileskin<<pX[i]<<" "<<pY[i]<<" "<<pZ[i]<<" "<<P[i]<<" "<<
					  fxP[i]<<" "<<fyP[i]<<" "<<fzP[i]<<" "<<fxV[i]<<" "<<fyV[i]<<" "<<fzV[i]<<" "<<
					  vx[i]<<" "<<vy[i]<<" "<<vz[i]<<" "<<vxDef[i]<<" "<<vyDef[i]<<" "<<vzDef[i]<<endl;
        }
        fileskin.close();
    }

    void printSorted(const int ID, const int stepNumber)
    {
        ofstream fileskin;
        char buf[500];
        sprintf(buf, "skinDistribSorted_%1d_%07d.txt", ID, stepNumber);
        string filename(buf);
        fileskin.open(filename, ios::trunc);

        for (int j=0; j<nMapped; j++) {
            Real Fx(0.), Fy(0.), Fz(0.), FxP(0.), FyP(0.), FzP(0.), FxV(0.), FyV(0.), FzV(0.), p(0.);
            Real Vx(0.), Vy(0.), Vz(0.), VxDef(0.), VyDef(0.), VzDef(0.), Px(0.), Py(0.), Pz(0.);
            int cnt(0);

            for(int k = 0; k < Ndata; k++) {
                if (j not_eq gridMap[k]) continue;

                cnt++;
                 VxDef += vxDef[k]; VyDef += vyDef[k]; VzDef += vzDef[k];
                 Vx += vx[k];    Vy += vy[k];    Vz += vz[k];
                 Fx += fX[k];    Fy += fY[k];    Fz += fZ[k];
                FxP += fxP[k];  FyP += fyP[k];  FzP += fzP[k];
                FxV += fxV[k];  FyV += fxV[k];  FzV += fzV[k];
                 Px += pX[k];    Py += pY[k];    Pz += pZ[k]; p += P[k];
            }
            //TMP
            const Real ds = .1*M_PI/(Real)Ndata;
            Fx  /= ds; Fy  /= ds; Fz  /= ds;
            FxP /= ds; FyP /= ds; FzP /= ds;
            FxV /= ds; FyV /= ds; FzV /= ds;
            Vx /= (Real)cnt;
            Vy /= (Real)cnt;
            Vz /= (Real)cnt;
            Px /= (Real)cnt;
            Py /= (Real)cnt;
            Pz /= (Real)cnt;
            p  /= (Real)cnt;

            fileskin<<Px<<" "<<Py<<" "<<Pz<<" "<<p<<" "<<Fx<<" "<<Fy<<" "<<Fz<<" "
                    <<FxP<<" "<<FyP<<" "<<FzP<<" "<<FxV<<" "<<FyV<<" "<<FzV<<" "
                    <<Vx<<" "<<Vy<<" "<<Vz<<" "<<VxDef<<" "<<VyDef<<" "<<VzDef<<" "<<endl;
        }
        fileskin.close();
    }
};

template <> inline void FluidBlock::Write<StreamerGridPoint>(ofstream& output, StreamerGridPoint streamer) const
{
	output.write((const char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}

template <> inline void FluidBlock::Read<StreamerGridPoint>(ifstream& input, StreamerGridPoint streamer)
{
	input.read((char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}

// VP Streamers
struct FluidVPStreamer
{
	static const int channels = 8;
	
	FluidBlock * ref;
	FluidVPStreamer(FluidBlock& b): ref(&b) {}
	FluidVPStreamer(): ref(NULL) {}
	
	template<int channel>
	static inline Real operate(const FluidElement& input) { abort(); return 0; }
	
	inline Real operate(const int ix, const int iy, const int iz) const
	{
		cout << "You must not call this operate method of FluidVPStreamer" << endl;
		abort();
		return 0;
	}
	
	const char * name() { return "FluidVPStreamer" ; }
};

template<> inline Real FluidVPStreamer::operate<0>(const FluidElement& e) { return e.u; }
template<> inline Real FluidVPStreamer::operate<1>(const FluidElement& e) { return e.v; }
template<> inline Real FluidVPStreamer::operate<2>(const FluidElement& e) { return e.w; }
template<> inline Real FluidVPStreamer::operate<3>(const FluidElement& e) { return e.chi; }
template<> inline Real FluidVPStreamer::operate<4>(const FluidElement& e) { return e.p; }
template<> inline Real FluidVPStreamer::operate<5>(const FluidElement& e) { return e.tmpU; }
template<> inline Real FluidVPStreamer::operate<6>(const FluidElement& e) { return e.tmpV; }
template<> inline Real FluidVPStreamer::operate<7>(const FluidElement& e) { return e.tmpW; }

/*
struct TmpVPStreamer
{
	static const int channels = 1;
	
	FluidBlock * ref;
	TmpVPStreamer(FluidBlock& b): ref(&b) {}
	TmpVPStreamer(): ref(NULL) {  }
	
	template<int channel>
	static inline Real operate(const FluidElement& input) { abort(); return 0; }
	
	inline Real operate(const int ix, const int iy, const int iz) const
	{
		cout << "You must not call this operate method of TmpVPStreamer" << endl;
		abort();
		return 0;
	}
	
	const char * name() { return "TmpVPStreamer" ; }
};

template<> inline Real TmpVPStreamer::operate<0>(const FluidElement& e) { return e.tmp; }
*/

struct ScalarStreamer
{
	static const int channels = 1;
	
	void operate(Real input, Real output[1])
	{
		output[0] = input;
	}
};

struct ScalarBlock
{
	static const int sizeX = _BS_/2;
	static const int sizeY = _BS_/2;
	static const int sizeZ = _BS_/2;
	typedef Real ElementType;
	
	Real  data[sizeZ][sizeY][sizeX];
	Real  tmp[sizeZ][sizeY][sizeX];
	
	void clear_data()
	{
		const int N = sizeX*sizeY*sizeZ;
		Real * const e = &data[0][0][0];
		for(int i=0; i<N; ++i) e[i] = 0;
	}
	
	void clear_tmp()
	{
		const int N = sizeX*sizeY*sizeZ;
		Real * const e = &tmp[0][0][0];
		for(int i=0; i<N; ++i) e[i] = 0;
	}
	
	void clear()
	{
		clear_data();
		clear_tmp();
	}
	
	inline Real& operator()(int ix, int iy=0, int iz=0)
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);
		
		return data[iz][iy][ix];
	}
	
	template <typename Streamer>
	inline void Write(ofstream& output, Streamer streamer) const
	{
		for(int iz=0; iz<sizeZ; iz++)
		for(int iy=0; iy<sizeY; iy++)
		for(int ix=0; ix<sizeX; ix++)
		streamer.operate(data[iz][iy][ix], output);
	}
	
	template <typename Streamer>
	inline void Read(ifstream& input, Streamer streamer)
	{
		for(int iz=0; iz<sizeZ; iz++)
		for(int iy=0; iy<sizeY; iy++)
		for(int ix=0; ix<sizeX; ix++)
		streamer.operate(input, data[iz][iy][ix]);
	}
};

template <> inline void ScalarBlock::Write<StreamerGridPoint>(ofstream& output, StreamerGridPoint streamer) const
{
	output.write((const char *)&data[0][0][0], sizeof(Real)*sizeX*sizeY*sizeZ);
}

template <> inline void ScalarBlock::Read<StreamerGridPoint>(ifstream& input, StreamerGridPoint streamer)
{
	input.read((char *)&data[0][0][0], sizeof(Real)*sizeX*sizeY*sizeZ);
}


struct StreamerSerialization
{ //TODO: why save also the tmp fields?? why chi field? (would reduce save time&size by ~50%)
	//static const int NCHANNELS = 10;
	static const int NCHANNELS = 5;
	FluidBlock& ref;
	
	StreamerSerialization(FluidBlock& b): ref(b) {}
	
	void operate(const int ix, const int iy, const int iz, Real output[NCHANNELS]) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		output[0]  = input.u;
		output[1]  = input.v;
		output[2]  = input.w;
		output[3]  = input.chi;
		output[4]  = input.p;
		output[5]  = input.tmpU;
		output[6]  = input.tmpV;
		output[7]  = input.tmpW;
	}
	
	void operate(const Real input[NCHANNELS], const int ix, const int iy, const int iz) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		output.u    = input[0];
		output.v    = input[1];
		output.w    = input[2];
		output.chi  = input[3];
		output.p    = input[4];
		output.tmpU = input[5];
		output.tmpV = input[6];
		output.tmpW = input[7];
	}
	
	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0: *ovalue  = input.u; break;
			case 1: *ovalue  = input.v; break;
			case 2: *ovalue  = input.w; break;
			case 3: *ovalue  = input.chi; break;
			case 4: *ovalue  = input.p; break;
			case 5: *ovalue  = input.tmpU; break;
			case 6: *ovalue  = input.tmpV; break;
			case 7: *ovalue  = input.tmpW; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0:  output.u    = ivalue; break;
			case 1:  output.v    = ivalue; break;
			case 2:  output.w    = ivalue; break;
			case 3:  output.chi  = ivalue; break;
			case 4:  output.p    = ivalue; break;
			case 5:  output.tmpU = ivalue; break;
			case 6:  output.tmpV = ivalue; break;
			case 7:  output.tmpW = ivalue; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	static const char * getAttributeName() { return "Tensor"; }
};

struct StreamerHDF5
{
	static const int NCHANNELS = 9;
	
	FluidBlock& ref;
	
	StreamerHDF5(FluidBlock& b): ref(b) {}
	
	void operate(const int ix, const int iy, const int iz, Real output[NCHANNELS]) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		output[0] = input.u;
		output[1] = input.v;
		output[2] = input.w;
		output[3] = input.chi;
		output[4] = input.p;
		output[5] = input.tmpU;
		output[6] = input.tmpV;
		output[7] = input.tmpW;
		output[8] = 0;
	}
	
	void operate(const Real input[NCHANNELS], const int ix, const int iy, const int iz) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		output.u    = input[0];
		output.v    = input[1];
		output.w    = input[2];
		output.chi  = input[3];
		output.p    = input[4];
		output.tmpU = input[5];
		output.tmpV = input[6];
		output.tmpW = input[7];
	}
	
	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0: *ovalue  = input.u; break;
			case 1: *ovalue  = input.v; break;
			case 2: *ovalue  = input.w; break;
			case 3: *ovalue  = input.chi; break;
			case 4: *ovalue  = input.p; break;
			case 5: *ovalue  = input.tmpU; break;
			case 6: *ovalue  = input.tmpV; break;
			case 7: *ovalue  = input.tmpW; break;
			case 8: *ovalue  = 0; 		   break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0:  output.u    = ivalue; break;
			case 1:  output.v    = ivalue; break;
			case 2:  output.w    = ivalue; break;
			case 3:  output.chi  = ivalue; break;
			case 4:  output.p    = ivalue; break;
			case 5:  output.tmpU = ivalue; break;
			case 6:  output.tmpV = ivalue; break;
			case 7:  output.tmpW = ivalue; break;
			case 8:  	 				 ; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	static const char * getAttributeName() { return "Tensor"; }
};

struct StreamerScalarHDF5
{
	static const int NCHANNELS = 1;
	
	ScalarBlock& ref;
	
	StreamerScalarHDF5(ScalarBlock& b): ref(b) {}
	
	void operate(const int ix, const int iy, const int iz, Real output[NCHANNELS]) const
	{
		const Real& input = ref.data[iz][iy][ix];
		output[0] = input;
	}
	
	void operate(const Real input[NCHANNELS], const int ix, const int iy, const int iz) const
	{
		Real& output = ref.data[iz][iy][ix];
		output  = input[0];
	}
	
	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const Real& input = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0: *ovalue = input; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		Real& output = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0:  output = ivalue; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	static const char * getAttributeName() { return "Scalar"; }
};

/*
template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBottomWall : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
    ElementTypeBlock pDirichlet;
    
	BlockLabBottomWall(): BlockLab<BlockType,allocator>()
    {
        pDirichlet.chi = 0;
		pDirichlet.u = 0;
		pDirichlet.v = 0;
		pDirichlet.w = 0;
        pDirichlet.p = 0;
		pDirichlet.tmpU = 0;
		pDirichlet.tmpV = 0;
		pDirichlet.tmpW = 0;
    }
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		// keep periodicity in x,z direction
		if (info.index[1]==0)		   bc.template applyBC_mixedBottom<1,0>(pDirichlet);
		if (info.index[1]==this->NY-1) bc.template applyBC_mixedTop<1,1>(pDirichlet);
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabPipe : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabPipe(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		if (info.index[1]==0)		   bc.template applyBC_mixedBottom<1,0>();
		if (info.index[1]==this->NY-1) bc.template applyBC_mixedBottom<1,1>();
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabVortex : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabVortex(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		if (info.index[0]==0)		   bc.template applyBC_vortex<0,0>(info);
		if (info.index[0]==this->NX-1) bc.template applyBC_vortex<0,1>(info);
		if (info.index[1]==0)		   bc.template applyBC_vortex<1,0>(info);
		if (info.index[1]==this->NY-1) bc.template applyBC_vortex<1,1>(info);
		if (info.index[2]==0)		   bc.template applyBC_vortex<2,0>(info);
		if (info.index[2]==this->NZ-1) bc.template applyBC_vortex<2,1>(info);
	}
};
*/
typedef Grid<FluidBlock, std::allocator> FluidGrid;
typedef Grid<ScalarBlock, std::allocator> ScalarGrid;
typedef GridMPI<FluidGrid> FluidGridMPI;

//#ifdef _MIXED_
//typedef BlockLabBottomWall<FluidBlock, std::allocator> Lab;
//#endif // _MIXED_

//#ifdef _PERIODIC_
//typedef  Lab;
//#endif // _PERIODIC_

//#ifdef _VORTEX_
//typedef BlockLabVortex<FluidBlock, std::allocator> Lab;
//#endif // _VORTEX_

//#ifdef _PIPE_
//typedef BlockLabPipe<FluidBlock, std::allocator> Lab;
//#endif // _PIPE_

typedef BlockLabMPI<BlockLab<FluidBlock, std::allocator>> LabMPI;


struct Layer
{
	const int sizeX;
	const int sizeY;
	const int sizeZ;
	const int nDim;

	Real * data;

	Layer(const int sizeX, const int sizeY, const int sizeZ, const int nDim) : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), nDim(nDim)
	{
		data = new Real[nDim*sizeX*sizeY*sizeZ];
	}

	~Layer()
	{
		delete [] data;
	}

	inline Real& operator()(int ix=0, int iy=0, int iz=0, int dim=0)
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);

		return data[dim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
	}

	inline Real read(int ix=0, int iy=0, int iz=0, int dim=0) const
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);

		return data[dim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
	}

	const Layer& operator=(const Real val)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iz = 0; iz<sizeZ; iz++)
				for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] = val;

		return *this;
	}

	const Layer& operator=(const Layer& l)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iz = 0; iz<sizeZ; iz++)
			for(int iy = 0; iy<sizeY; iy++)
				for(int ix = 0; ix<sizeX; ix++)
					data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] = l.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];

		return *this;
	}

	template<int dim>
	void clear(Real val)
	{
		for(int iz = 0; iz<sizeZ; iz++)
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[dim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] = val;
	}

	Real getH0() const
	{
		return 1./(Real)sizeX;
	}

	Real getH1() const
	{
		return 1./(Real)sizeY;
	}

	Real getH2() const
	{
		return 1./(Real)sizeZ;
	}

	const vector<Real> operator -(const Layer& layer)
	{
		vector<Real> result;

		//compute linf distance
		{
			Real LInf_diff = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iz = 0; iz<sizeZ; iz++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						LInf_diff = max(LInf_diff, (Real)fabs(data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] - layer.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix]));

			result.push_back(LInf_diff);
		}

		//compute linf distance
		{
			Real L2error = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iz = 0; iz<sizeZ; iz++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						L2error += pow(data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] - layer.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix], 2);

			result.push_back(sqrt((Real)L2error/(sizeY*sizeX)));
		}

		return result;
	}


	void difference(const Layer& input, Layer& difference)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iz = 0; iz<sizeZ; iz++)
			for(int iy = 0; iy<sizeY; iy++)
				for(int ix = 0; ix<sizeX; ix++)
					difference.data[idim*sizeX*sizeY + iy*sizeX + ix]= data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix] - input.data[idim*sizeX*sizeY*sizeZ + iz*sizeX*sizeY + iy*sizeX + ix];
	}

	template<int iDim>
	Real * getPlane()
	{
		return (Real*)&data[iDim*sizeZ*sizeX*sizeY];
	}
};


#include <xmmintrin.h>
template <typename LayerT>
LayerT * allocate()
{
	void * data = _mm_malloc(sizeof(LayerT), 64);
	return new (data) LayerT;
};

#endif
