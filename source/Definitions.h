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
//#define __2Leads_
//#define __DumpWakeStefan 9
#define __useSkin_
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
#define __NpLatLine 20
#define __ExploreHalfWake
#include <mpi.h>
#include <omp.h>

#ifndef _SP_COMP_
typedef double Real;
#else // _SP_COMP_
typedef float Real;
#endif // _SP_COMP_

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

struct FluidElement
{
    Real u, v, w, chi, p, tmpU, tmpV, tmpW;
    FluidElement(): u(0),v(0),w(0),chi(0),p(0),tmpU(0),tmpV(0),tmpW(0) {}
    void clear() { u = v = w = chi = p = tmpU = tmpV = tmpW = 0; }
};

struct DumpElement
{
    Real u, v, w, chi;
    DumpElement() : u(0), v(0), w(0), chi(0) {}
    void clear() { u = v = w = chi = 0; }
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
    printf("Abort in StreamerGridPoint in!\n");
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
    printf("Abort in StreamerGridPoint out!\n");
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
	{ output[0] = input.p; }

  static void operate(const Real input[1], FluidElement& output)
  { output.p = input[0]; }
};

struct ObstacleBlock
{
    static const int sizeX = _BSX_;
    static const int sizeY = _BSY_;
    static const int sizeZ = _BSZ_;
    Real chi[sizeX][sizeY][sizeZ];
    Real udef[sizeX][sizeY][sizeZ][3];
    int sectionMarker[sizeX][sizeY][sizeZ];
    int hinge2Index;
    double hinge2LabFrame[3] = {0.0};

    void clear()
    {
        memset(chi, 0, sizeof(Real)*sizeX*sizeY*sizeZ);
        memset(udef, 0, sizeof(Real)*sizeX*sizeY*sizeZ*3);
	memset(sectionMarker, 0, sizeof(int)*sizeX*sizeY*sizeZ);
	hinge2Index = -100;
    }

    void saveHinge2Loc(double hinge2Location[3])
    {
	    hinge2LabFrame[0] = hinge2Location[0];
	    hinge2LabFrame[1] = hinge2Location[1];
	    hinge2LabFrame[2] = hinge2Location[2];
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

struct DumpBlock
{
  //these identifiers are required by cubism!
  static const int sizeX = _BSX_;
  static const int sizeY = _BSY_;
  static const int sizeZ = _BSZ_;
	typedef DumpElement ElementType;
	typedef DumpElement element_type;
	DumpElement data[sizeZ][sizeY][sizeX];
  //required from Grid.h
  void clear()
  {
    DumpElement * entry = &data[0][0][0];
    const int N = sizeX*sizeY*sizeZ;
    for(int i=0; i<N; ++i) entry[i].clear();
  }
  DumpElement& operator()(int ix, int iy=0, int iz=0)
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
    Real *pX, *pY, *pZ, *P, *fX, *fY, *fZ;
    Real *fxP, *fyP, *fzP, *fxV, *fyV, *fzV;
    Real *vx, *vy, *vz, *vxDef, *vyDef, *vzDef;
    Real *chi, *thrust, *pThrust, *pDef;
    int Ndata, nAlloc, nMapped, *gridMap;
    vector<surfData*> Set;

    surfacePoints() :
    	Ndata(0), nAlloc(0), nMapped(0), pX(nullptr), pY(nullptr), pZ(nullptr), P(nullptr),
    fX(nullptr), fY(nullptr), fZ(nullptr), fxP(nullptr), fyP(nullptr), fzP(nullptr), fxV(nullptr), fyV(nullptr), fzV(nullptr),
    vx(nullptr), vy(nullptr), vz(nullptr), vxDef(nullptr), vyDef(nullptr), vzDef(nullptr), gridMap(nullptr), chi(nullptr), thrust(nullptr), pThrust(nullptr), pDef(nullptr)
    { }

    ~surfacePoints()
    {
        if(pX      not_eq nullptr){delete[] pX;      pX=nullptr;     }
        if(pY      not_eq nullptr){delete[] pY;      pY=nullptr;     }
        if(pZ      not_eq nullptr){delete[] pZ;      pZ=nullptr;     }
        if(P       not_eq nullptr){delete[] P;       P=nullptr;	     }
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
        if(chi	   not_eq nullptr){delete[] chi;     chi=nullptr;    }
        if(thrust  not_eq nullptr){delete[] thrust;  thrust=nullptr; }
        if(pThrust not_eq nullptr){delete[] pThrust; pThrust=nullptr;}
        if(pDef    not_eq nullptr){delete[] pDef;    pDef=nullptr;   }
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
    		for(int i=0; i<blocksPerThread.size(); i++)
    			for(int j=0; j<blocksPerThread[i].Set.size(); j++)
    				_add(blocksPerThread[i].Set[j]);

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
            if(chi     not_eq nullptr){delete[] chi;     chi=nullptr;    }
            if(thrust  not_eq nullptr){delete[] thrust;  thrust=nullptr; }
            if(pThrust not_eq nullptr){delete[] pThrust; pThrust=nullptr;}
            if(pDef    not_eq nullptr){delete[] pDef;    pDef=nullptr;   }

            pX      = new Real[nAlloc]; pY      = new Real[nAlloc]; pZ      = new Real[nAlloc];
            fX      = new Real[nAlloc]; fY      = new Real[nAlloc]; fZ      = new Real[nAlloc];
            fxP     = new Real[nAlloc]; fyP     = new Real[nAlloc]; fzP     = new Real[nAlloc];
            fxV     = new Real[nAlloc]; fyV     = new Real[nAlloc]; fzV     = new Real[nAlloc];
            vx      = new Real[nAlloc]; vy      = new Real[nAlloc]; vz      = new Real[nAlloc];
            vxDef   = new Real[nAlloc]; vyDef   = new Real[nAlloc]; vzDef   = new Real[nAlloc];
            P       = new Real[nAlloc]; gridMap = new  int[nAlloc]; chi     = new Real[nAlloc];
	    thrust  = new Real[nAlloc]; pThrust = new Real[nAlloc]; pDef    = new Real[nAlloc];
        }

        #ifndef NDEBUG
    		int checksum = 0;
    		for(int i=0; i<blocksPerThread.size(); i++) checksum += blocksPerThread[i].Ndata;
    		assert(checksum==Ndata && Ndata <= nAlloc);
        #endif
    }
    /*
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
    */
    void print(const int ID, const int stepNumber, const int rank)
    {
        ofstream fileskin;
        char buf[500];
        sprintf(buf, "skinDistrib_%1d_%07d_rank%03d.txt", ID, stepNumber+1, rank);
        string filename(buf);
        fileskin.open(filename, ios::trunc);

            fileskin<< "x,y,z,chi,thrust,pDef,pThrust,press,fxP,fyP,fzP,fxV,fyV,fzV,vx,vy,vz,vxDef,vyDef,vzDef"<<endl;

        for(int i=0; i<Ndata; ++i) {
            fileskin<< pX[i]<<","<< pY[i]<<","<< pZ[i]<<","
                    <<chi[i]<<","<< thrust[i]<<","<< pDef[i]<<","<< pThrust[i]<< ","
                    <<P[i]<<","
                    <<fxP[i]<<","<<fyP[i]<<","<<fzP[i]<<","
                    <<fxV[i]<<","<<fyV[i]<<","<<fzV[i]<<","
                    << vx[i]<<","<< vy[i]<<","<< vz[i]<<","
                    <<vxDef[i]<<","<<vyDef[i]<<","<<vzDef[i]<<endl;
        }
        fileskin.close();
    }

    /*
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
    */
};

struct StateReward
{
    bool bRestart, bForgiving;
    int info, NpLatLine, stepId;
    Real t_next_comm, avg_wght, GoalDX;
    Real Xrel, Xabs, Xpov, Yrel, Yabs, Ypov, Theta, VxAvg, VyAvg, AvAvg;
    Real thExp, vxExp, vyExp, avExp, VxInst, VyInst, AvInst;
    Real Dist, Quad, RelAng, VX, VY, AV, ThetaAvg, ThetaVel;
    Real PoutBnd, Pout, defPowerBnd, defPower, EffPDefBnd, EffPDef, Pthrust, Pdrag, ToD;
    Real battery, ext_X, ext_Y, ext_Z;;

    vector<Real> VelNAbove, VelTAbove, VelNBelow, VelTBelow;
    vector<Real> FPAbove, FVAbove, FPBelow, FVBelow;
    vector<Real> PXAbove, PYAbove, PXBelow, PYBelow;
    vector<Real> raySight;

    StateReward(const int _NpLatLine = 0) :
    bRestart(false), bForgiving(true), NpLatLine(_NpLatLine), info(1),
    t_next_comm(1e6), GoalDX(0), battery(1)
    {
        avg_wght = Xrel = Xabs = Yrel = Yabs = Theta = VxAvg = VyAvg = AvAvg = 0;
        thExp = vxExp = vyExp = avExp = VxInst = VyInst = AvInst = 0;
        Dist = Quad = RelAng = VX = VY = AV = ThetaAvg = ThetaVel = 0;
        Pout = PoutBnd = defPower = defPowerBnd = EffPDef = EffPDefBnd = Pthrust = Pdrag = ToD = 0;

        VelNAbove.resize(NpLatLine); VelTAbove.resize(NpLatLine);
        VelNBelow.resize(NpLatLine); VelTBelow.resize(NpLatLine);

        FPAbove.resize(_NpLatLine); FVAbove.resize(_NpLatLine);
        FPBelow.resize(_NpLatLine); FVBelow.resize(_NpLatLine);

        PXAbove.resize(NpLatLine); PYAbove.resize(NpLatLine);
        PXBelow.resize(NpLatLine); PYBelow.resize(NpLatLine);
        raySight.resize(2*NpLatLine);

    }

    void set_NpLatLine(const int _NpLatLine)
    {
      NpLatLine = _NpLatLine;
      VelNAbove.resize(NpLatLine); VelTAbove.resize(NpLatLine);
      VelNBelow.resize(NpLatLine); VelTBelow.resize(NpLatLine);

      FPAbove.resize(_NpLatLine); FVAbove.resize(_NpLatLine);
      FPBelow.resize(_NpLatLine); FVBelow.resize(_NpLatLine);

      PXAbove.resize(NpLatLine); PYAbove.resize(NpLatLine);
      PXBelow.resize(NpLatLine); PYBelow.resize(NpLatLine);
      raySight.resize(2*NpLatLine);
    }

    void updateStepId(const int _stepId) {stepId=_stepId;}

    void updateAverages(const Real _dt, const Real _th, const Real _vx, const Real _vy,
                        const Real _av, const Real _pO1, const Real _pO2, const Real _pW1,
                        const Real _pW2, const Real _eff1, const Real _eff2, const Real _pT,
                        const Real _pD, const Real _T, const Real _D)
    {
        if (_dt<=0) return;
    		const Real    W = avg_wght + _dt;
    		const Real _ToD = (_D<1e-9) ? 0 : _T/_D;
    		const Real _1oW = 1./W;

    		VxAvg = ( VxAvg * avg_wght + _vx * _dt )*_1oW;
    		VyAvg = ( VyAvg * avg_wght + _vy * _dt )*_1oW;
    		AvAvg = ( AvAvg * avg_wght + _av * _dt )*_1oW;
    		ThetaAvg = ( ThetaAvg * avg_wght + _th            * _dt )*_1oW;
    		ThetaVel = ( ThetaVel * avg_wght + atan2(_vy,_vx) * _dt )*_1oW;
    		Pout =    (Pout    * avg_wght + _pO1 * _dt)*_1oW;
    		PoutBnd = (PoutBnd * avg_wght + _pO2 * _dt)*_1oW;
    		defPower =    (defPower    * avg_wght + _pW1 * _dt)*_1oW;
    		defPowerBnd = (defPowerBnd * avg_wght + _pW2 * _dt)*_1oW;
    		EffPDef =    (EffPDef    * avg_wght + _eff1 * _dt)*_1oW;
    		EffPDefBnd = (EffPDefBnd * avg_wght + _eff2 * _dt)*_1oW;
    		Pthrust = ( Pthrust * avg_wght + _pT  * _dt )*_1oW;
    		Pdrag =   ( Pdrag   * avg_wght + _pD  * _dt )*_1oW;
    		ToD =     ( ToD     * avg_wght + _ToD * _dt )*_1oW;
    		battery += defPowerBnd*_dt/1e-3; // roughly 30 swimming periods of solo guy (3.5e-5 ~ mean power per one swimming period)
    		avg_wght += _dt;

    		thExp = (1.-_dt) * thExp + _dt * _th;
    		vxExp = (1.-_dt) * vxExp + _dt * _vx;
    		vyExp = (1.-_dt) * vyExp + _dt * _vy;
    		avExp = (1.-_dt) * avExp + _dt * _av;
    }

    void resetAverage()
    {
        VxAvg = VyAvg = AvAvg = ThetaAvg = ThetaVel = Pout = PoutBnd = 0;
        defPower = defPowerBnd = EffPDef = EffPDefBnd = 0;
        Pthrust = Pdrag = ToD = avg_wght = 0;
    }

    void updateInstant(const Real _xR, const Real _xA, const Real _yR, const Real _yA,
                       const Real _th, const Real _vx, const Real _vy, const Real _av)
    {
        Xrel = _xR; Xabs = _xA; Yrel = _yR; Yabs = _yA; Theta= _th;
        VxInst=_vx; VyInst=_vy; AvInst=_av;
        if (Xrel<0.05 || Yrel<0.025)  bRestart = true;
        if (ext_X>0 && ext_X-Xrel<0.2)  bRestart = true;
        if (ext_Y>0 && ext_Y-Yrel<.025) bRestart = true;
    }

    void finalizePos(const Real  xFOR,   const Real yFOR, const Real thFOR,
                     const Real vxFOR,  const Real vyFOR, const Real avFOR,
                     const Real lscale, const Real tscale)
    {
      //velocity of reference from fish pov
      VX = (VxInst-vxFOR)*std::cos(Theta) + (VyInst-vyFOR)*std::sin(Theta);
      VY = (VyInst-vyFOR)*std::cos(Theta) - (VxInst-vxFOR)*std::sin(Theta);
      AV = (AvInst-avFOR);
      //velocity of fish in reference pov
      VxAvg = VxAvg*std::cos(thFOR) + VyAvg*std::sin(thFOR);
      VyAvg = VyAvg*std::cos(thFOR) - VxAvg*std::sin(thFOR);
      AvAvg = AvAvg;
      //position in reference frame
      Xpov = (Xrel-xFOR)*std::cos(thFOR) + (Yrel-yFOR)*std::sin(thFOR);
      Ypov = (Yrel-yFOR)*std::cos(thFOR) - (Xrel-xFOR)*std::sin(thFOR);
      RelAng = Theta - thFOR;

      Dist = std::sqrt(std::pow(Xrel-xFOR,2) + std::pow(Yrel-yFOR,2));
      Quad = std::atan2(Ypov,Xpov) - (Theta-thFOR);
    }

    void finalize(const Real xFOR,   const Real yFOR,  const Real thFOR,
                  const Real vxFOR,  const Real vyFOR, const Real avFOR,
                  const Real lscale, const Real tscale)
    {

        //EffPDefBnd, EffPDef, ToD are non dimensional already
    }

    bool checkFail(const Real xFOR,const Real yFOR,const Real thFOR,const Real ls)
    {
        bRestart = false;
        if (Xrel<0.05 || Yrel<0.025)  bRestart = true;
        if (ext_X>0 && ext_X-Xrel<0.2)  bRestart = true;
        if (ext_Y>0 && ext_Y-Yrel<.025) bRestart = true;
        if (bRestart) { printf("Out of bounds\n"); return true; }

        const Real _Xrel = (Xrel-xFOR)*std::cos(thFOR) + (Yrel-yFOR)*std::sin(thFOR);
        const Real _Yrel = (Yrel-yFOR)*std::cos(thFOR) - (Xrel-xFOR)*std::sin(thFOR);
        const Real _thRel= Theta - thFOR;
        const Real _Dist = std::sqrt(std::pow(Xrel-xFOR,2) + std::pow(Yrel-yFOR,2));

        if(not bForgiving) {
            bRestart = _Dist<0.25*ls;
            if(bRestart) {printf("Too close\n"); return bRestart;}
            //at DX=1, allowed DY=.5, at DX=2.5 allowed DY=.75
            bRestart = std::fabs(_Yrel) > _Xrel/6. + 7*ls/12.;
            if(bRestart) {printf("Too much vertical distance\n"); return bRestart;}

            #ifdef __ExploreHalfWake
              bRestart = _Yrel < -.1*ls;
              if(bRestart) {printf("Wrong half of the wake\n"); return bRestart;}
            #endif

            bRestart = std::fabs(_thRel)>1.5708;
            if(bRestart) {printf("Too different inclination\n"); return bRestart;}

            bRestart = _Xrel<ls || _Xrel >2.5*ls;
            if(bRestart) {printf("Too far from horizontal goal\n"); return bRestart;}
        } else {
            bRestart = _Dist<0.25*ls;
            if(bRestart) {printf("Too close\n"); return bRestart;}

            bRestart = std::fabs(_Yrel)>ls;
            if(bRestart) {printf("Too much vertical distance\n"); return bRestart;}

            bRestart = std::fabs(_thRel)>M_PI;
            if(bRestart) {printf("Too different inclination\n"); return bRestart;}

            bRestart = std::fabs(_Xrel-GoalDX*ls)>ls;
            if(bRestart) {printf("Too far from horizontal goal\n"); return bRestart;}
        }

        return bRestart;
    }

    struct skinForcesVels
    {
      skinForcesVels(const int _nDest) : nDest(_nDest), data(_alloc(7*_nDest))
      {
        memset(data, 0, sizeof(Real)*7*nDest);
      }

      virtual ~skinForcesVels() { _dealloc(data); }

      inline void storeNearest(const Real fxP, const Real fyP, const Real fxV,
                      const Real fyV, const Real vx, const Real vy, const int i)
      {
        data[i+0*nDest] += fxP; data[i+1*nDest] += fyP;
        data[i+2*nDest] += fxV; data[i+3*nDest] += fyV;
        data[i+4*nDest] += vx;  data[i+5*nDest] += vy;
        data[i+6*nDest] += 1.;
      }

      inline Real fxP(const int i) { return data[i+0*nDest]; }
      inline Real fyP(const int i) { return data[i+1*nDest]; }
      inline Real fxV(const int i) { return data[i+2*nDest]; }
      inline Real fyV(const int i) { return data[i+3*nDest]; }
      inline Real vx (const int i) { return data[i+4*nDest]; }
      inline Real vy (const int i) { return data[i+5*nDest]; }

      void synchronize(const MPI_Comm comm)
      {
        int rank;
        MPI_Comm_rank(comm, &rank);
        #ifndef _SP_COMP_
        MPI_Allreduce(MPI_IN_PLACE, data, 7*nDest, MPI_DOUBLE, MPI_SUM, comm);
        #else // _SP_COMP_
        MPI_Allreduce(MPI_IN_PLACE, data, 7*nDest, MPI_FLOAT,  MPI_SUM, comm);
        #endif //

        for(int i=0; i<nDest; ++i)
          if(data[i+6*nDest]) {
            //data[i+0*nDest] /= data[i+6*nDest];
            //data[i+1*nDest] /= data[i+6*nDest];
            //data[i+2*nDest] /= data[i+6*nDest];
            //data[i+3*nDest] /= data[i+6*nDest];
              data[i+4*nDest] /= data[i+6*nDest];
              data[i+5*nDest] /= data[i+6*nDest];
          } else if (0)//(!rank)
          {
            printf("Worryingly, some of the entries in skinForcesVels probably got no data\n");
            fflush(0);
          }
      }
      void print(const MPI_Comm comm, const int stepNumber)
      {
          int rank;
          MPI_Comm_rank(comm, &rank);
          if(rank) return;
          ofstream fileskin;
          char buf[500];
          sprintf(buf, "midplaneData_%07d.txt", stepNumber);
          string filename(buf);
          fileskin.open(filename, ios::trunc);
          int k=0;
          for(int i=0; i<nDest; ++i)
              fileskin<<fxP(i)<<"\t"<<fyP(i)<<"\t"
                      <<fxV(i)<<"\t"<<fyV(i)<<"\t"
                      <<vx(i)<<"\t"<<vy(i)<<std::endl;
          fileskin.close();
      }

     private:
      const int nDest;
      Real*const data;
      Real * _alloc(const int N) { return new Real[N]; }

      void _dealloc(Real * ptr) {
          if(ptr not_eq nullptr) {
              delete [] ptr;
              ptr=nullptr;
          }
      }
    };

    void nearestGridPoints(const surfacePoints*const surface, const int Nskin,
      const Real*const xU, const Real*const yU, const Real*const xL, const Real*const yL,
      const Real*const nxU,const Real*const nyU,const Real*const nxL,const Real*const nyL,
      const Real zObst, const Real h, const MPI_Comm comm)
    {
        const int Nsurf = surface->Ndata;
        skinForcesVels data(Nskin*2);
        const Real*const xS    = surface->pX;
        const Real*const yS    = surface->pY;
        const Real*const zS    = surface->pZ;
        const Real*const fxP   = surface->fxP;
        const Real*const fyP   = surface->fyP;
        const Real*const fxV   = surface->fxV;
        const Real*const fyV   = surface->fyV;
        const Real*const vx    = surface->vx;
        const Real*const vy    = surface->vy;
        const Real*const vxDef = surface->vxDef;
        const Real*const vyDef = surface->vyDef;
        const Real eps = std::numeric_limits<Real>::epsilon();

        for (int j=0; j<Nsurf; j++) {
          if(std::fabs(zS[j]-zObst)>h+eps) continue;

          for (int i=0; i<Nskin; i++) {
            if(std::fabs(yS[j]-yU[i])<=h+eps && std::fabs(xS[j]-xU[i])<=h+eps)
            {
              const Real sensorVX = vx[j]-vxDef[j]-VxInst+(yS[j]-Yrel)*AvInst;
              const Real sensorVY = vy[j]-vyDef[j]-VyInst-(xS[j]-Xrel)*AvInst;
              #ifdef __DumpWakeStefan
              #warning "RL velocity sensor states are hijacked"
              data.storeNearest(fxP[j], fyP[j], fxV[j], fyV[j], vxDef[j], vyDef[j], i);
              #else
              data.storeNearest(fxP[j], fyP[j], fxV[j], fyV[j], sensorVX, sensorVY, i);
              //data.storeNearest(xU[i], yU[i], xS[j], yS[j], i); //for debug
              #endif
            }

            if(std::fabs(yS[j]-yL[i])<=h+eps && std::fabs(xS[j]-xL[i])<=h+eps)
            {
              const Real sensorVX = vx[j]-vxDef[j]-VxInst+(yS[j]-Yrel)*AvInst;
              const Real sensorVY = vy[j]-vyDef[j]-VyInst-(xS[j]-Xrel)*AvInst;
              #ifdef __DumpWakeStefan
              data.storeNearest(fxP[j], fyP[j], fxV[j], fyV[j], vxDef[j], vyDef[j], i+Nskin);
              #else
              data.storeNearest(fxP[j], fyP[j], fxV[j], fyV[j], sensorVX, sensorVY, i+Nskin);
              //data.storeNearest(xL[i], yL[i], xS[j], yS[j], i+Nskin); //for debug
              #endif
            }
          }
        }

        /*
        //we know the spacing between non zero grad chi points: exploit it
        const Real nearZ  = (std::round(zObst/h -.5) +.5)*h; //cubiz be cell centered
        //for each skin point, compute what is the closest non zero gradchi point
        #pragma omp parallel for
        for (int i=0; i<Nskin; i++) {
          const int k = i+Nskin; //index for lower skin
          const Real nearXU = (std::round(xU[i]/h -.5) +.5)*h;
          const Real nearYU = (std::round(yU[i]/h -.5) +.5)*h;
          const Real nearXL = (std::round(xL[i]/h -.5) +.5)*h;
          const Real nearYL = (std::round(yL[i]/h -.5) +.5)*h;

          for (int j=0; j<Nsurf; j++) {
            if(std::fabs(zS[j]-nearZ )>.5*h) continue;

            if(std::fabs(yS[j]-nearYL)<=.5*h && std::fabs(xS[j]-nearXL)<=.5*h){
              const Real sensorVX = vx[j]-vxDef[j]-VxInst+(yS[j]-Yrel)*AvInst;
              const Real sensorVY = vy[j]-vyDef[j]-VyInst-(xS[j]-Xrel)*AvInst;
              //data.storeNearest(fx[j], fy[j], sensorVX, sensorVY, k);
              data.storeNearest(xL[i], yL[i], xS[j], yS[j], k);
            }

            if(std::fabs(yS[j]-nearYU)<=.5*h && std::fabs(xS[j]-nearXU)<=.5*h){
              const Real sensorVX = vx[j]-vxDef[j]-VxInst+(yS[j]-Yrel)*AvInst;
              const Real sensorVY = vy[j]-vyDef[j]-VyInst-(xS[j]-Xrel)*AvInst;
              //data.storeNearest(fx[j], fy[j], sensorVX, sensorVY, i);
              data.storeNearest(xU[i], yU[i], xS[j], yS[j], i);
            }
          }
        }
        */

        data.synchronize(comm);
        //data.print(comm,stepId);

        /*
          int rank;
          MPI_Comm_rank(comm, &rank);
          if(!rank) {
            ofstream fileskin;
            char buf[500];
            sprintf(buf, "skinPoints_%07d.txt", stepId);
            string filename(buf);
            fileskin.open(filename, ios::trunc);
            for (int j=0;       j<Nskin; j++)
                fileskin<<xU[j]<<"\t"<<yU[j]<<std::endl;
            for (int j=Nskin-1; j>=0;    j--)
                fileskin<<xL[j]<<"\t"<<yL[j]<<std::endl;
            fileskin.close();
          }
        */
        vector<Real> NxAbove(NpLatLine), NyAbove(NpLatLine);
        vector<Real> NxBelow(NpLatLine), NyBelow(NpLatLine);
        //now, feed the sensors
        for (int k=0; k<NpLatLine; k++)
        {
            const int first = k   *(Real)Nskin/(Real)NpLatLine;
            const int last = (k+1)*(Real)Nskin/(Real)NpLatLine;
            Real VelXAbove=0,VelXBelow=0,VelYAbove=0,VelYBelow=0;
            Real  FPxAbove=0, FPxBelow=0, FPyAbove=0, FPyBelow=0;
            Real  FVxAbove=0, FVxBelow=0, FVyAbove=0, FVyBelow=0;

            for (int j=first; j<last; j++) {
                FPxAbove+=data.fxP(j); FPxBelow+=data.fxP(j+Nskin);
                FPyAbove+=data.fyP(j); FPyBelow+=data.fyP(j+Nskin);

                FVxAbove+=data.fxV(j); FVxBelow+=data.fxV(j+Nskin);
                FVyAbove+=data.fyV(j); FVyBelow+=data.fyV(j+Nskin);

                VelXAbove+=data.vx(j); VelXBelow+=data.vx(j+Nskin);
                VelYAbove+=data.vy(j); VelYBelow+=data.vy(j+Nskin);
            }
            const Real fac = 1./(Real)(last-first);
            VelXAbove*=fac; VelYAbove*=fac;
            VelXBelow*=fac; VelYBelow*=fac;

            const int mid = 0.5*(first+last);
            PXAbove[k] = xU[mid]; PYAbove[k] = yU[mid];
            PXBelow[k] = xL[mid]; PYBelow[k] = yL[mid];

            const Real nxAbove = nxU[mid]; // ^            ^
            const Real nyAbove = nyU[mid]; //   `        /
            const Real txAbove = nyU[mid]; //   n `    / t
            const Real tyAbove =-nxU[mid]; //       `
            NxAbove[k] = nxAbove;
            NyAbove[k] = nyAbove;
            const Real nxBelow = nxL[mid]; //    /`
            const Real nyBelow = nyL[mid]; // n /    `  t
            const Real txBelow =-nyL[mid]; //  /        `
            const Real tyBelow = nxL[mid]; // v            v
            NxBelow[k] = nxBelow;
            NyBelow[k] = nyBelow;

            #ifdef __DumpWakeStefan //i want to keep x and y
            VelNAbove[k] = VelXAbove;
            VelTAbove[k] = VelYAbove;

            VelNBelow[k] = VelXBelow;
            VelTBelow[k] = VelYBelow;
            #else
            VelNAbove[k] = VelXAbove*nxAbove + VelYAbove*nyAbove;
            VelTAbove[k] = VelXAbove*txAbove + VelYAbove*tyAbove;

            VelNBelow[k] = VelXBelow*nxBelow + VelYBelow*nyBelow;
            VelTBelow[k] = VelXBelow*txBelow + VelYBelow*tyBelow;
            #endif

            FPAbove[k] = FPxAbove*nxAbove + FPyAbove*nyAbove;
            FVAbove[k] = FVxAbove*txAbove + FVyAbove*tyAbove;

            FPBelow[k] = FPxBelow*nxBelow + FPyBelow*nyBelow;
            FVBelow[k] = FVxBelow*txBelow + FVyBelow*tyBelow;
        }

        if(0){
            ofstream fileskin;
            char buf[500];
            sprintf(buf, "sensorDistrib_%07d.txt", stepId);
            string filename(buf);
            fileskin.open(filename, ios::trunc);
            int k=0;
            for(int i=0; i<NpLatLine; ++i)
                fileskin<<  PXAbove[i]<<"\t"<<  PYAbove[i]<<"\t"
                        <<  NxAbove[i]<<"\t"<<  NyAbove[i]<<"\t"
                        <<VelNAbove[i]<<"\t"<<VelTAbove[i]<<"\t"
                        <<  FPAbove[i]<<"\t"<<  FVAbove[i]<<"\t"
                        << raySight[k++] << std::endl;
            for(int i=0; i<NpLatLine; ++i)
                fileskin<<  PXBelow[i]<<"\t"<<  PYBelow[i]<<"\t"
                        <<  NxBelow[i]<<"\t"<<  NyBelow[i]<<"\t"
                        <<VelNBelow[i]<<"\t"<<VelTBelow[i]<<"\t"
                        <<  FPBelow[i]<<"\t"<<  FVBelow[i]<<"\t"
                        << raySight[k++] << std::endl;
            fileskin.close();
        }
    }

    void save(const int step_id, string filename)
    {
        ofstream savestream;
        savestream.setf(std::ios::scientific);
        savestream.precision(std::numeric_limits<double>::digits10 + 1);
        string fullFileName = filename==string() ? "restart_IF2D_Stefan" : filename;

        savestream.open(fullFileName+"_save_data.txt");

        savestream << bRestart << "\t" << info << "\t" << avg_wght << "\t" << t_next_comm << "\t"
            << Xrel << "\t" << Xabs << "\t" << Yrel << "\t" << Yabs << "\t"
            << Theta << "\t" << VxAvg << "\t" << VyAvg<< "\t" << AvAvg << "\t"
            << thExp << "\t" << vxExp << "\t" << vyExp<< "\t" << avExp << "\t"
            << VxInst << "\t" << VyInst<< "\t" << AvInst << "\t"
            << Dist << "\t" << Quad << "\t" << RelAng<< "\t"
            << VX << "\t" << VY << "\t" << AV << "\t"
            << ThetaAvg << "\t" << ThetaVel << "\t" << PoutBnd << "\t" << Pout << "\t"
            << defPowerBnd << "\t" << defPower << "\t" << EffPDefBnd<< "\t" << EffPDef << "\t"
            << Pthrust << "\t" << Pdrag << "\t" << ToD << std::endl;

        for (int i=0; i<NpLatLine; i++) {
            savestream <<
            PXAbove[i] << "\t" << PYAbove[i] << "\t" <<
            PXBelow[i] << "\t" << PYBelow[i] << "\t" <<
            VelNAbove[i] << "\t" << VelTAbove[i] << "\t" <<
            VelNBelow[i] << "\t" << VelTBelow[i] << "\t" <<
            FPAbove[i] << "\t" << FVAbove[i] << "\t" <<
            FPBelow[i] << "\t" << FVBelow[i] << std::endl;
        }

        savestream.close();
    }

    void restart(string filename)
    {
        ifstream restartstream;
        string fullFileName = filename;
        restartstream.open(fullFileName+"_save_data.txt");

        restartstream >> bRestart >> info >> avg_wght >> t_next_comm >>
        Xrel >> Xabs >> Yrel >> Yabs >>
        Theta >> VxAvg >> VyAvg >> AvAvg >>
        thExp >> vxExp >> vyExp >> avExp >>
        VxInst >> VyInst >> AvInst >>
        Dist >> Quad >> RelAng >>
        VX >> VY >> AV >>
        ThetaAvg >> ThetaVel >> PoutBnd >> Pout >>
        defPowerBnd >> defPower >> EffPDefBnd>> EffPDef >>
        Pthrust >> Pdrag >> ToD;

        for (int i=0; i<NpLatLine; i++) {
            restartstream >>
            PXAbove[i] >> PYAbove[i] >>
            PXBelow[i] >> PYBelow[i] >>
            VelNAbove[i] >> VelTAbove[i] >>
            VelNBelow[i] >> VelTBelow[i] >>
            FPAbove[i] >> FVAbove[i] >>
            FPBelow[i] >> FVBelow[i];
        }

        restartstream.close();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0)
        {
            cout << bRestart << "\t" << info << "\t" << avg_wght << "\t" << t_next_comm << "\t"
            << Xrel << "\t" << Xabs << "\t" << Yrel << "\t" << Yabs << "\t"
            << Theta << "\t" << VxAvg << "\t" << VyAvg<< "\t" << AvAvg << "\t"
            << thExp << "\t" << vxExp << "\t" << vyExp<< "\t" << avExp << "\t"
            << VxInst << "\t" << VyInst<< "\t" << AvInst << "\t"
            << Dist << "\t" << Quad << "\t" << RelAng<< "\t"
            << VX << "\t" << VY << "\t" << AV << "\t"
            << ThetaAvg << "\t" << ThetaVel << "\t" << PoutBnd << "\t" << Pout << "\t"
            << defPowerBnd << "\t" << defPower << "\t" << EffPDefBnd<< "\t" << EffPDef << "\t"
            << Pthrust << "\t" << Pdrag << "\t" << ToD << std::endl;

            for (int i=0; i<NpLatLine; i++) {
            cout << PXAbove[i] << "\t" << PYAbove[i] << "\t" <<
                    PXBelow[i] << "\t" << PYBelow[i] << "\t" <<
                    VelNAbove[i] << "\t" << VelTAbove[i] << "\t" <<
                    VelNBelow[i] << "\t" << VelTBelow[i] << "\t" <<
                    FPAbove[i] << "\t" << FVAbove[i] << "\t" <<
                    FPBelow[i] << "\t" << FVBelow[i] << std::endl;
            }
        }
    }

    void print(const int ID, const int stepNumber, const Real time)
    {
	//int rank;
        //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        //if (rank) return;
        {
            ofstream fileskin;
            char buf[500];
            sprintf(buf, "sensorDistrib_%1d_%07d.txt", ID, stepNumber);
            string filename(buf);
            fileskin.open(filename, ios::trunc);
            int k=0;
            for(int i=0; i<NpLatLine; ++i)
                fileskin<<PXAbove[i]<<"\t"<<PYAbove[i]<<"\t"<<VelNAbove[i] << "\t" << VelTAbove[i] << "\t" << FPAbove[i] << "\t" << FVAbove[i] << "\t" << raySight[k++] << std::endl;
            for(int i=0; i<NpLatLine; ++i)
                fileskin<<PXBelow[i]<<"\t"<<PYBelow[i]<<"\t"<<VelNBelow[i] << "\t" << VelTBelow[i] << "\t" << FPBelow[i] << "\t" << FVBelow[i] << "\t" << raySight[k++] << std::endl;
            fileskin.close();
        }
        {
            ofstream fileskin;
            char buf[500];
            sprintf(buf, "avgSensors_%1d.txt",ID);
            string filename(buf);
            fileskin.open(filename, ios::app);

            fileskin<< avg_wght << "\t" << t_next_comm  << "\t"
                    << Xrel << "\t" << Xabs << "\t" << Yrel << "\t" << Yabs << "\t"
                    << Theta << "\t" << VxAvg << "\t" << VyAvg<< "\t" << AvAvg << "\t"
                    << thExp << "\t" << vxExp << "\t" << vyExp<< "\t" << avExp << "\t"
                    << VxInst << "\t" << VyInst<< "\t" << AvInst << "\t"
                    << Dist << "\t" << Quad << "\t" << RelAng<< "\t"
                    << VX << "\t" << VY << "\t" << AV << "\t"
                    << ThetaAvg << "\t" << ThetaVel << "\t" << PoutBnd << "\t" << Pout << "\t"
                    << defPowerBnd << "\t" << defPower << "\t" << EffPDefBnd<< "\t" << EffPDef << "\t"
                    << Pthrust << "\t" << Pdrag << "\t" << ToD << std::endl;
            fileskin.close();
        }
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
struct ChiStreamer
{
	static const int channels = 1;

	FluidBlock* ref;
	ChiStreamer(FluidBlock& b): ref(&b) {}
  ChiStreamer(): ref(NULL) {}
	inline Real operate(const int ix, const int iy, const int iz) const
	{
    const FluidElement& input = ref->data[iz][iy][ix];
		return input.chi;
	}
  template<int channel>
	static inline Real operate(const FluidElement& input)
  {
    printf("Abort in ChiStreamer\n");
    abort(); return 0;
  }

	const char * name() { return "ChiStreamer"; }
};
/*
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

	const char * name() { return "StreamerGridPointIterative" ; }
};

template<> inline Real FluidVPStreamer::operate<0>(const FluidElement& e)
{ return e.u; }
template<> inline Real FluidVPStreamer::operate<1>(const FluidElement& e)
{ return e.v; }
template<> inline Real FluidVPStreamer::operate<2>(const FluidElement& e)
{ return e.w; }
template<> inline Real FluidVPStreamer::operate<3>(const FluidElement& e)
{ return e.chi; }
template<> inline Real FluidVPStreamer::operate<4>(const FluidElement& e)
{ return e.p; }
template<> inline Real FluidVPStreamer::operate<5>(const FluidElement& e)
{ return e.tmpU; }
template<> inline Real FluidVPStreamer::operate<6>(const FluidElement& e)
{ return e.tmpV; }
template<> inline Real FluidVPStreamer::operate<7>(const FluidElement& e)
{ return e.tmpW; }

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

*/
struct StreamerHDF5Dump
{
	static const int NCHANNELS = 4;

	DumpBlock& ref;

	StreamerHDF5Dump(DumpBlock& b): ref(b) {}

	void operate(const int ix, const int iy, const int iz, Real output[NCHANNELS]) const
	{
		const DumpElement& input = ref.data[iz][iy][ix];
		output[0] = input.u;
		output[1] = input.v;
		output[2] = input.w;
		output[3] = input.chi;
	}

	void operate(const Real input[NCHANNELS], const int ix, const int iy, const int iz) const
	{
		DumpElement& output = ref.data[iz][iy][ix];
		output.u    = input[0];
		output.v    = input[1];
		output.w    = input[2];
		output.chi  = input[3];
	}

	inline void dump(const int ix, const int iy, const int iz, float* const ovalue, const int field) const
	{
		const DumpElement& input = ref.data[iz][iy][ix];
		switch(field) {
			case 0: *ovalue  = input.u; break;
			case 1: *ovalue  = input.v; break;
			case 2: *ovalue  = input.w; break;
			case 3: *ovalue  = input.chi; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}

  template<int field>
	inline void load(const Real ivalue, const int ix, const int iy, const int iz) const
	{
		DumpElement& output = ref.data[iz][iy][ix];
		switch(field) {
			case 0:  output.u    = ivalue; break;
			case 1:  output.v    = ivalue; break;
			case 2:  output.w    = ivalue; break;
			case 3:  output.chi  = ivalue; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}

	static const char * getAttributeName() { return "Vector"; }
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

	inline void dump(const int ix, const int iy, const int iz, float* const ovalue, const int field) const
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

  template<int field>
	inline void load(const Real ivalue, const int ix, const int iy, const int iz) const
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
/*
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
*/
template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabOpen: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;
  public:
    //virtual inline std::string name() const { return "BlockLabOpen"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}
    BlockLabOpen(): BlockLab<BlockType,allocator>(){}
    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator>
                bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
        //if (info.index[0]==0)           bc.template applyBC_dirichlet_inflow<0,0>();
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};

typedef Grid<FluidBlock, std::allocator> FluidGrid;
typedef GridMPI<FluidGrid> FluidGridMPI;

typedef Grid<DumpBlock, std::allocator> DumpGrid;
typedef GridMPI<DumpGrid> DumpGridMPI;

//#define _OPEN_BC_
//#ifndef _OPEN_BC_
//typedef  Lab;
//#else
typedef  BlockLabOpen<FluidBlock, std::allocator> Lab;
//#endif

//#ifdef _MIXED_
//typedef BlockLabBottomWall<FluidBlock, std::allocator> Lab;
//#endif // _MIXED_

//#ifdef _VORTEX_
//typedef BlockLabVortex<FluidBlock, std::allocator> Lab;
//#endif // _VORTEX_

//#ifdef _PIPE_
//typedef BlockLabPipe<FluidBlock, std::allocator> Lab;
//#endif // _PIPE_

typedef BlockLabMPI<BlockLabOpen<FluidBlock, std::allocator>> LabMPI;

/*
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
*/
#endif
