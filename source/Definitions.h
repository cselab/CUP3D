//
//  DataStructures.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_DataStructures_h
#define CubismUP_3D_DataStructures_h

#include "common.h"
#include "Layer.h"
#include "BoundaryConditions.h"

#ifndef _BS_
#define _BS_ 32
#endif // _BS_

struct FluidElement
{
    Real rho, u, v, w, chi, p, pOld;
	Real tmpU, tmpV, tmpW, tmp;
	Real divU;
	Real x, y, z;
    
    FluidElement() : rho(0), u(0), v(0), w(0), chi(0), p(0), pOld(0), divU(0), tmpU(0), tmpV(0), tmpW(0), tmp(0), x(0), y(0), z(0) {}
    
    void clear()
    {
        rho = u = v = w = chi = p = pOld = 0;
		tmpU = tmpV = tmpW = tmp = 0;
		divU = 0;
		x = y = z = 0;
    }
};

struct FluidVTKStreamer
{
	static const int channels = 9;
	
	void operate(FluidElement input, Real output[9])
	{
		output[0] = input.rho;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.w;
		output[4] = input.p;
		output[5] = input.chi;
		output[6] = input.divU;
		output[7] = input.pOld;
		output[8] = input.tmp;
	}
};

// this is used for serialization - important that ALL the quantities are streamed
struct StreamerGridPoint
{
	static const int channels = 12;
	
	void operate(const FluidElement& input, Real output[12]) const
	{
		abort();
		output[0] = input.rho;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.w;
		output[4] = input.chi;
		output[5] = input.p;
		output[6] = input.pOld;
		output[7] = input.tmpU;
		output[8] = input.tmpV;
		output[8] = input.tmpW;
		output[10] = input.tmp;
		output[11] = input.divU;
	}
	
	void operate(const Real input[12], FluidElement& output) const
	{
		abort();
		output.rho  = input[0];
		output.u    = input[1];
		output.v    = input[2];
		output.w    = input[3];
		output.chi  = input[4];
		output.p    = input[5];
		output.pOld = input[6];
		output.tmpU = input[7];
		output.tmpV = input[8];
		output.tmpW = input[9];
		output.tmp  = input[10];
		output.divU = input[11];
	}
};

struct StreamerGridPointASCII
{
	void operate(const FluidElement& input, ofstream& output) const
	{
		output << input.rho << " " << input.u << " " << input.v << " " << input.w << " " << input.chi << " " << input.p << " " << input.pOld << " " << input.tmpU << " " << input.tmpV << " " << input.tmpW << " " << input.tmp << " " << input.divU;
	}
	
	void operate(ifstream& input, FluidElement& output) const
	{
		input >> output.rho;
		input >> output.u;
		input >> output.v;
		input >> output.w;
		input >> output.chi;
		input >> output.p;
		input >> output.pOld;
		input >> output.tmpU;
		input >> output.tmpV;
		input >> output.tmpW;
		input >> output.tmp;
		input >> output.divU;
	}
};

struct StreamerDiv
{
	static const int channels = 1;
	static void operate(const FluidElement& input, Real output[1])
	{
		output[0] = input.divU;
    }
    
    static void operate(const Real input[1], FluidElement& output)
    {
        output.divU = input[0];
    }
    
    static void setTmp(const Real input[1], FluidElement& output)
    {
        output.tmp = input[0];
    }
};

struct FluidBlock
{
    //these identifiers are required by cubism!
    static const int sizeX = _BS_;
    static const int sizeY = _BS_;
    static const int sizeZ = _BS_;
    typedef FluidElement ElementType;
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
	static const int channels = 6;
	
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

template<> inline Real FluidVPStreamer::operate<0>(const FluidElement& e) { return e.rho; }
template<> inline Real FluidVPStreamer::operate<1>(const FluidElement& e) { return e.u; }
template<> inline Real FluidVPStreamer::operate<2>(const FluidElement& e) { return e.v; }
template<> inline Real FluidVPStreamer::operate<3>(const FluidElement& e) { return e.w; }
template<> inline Real FluidVPStreamer::operate<4>(const FluidElement& e) { return e.p; }
template<> inline Real FluidVPStreamer::operate<5>(const FluidElement& e) { return e.chi; }

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
{
	static const int NCHANNELS = 12;
	
	FluidBlock& ref;
	
	StreamerSerialization(FluidBlock& b): ref(b) {}
	
	void operate(const int ix, const int iy, const int iz, Real output[12]) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		output[0]  = input.rho;
		output[1]  = input.u;
		output[2]  = input.v;
		output[3]  = input.w;
		output[4]  = input.chi;
		output[5]  = input.p;
		output[6]  = input.pOld;
		output[7]  = input.tmpU;
		output[8]  = input.tmpV;
		output[9]  = input.tmpW;
		output[10] = input.tmp;
		output[11] = input.divU;
	}
	
	void operate(const Real input[12], const int ix, const int iy, const int iz) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		output.rho  = input[0];
		output.u    = input[1];
		output.v    = input[2];
		output.w    = input[3];
		output.chi  = input[4];
		output.p    = input[5];
		output.pOld = input[6];
		output.tmpU = input[7];
		output.tmpV = input[8];
		output.tmpW = input[9];
		output.tmp  = input[10];
		output.divU = input[11];
	}
	
	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0: *ovalue  = input.rho; break;
			case 1: *ovalue  = input.u; break;
			case 2: *ovalue  = input.v; break;
			case 3: *ovalue  = input.w; break;
			case 4: *ovalue  = input.chi; break;
			case 5: *ovalue  = input.p; break;
			case 6: *ovalue  = input.pOld; break;
			case 7: *ovalue  = input.tmpU; break;
			case 8: *ovalue  = input.tmpV; break;
			case 9: *ovalue  = input.tmpW; break;
			case 10: *ovalue = input.tmp; break;
			case 11: *ovalue = input.divU; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0:  output.rho  = ivalue; break;
			case 1:  output.u    = ivalue; break;
			case 2:  output.v    = ivalue; break;
			case 3:  output.w    = ivalue; break;
			case 4:  output.chi  = ivalue; break;
			case 5:  output.p    = ivalue; break;
			case 6:  output.pOld = ivalue; break;
			case 7:  output.tmpU = ivalue; break;
			case 8:  output.tmpV = ivalue; break;
			case 9:  output.tmpW = ivalue; break;
			case 10: output.tmp  = ivalue; break;
			case 11: output.divU = ivalue; break;
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
		
		output[0] = input.rho;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.w;
		output[4] = input.chi;
		output[5] = input.p;
		output[6] = input.tmp;
		output[7] = 0;
		output[8] = 0;
	}
	
	void operate(const Real input[NCHANNELS], const int ix, const int iy, const int iz) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		output.rho  = input[0];
		output.u    = input[1];
		output.v    = input[2];
		output.w    = input[3];
		output.chi  = input[4];
		output.p    = input[5];
		output.tmp  = input[6];
	}
	
	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0: *ovalue  = input.rho; break;
			case 1: *ovalue  = input.u; break;
			case 2: *ovalue  = input.v; break;
			case 3: *ovalue  = input.w; break;
			case 4: *ovalue  = input.chi; break;
			case 5: *ovalue  = input.p; break;
			case 6: *ovalue  = input.tmp; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0:  output.rho  = ivalue; break;
			case 1:  output.u    = ivalue; break;
			case 2:  output.v    = ivalue; break;
			case 3:  output.w    = ivalue; break;
			case 4:  output.chi  = ivalue; break;
			case 5:  output.p    = ivalue; break;
			case 6:  output.tmp = ivalue; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	static const char * getAttributeName() { return "Tensor"; }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBottomWall : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
    ElementTypeBlock pDirichlet;
    
	BlockLabBottomWall(): BlockLab<BlockType,allocator>()
    {
        pDirichlet.rho = 1;
        pDirichlet.chi = 0;
		pDirichlet.u = 0;
		pDirichlet.v = 0;
		pDirichlet.w = 0;
        pDirichlet.p = 0;
        pDirichlet.pOld = 0;
        pDirichlet.divU = 0;
        pDirichlet.tmp = 1;
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

typedef Grid<FluidBlock, std::allocator> FluidGrid;
typedef Grid<ScalarBlock, std::allocator> ScalarGrid;

#ifdef _MIXED_
typedef BlockLabBottomWall<FluidBlock, std::allocator> Lab;
#endif // _MIXED_

#ifdef _PERIODIC_
typedef BlockLab<FluidBlock, std::allocator> Lab;
#endif // _PERIODIC_

#ifdef _VORTEX_
typedef BlockLabVortex<FluidBlock, std::allocator> Lab;
#endif // _VORTEX_

#ifdef _PIPE_
typedef BlockLabPipe<FluidBlock, std::allocator> Lab;
#endif // _PIPE_




#endif
