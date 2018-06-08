//
//  CubismUP_3D
//
//  Written by Fabian Wermelinger.
//  Copyright (c) 2018 ETHZ. All rights reserved.
//
//  This algorithm uses the cyclic convolution method described in Eastwood and
//  Brownrigg (1979) for unbounded domains.  (In a cell-centered data
//  formulation).
//  WARNING: This implementation only works with a 1D domain decomposition
//  along the x-coordinate.
#ifndef POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV
#define POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV

#include <fstream>
#include <sstream>

#include <cassert>
#include <cstring>
#include <fftw3-mpi.h>

#include "Definitions.h"

#ifndef _FLOAT_PRECISION_
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#define MPIREAL MPI_DOUBLE
#else
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#define MPIREAL MPI_FLOAT
#endif

using namespace std;

#include <BlockInfo.h>


template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_MPI
{
    typedef typename TGrid::BlockType BlockType;
    TGrid& grid;

    // FFTW
    //mycomplex local_rhs, local_work;
    myplan fwd, bwd;

    // MPI
    MPI_Comm m_comm;
    int m_rank, m_size;
    Real* m_buffer; // MPI send/recv buffer
    mycomplex m_kernel; // Fourier transform of Green's function

    // helper types
    MPI_Datatype m_arylo_t;
    MPI_Datatype m_aryhi_t;

    // grid dimensions
    const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
    const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
    const size_t mybpd[3] = {
        static_cast<size_t>(grid.getResidentBlocksPerDimension(0)),
        static_cast<size_t>(grid.getResidentBlocksPerDimension(1)),
        static_cast<size_t>(grid.getResidentBlocksPerDimension(2))
    };
    const size_t gsize[3] = {
        static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]),
        static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]),
        static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])
    };
    const size_t gsize_0[3] = {
        2*gsize[0] - 1, // x-size for zero-padded arrays. (w/o periodic copy)
        2*gsize[1] - 1, // y-size for zero-padded arrays. (w/o periodic copy)
        2*gsize[2] - 1  // z-size for zero-padded arrays. (w/o periodic copy)
    };
    const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
    const size_t nz_hat = gsize[2]/2+1;
    const double norm_factor = 1./(gsize_0[0]*gsize_0[1]*gsize_0[2]);
    const double h = grid.getBlocksInfo().front().h_gridpoint;
    ptrdiff_t alloc_local=0;
    ptrdiff_t src_local_n0=0, src_local_0_start=0, src_local_n1=0, src_local_1_start=0;
    ptrdiff_t kern_local_n0=0, kern_local_0_start=0, kern_local_n1=0, kern_local_1_start=0;

    void _show_some() // testing some
    {
        ostringstream fname;
        fname << "rank" << m_rank << ".dat";
        ofstream fout(fname.str());
        for (size_t i = 0; i < src_local_n0; ++i)
        {
            for (size_t j = 0; j < gsize_0[1]; ++j)
            {
                for (size_t k = 0; k < 2*(gsize_0[2]/2+1); ++k)
                {
                    const size_t idx = (i*gsize_0[1] + j)* 2*(gsize_0[2]/2+1) + k;
                    fout << data[idx] << '\t';
                }
                fout << endl;
            }
            fout << endl;
            fout << endl;
        }
        fout.close();
    }

    // setup Green's function and compute transform, done once
    _initialize_green()
    {
    }


protected:

    void _fftw2cub(const Real * const out) const
    {
#pragma omp parallel for
        for(size_t i=0; i<local_infos.size(); ++i) {
            const BlockInfo info = local_infos[i];
            BlockType& b = *(BlockType*)info.ptrBlock;
            const size_t offset = _offset(info);

            for(int ix=0; ix<BlockType::sizeX; ix++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int iz=0; iz<BlockType::sizeZ; iz++) {
                        const size_t src_index = _dest(offset, iz, iy, ix);
                        assert(src_index>=0 && src_index<myN[0]*myN[1]*myN[2]);
                        b(ix,iy,iz).p = out[src_index];
                    }
        }
    }

    void _mpi_send_src()
    {
        // distribute the source data to the zero-padded arrays.
        vector<MPI_Request> recvreq(m_size);
        if (m_rank < m_size/2)
        {
            int src_lo = m_rank*2 + 0;
            int src_hi = m_rank*2 + 1;
            MPI_Irecv(data, 1, m_arylo_t, src_lo, 0, m_comm, &recvreq[src_lo]);
            MPI_Irecv(data, 1, m_aryhi_t, src_hi, 1, m_comm, &recvreq[src_hi]);
        }

        int dst = m_rank/2;
        MPI_Send(m_buffer, (int)myN[0]*myN[1]*myN[2], MPIREAL, dst, m_rank%2, m_comm);

        MPI_Barrier(m_comm);
    }

    void _mpi_recv_sol()
    {
        // receive back solution
        vector<MPI_Request> sendreq(m_size);
        if (m_rank < m_size/2)
        {
            int dst_lo = m_rank*2 + 0;
            int dst_hi = m_rank*2 + 1;
            MPI_Isend(data, 1, m_arylo_t, dst_lo, 0, m_comm, &sendreq[dst_lo]);
            MPI_Isend(data, 1, m_aryhi_t, dst_hi, 1, m_comm, &sendreq[dst_hi]);
        }

        int src = m_rank/2;
        MPI_Recv(m_buffer, (int)myN[0]*myN[1]*myN[2], MPIREAL, src, m_rank%2, m_comm, MPI_STATUSES_IGNORE);
    }

    void _solve()
    {
        mycomplex *const in_out = (mycomplex *) data;
#if 0
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;

#pragma omp parallel for
        for(int j = 0; j<local_n1; ++j)
            for(int i = 0; i<gsize[0]; ++i)
                for(int k = 0; k<nz_hat; ++k) {
                    const int linidx = (j*gsize[0] +i)*nz_hat + k;
                    assert(linidx >=0 && linidx<nx*local_n1*nz_hat);
                    assert(linidx < alloc_local);

                    const Real denom = 32.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*(local_1_start+j)/ny) + cos(2.*M_PI*k/nz)) -
                    2.*(cos(4.*M_PI*i/nx) + cos(4.*M_PI*(local_1_start+j)/ny) + cos(4.*M_PI*k/nz)) - 90.;

                    const Real inv_denom = (denom==0)? 0.:1./denom;
                    const Real fatfactor = 12. * inv_denom * factor;

                    in_out[linidx][0] *= fatfactor;
                    in_out[linidx][1] *= fatfactor;
                }
#else
        const Real waveFactX = 2.0*M_PI/(gsize[0]*h);
        const Real waveFactY = 2.0*M_PI/(gsize[1]*h);
        const Real waveFactZ = 2.0*M_PI/(gsize[2]*h);
        const long nKx = static_cast<long>(gsize[0]);
        const long nKy = static_cast<long>(gsize[1]);
        const long nKz = static_cast<long>(gsize[2]);
        const long shifty = static_cast<long>(src_local_1_start);
#pragma omp parallel for
        for(long j = 0; j<static_cast<long>(src_local_n1); ++j)
            for(long i = 0; i<static_cast<long>(gsize[0]); ++i)
                for(long k = 0; k<static_cast<long>(nz_hat);   ++k) {
                    const size_t linidx = (j*gsize[0] +i)*nz_hat + k;
                    const long kx = (i <= nKx/2) ? i : -(nKx-i);
                    const long l = shifty + j; //memory index plus shift due to decomp
                    const long ky = (l <= nKy/2) ? l : -(nKy-l);
                    const long kz = (k <= nKz/2) ? k : -(nKz-k);

                    const Real rkx = kx*waveFactX, rky = ky*waveFactY, rkz = kz*waveFactZ;
                    const Real kinv = kx || ky || kz ? -1/(rkx*rkx+rky*rky+rkz*rkz) : 0;
                    in_out[linidx][0] *= kinv*norm_factor;
                    in_out[linidx][1] *= kinv*norm_factor;
                }
#endif

        //this is sparta!
        if (src_local_1_start == 0) in_out[0][0] = in_out[0][1] = 0;
    }

public:
    Real * data = nullptr;

    PoissonSolverScalarFFTW_MPI(const int desired_threads, TGrid& g): grid(g)
    {
        if (TStreamer::channels != 1) {
            cout << "PoissonSolverScalar_MPI(): Error: TStreamer::channels is "
            << TStreamer::channels << " (should be 1).\n";
            abort();
        }
        m_comm = grid.getCartComm();
        MPI_Comm_rank(m_comm, &m_rank);
        MPI_Comm_size(m_comm, &m_size);

        {
            int supported_threads;
            MPI_Query_thread(&supported_threads);
            if (supported_threads<MPI_THREAD_FUNNELED) {
                cout << "MPI implementation does not support threads.\n";
                abort();
            }
        }
#ifndef _FLOAT_PRECISION_
        const int retval = fftw_init_threads();
#else
        const int retval = fftwf_init_threads();
#endif
        if(retval==0) {
            cout << "FFTWBase::setup(): Call to fftw_init_threads() returned zero.\n";
            abort();
        }

#ifndef _FLOAT_PRECISION_
        fftw_plan_with_nthreads(desired_threads);
        fftw_mpi_init();
        alloc_local = fftw_mpi_local_size_3d_transposed(
                gsize_0[0], gsize_0[1], gsize_0[2]/2+1, m_comm,
                &src_local_n0, &src_local_0_start, &src_local_n1, &src_local_1_start);

        data = fftw_alloc_real(2*alloc_local); // source terms
        fwd = fftw_mpi_plan_dft_r2c_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                data, (mycomplex *)data, m_comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
        bwd = fftw_mpi_plan_dft_c2r_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                (mycomplex *)data, data, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
#else
        fftwf_plan_with_nthreads(desired_threads);
        fftwf_mpi_init();
        alloc_local = fftwf_mpi_local_size_3d_transposed(
                gsize_0[0], gsize_0[1], gsize_0[2]/2+1, m_comm,
                &src_local_n0, &src_local_0_start, &src_local_n1, &src_local_1_start);

        data = fftwf_alloc_real(2*alloc_local); // source terms
        fwd = fftwf_mpi_plan_dft_r2c_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                data, (mycomplex *)data, m_comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
        bwd = fftwf_mpi_plan_dft_c2r_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                (mycomplex *)data, data, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
#endif
        cout << "myN: " << myN[0] * myN[1] * myN[2] << endl;
        cout << "alloc_local: " << alloc_local << endl;
        cout << "2*alloc_local: " << 2*alloc_local << endl;
        cout << "local_n0: " << "rank " << m_rank << " " << src_local_n0 << endl;
        cout << "start_local_n0: " << src_local_0_start << endl;

        std::memset(data, 0, 2*alloc_local);   // ensure zeros

        m_buffer = new Real[myN[0] * myN[1] * myN[2]];

        // initialize MPI helper types (order iz,iy,ix, where z is the fastest
        // index).  These types will only be used by ranks < m_size/2
        if (m_size%2 != 0) { // this makes live easier
            cout << "PoissonSolverScalarFFTW_cyclicConvolution.h: ERROR: Require 1D domain decomposition and even number of processes.";
            abort();
        }
        int supersize[3] = { (int)src_local_n0, (int)gsize_0[1], 2*((int)gsize_0[2]/2+1) };
        int subsize[3] = { (int)myN[0], (int)myN[1], (int)myN[2] };
        int start_lo[3] = { 0, 0, 0 };
        int start_hi[3] = { (int)src_local_n0/2, 0, 0 };

        MPI_Type_create_subarray(3, supersize, subsize, start_lo, MPI_ORDER_C, MPIREAL, &m_arylo_t);
        MPI_Type_create_subarray(3, supersize, subsize, start_hi, MPI_ORDER_C, MPIREAL, &m_aryhi_t);
        MPI_Type_commit(&m_arylo_t);
        MPI_Type_commit(&m_aryhi_t);

        _initialize_green(); // setup fft of Green's function

        //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
        //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
        //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
        //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
    }

    ~PoissonSolverScalarFFTW_MPI()
    {
        delete m_buffer;
        MPI_Type_free(&m_arylo_t);
        MPI_Type_free(&m_aryhi_t);
        this->dispose();
    }

    void solve()
    {
        // Note: _cub2fftw() is called from outside via public member call (for
        // efficiency)

        _mpi_send_src();

// #ifndef _FLOAT_PRECISION_
//         fftw_execute(fwd);
// #else
//         fftwf_execute(fwd);
// #endif

        // _show_some();

//         _solve();

// #ifndef _FLOAT_PRECISION_
//         fftw_execute(bwd);
// #else
//         fftwf_execute(bwd);
// #endif

        _mpi_recv_sol();

        _fftw2cub(m_buffer);
    }

    void dispose()
    {
#ifndef _FLOAT_PRECISION_
        fftw_destroy_plan(fwd);
        fftw_destroy_plan(bwd);
        fftw_free(data);
        fftw_mpi_cleanup();
#else
        fftwf_destroy_plan(fwd);
        fftwf_destroy_plan(bwd);
        fftwf_free(data);
        fftwf_mpi_cleanup();
#endif
    }

    inline size_t _offset(const int blockID) const
    {
        const BlockInfo &info = local_infos[blockID];
        return _offset(info);
    }

    inline size_t _offset_ext(const BlockInfo &info) const
    {
        assert(local_infos[info.blockID].blockID == info.blockID);
        return _offset(local_infos[info.blockID]);
    }

    inline size_t _offset(const BlockInfo &info) const
    {
        // info must be a local BlockInfo! (obtained from
        // grid->getResidentBlocksInfo())
        const int myIstart[3] = {
            info.index[0]*bs[0],
            info.index[1]*bs[1],
            info.index[2]*bs[2]
        };
        return myIstart[2] + myN[2]*(myIstart[1] + myN[1]*myIstart[0]);
    }
    inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
    {
        return offset + z + myN[2]*(y + myN[1] * x);
    }
    inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
    {
        const size_t dest_index = _dest(offset, z, y, x);
        assert(dest_index>=0 && dest_index<myN[0]*myN[1]*myN[2]);
        m_buffer[dest_index] = rhs;
    }

    void _cub2fftw(Real * const out) const
    {
#pragma omp parallel for
        for(size_t i=0; i<local_infos.size(); ++i) {
            const BlockInfo info = local_infos[i];
            BlockType& b = *(BlockType*)info.ptrBlock;
            const size_t offset = _offset(info);

            for(int ix=0; ix<BlockType::sizeX; ix++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int iz=0; iz<BlockType::sizeZ; iz++) {
                        const size_t dest_index = _dest(offset, iz, iy, ix);
                        assert(dest_index>=0 && dest_index<myN[0]*myN[1]*myN[2]);
                        out[dest_index] = b(ix,iy,iz).p;
                        //TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
                    }
        }
    }
};

#undef MPIREAL

#endif /* POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV */
