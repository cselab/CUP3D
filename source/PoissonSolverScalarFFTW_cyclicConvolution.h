//
//  CubismUP_3D
//
//  Written by Fabian Wermelinger.
//  Copyright (c) 2018 ETHZ. All rights reserved.
//
//  This algorithm uses the cyclic convolution method described in Eastwood and
//  Brownrigg (1979) for unbounded domains.
//  WARNING: This implementation only works with a 1D domain decomposition
//  along the x-coordinate.
#ifndef POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV
#define POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV

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
    myplan fwd, bwd;

    // MPI
    MPI_Comm m_comm;
    int m_rank, m_size;
    Real* m_buffer; // MPI send/recv buffer
    Real* m_kernel; // Green's function

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
    const size_t nz_hat = gsize_0[2]/2+1;
    const double norm_factor = 1./(gsize_0[0]*gsize_0[1]*gsize_0[2]);
    const double h = grid.getBlocksInfo().front().h_gridpoint;
    ptrdiff_t alloc_local=0;
    ptrdiff_t src_local_n0=0, src_local_0_start=0, src_local_n1=0, src_local_1_start=0;

    // setup Green's function and compute transform, done once
    void _initialize_green()
    {
        myplan greenPlan;
        ptrdiff_t alloc;
        ptrdiff_t local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;
        Real* tf_buf;

#ifndef _FLOAT_PRECISION_
        alloc = fftw_mpi_local_size_3d_transposed(
                gsize_0[0], gsize_0[1], nz_hat, m_comm,
                &local_n0, &local_0_start, &local_n1, &local_1_start);

        tf_buf    = fftw_alloc_real(2*alloc);
        greenPlan = fftw_mpi_plan_dft_r2c_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                tf_buf, (mycomplex *)tf_buf, m_comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
#else
        alloc = fftwf_mpi_local_size_3d_transposed(
                gsize_0[0], gsize_0[1], nz_hat, m_comm,
                &local_n0, &local_0_start, &local_n1, &local_1_start);

        tf_buf    = fftwf_alloc_real(2*alloc);
        greenPlan = fftwf_mpi_plan_dft_r2c_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                tf_buf, (mycomplex *)tf_buf, m_comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
#endif
        assert(local_n0 == src_local_n0);
        assert(local_n1 == src_local_n1);
        assert(local_0_start == src_local_0_start);
        assert(local_1_start == src_local_1_start);

        _set_kernel(tf_buf);

#ifndef _FLOAT_PRECISION_
        fftw_execute(greenPlan);
#else
        fftwf_execute(greenPlan);
#endif

        m_kernel = new Real[alloc]; // FFT for this kernel is real

        const mycomplex *const G_hat = (mycomplex *) tf_buf;
#pragma omp parallel for
        for(long j = 0; j<static_cast<long>(src_local_n1); ++j)
            for(long i = 0; i<static_cast<long>(gsize_0[0]); ++i)
                for(long k = 0; k<static_cast<long>(nz_hat); ++k)
                {
                    const size_t linidx = (j*gsize_0[0] +i)*nz_hat + k;
                    m_kernel[linidx]  = G_hat[linidx][0]; // need real part only
                }

#ifndef _FLOAT_PRECISION_
        fftw_free(tf_buf);
        fftw_destroy_plan(greenPlan);
#else
        fftwf_free(tf_buf);
        fftwf_destroy_plan(greenPlan);
#endif
    }

    void _set_kernel(Real* const kern)
    {
        // This algorithm requires m_size >= 2.

        // This factor is due to the discretization of the convolution
        // integtal.  It is composed of (h*h*h) * (-1/[4*pi*h]), where h is the
        // uniform grid spacing.  The first factor is the discrete volume
        // element of the convolution integral; the second factor belongs to
        // Green's function on a uniform mesh.
        const double fac = -h*h / (4.0*M_PI);

        // octant 000
        if (m_rank < m_size/2)
        {
            for (size_t i = 0; i < (size_t)src_local_n0; ++i)
            {
                const double xi = src_local_0_start + i;
                for (size_t j = 0; j < myN[1]; ++j)
                {
                    const double yi = j;
                    for (size_t k = 0; k < myN[2]; ++k)
                    {
                        const size_t idx = k + 2*nz_hat*(j + gsize_0[1]*i);
                        const double zi = k;
                        const double r = std::sqrt(xi*xi + yi*yi + zi*zi);
                        if (r > 0.0)
                            kern[idx] = fac/r;
                        else
                            kern[idx] = fac;
                    }
                }
            }
        }

        // mirror
        // octant 100
        Real* const buf = new Real[src_local_n0 * myN[1] * myN[2]];
        if (m_rank < m_size/2)
        {
            const size_t off = (m_rank==0) ? 1 : 0;
            for (size_t i = off; i < (size_t)src_local_n0; ++i)
                for (size_t j = 0; j < myN[1]; ++j)
                    for (size_t k = 0; k < myN[2]; ++k)
                    {
                        const size_t src_idx = k + 2*nz_hat*(j + gsize_0[1]*i);
                        const size_t dst_idx = k + myN[2]*(j + myN[1]*(src_local_n0-i-1));
                        buf[dst_idx] = kern[src_idx]; // flip along x-dimension
                    }
            int dst_rank = m_size - m_rank - 1;
            int volume = (src_local_n0-off)*myN[1]*myN[2];
            MPI_Send(buf, volume, MPIREAL, dst_rank, m_rank, m_comm);
        }
        else
        {
            int src_rank = m_size - m_rank - 1;
            int volume = src_local_n0*myN[1]*myN[2];
            MPI_Recv(buf, volume, MPIREAL, src_rank, src_rank, m_comm, MPI_STATUSES_IGNORE);

            for (size_t i = 0; i < (size_t)src_local_n0; ++i)
                for (size_t j = 0; j < myN[1]; ++j)
                    for (size_t k = 0; k < myN[2]; ++k)
                    {
                        const size_t dst_idx = k + 2*nz_hat*(j + gsize_0[1]*i);
                        const size_t src_idx = k + myN[2]*(j + myN[1]*i);
                        kern[dst_idx] = buf[src_idx];
                    }
        }
        delete[] buf;

        // rotations
        // octant 010/110
        for (size_t i = 0; i < (size_t)src_local_n0; ++i)
            for (size_t j = 1; j < myN[1]; ++j)
                for (size_t k = 0; k < myN[2]; ++k)
                {
                    const size_t src_idx = k + 2*nz_hat*((myN[1]-j) + gsize_0[1]*i);
                    const size_t dst_idx = k + 2*nz_hat*((myN[1]+j-1) + gsize_0[1]*i);
                    kern[dst_idx] = kern[src_idx];
                }

        // octant 011/111
        for (size_t i = 0; i < (size_t)src_local_n0; ++i)
            for (size_t j = 1; j < myN[1]; ++j)
                for (size_t k = 1; k < myN[2]; ++k)
                {
                    const size_t src_idx = myN[2]-k + 2*nz_hat*((myN[1]-j) + gsize_0[1]*i);
                    const size_t dst_idx = myN[2]+k-1 + 2*nz_hat*((myN[1]+j-1) + gsize_0[1]*i);
                    kern[dst_idx] = kern[src_idx];
                }

        // octant 001/101
        for (size_t i = 0; i < (size_t)src_local_n0; ++i)
            for (size_t j = 0; j < myN[1]; ++j)
                for (size_t k = 1; k < myN[2]; ++k)
                {
                    const size_t src_idx = myN[2]-k + 2*nz_hat*(j + gsize_0[1]*i);
                    const size_t dst_idx = myN[2]+k-1 + 2*nz_hat*(j + gsize_0[1]*i);
                    kern[dst_idx] = kern[src_idx];
                }
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
        // distribute the source data to the zero-padded arrays.  This
        // algorithm requires m_size >= 2.
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
        // receive back solution.  This algorithm requires m_size >= 2.
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
        mycomplex* const rho_hat = (mycomplex *) data;
        const Real* const G_hat  = m_kernel;

        // perform convolution in frequency domain and normalize
#pragma omp parallel for
        for(long j = 0; j<static_cast<long>(src_local_n1); ++j)
            for(long i = 0; i<static_cast<long>(gsize_0[0]); ++i)
                for(long k = 0; k<static_cast<long>(nz_hat); ++k)
                {
                    const size_t linidx = (j*gsize_0[0] +i)*nz_hat + k;
                    rho_hat[linidx][0] *= G_hat[linidx] * norm_factor;
                    rho_hat[linidx][1] *= G_hat[linidx] * norm_factor;
                }
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

        if (m_size%2 != 0 || m_size < 2) { // this makes live easier
            cout << "PoissonSolverScalarFFTW_cyclicConvolution.h: ERROR: Require 1D domain decomposition and even number of processes with at least 2 processes." << endl;
            abort();
        }

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
                gsize_0[0], gsize_0[1], nz_hat, m_comm,
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
                gsize_0[0], gsize_0[1], nz_hat, m_comm,
                &src_local_n0, &src_local_0_start, &src_local_n1, &src_local_1_start);

        data = fftwf_alloc_real(2*alloc_local); // source terms
        fwd = fftwf_mpi_plan_dft_r2c_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                data, (mycomplex *)data, m_comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
        bwd = fftwf_mpi_plan_dft_c2r_3d(gsize_0[0], gsize_0[1], gsize_0[2],
                (mycomplex *)data, data, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
#endif

        // initialize MPI helper types (order iz,iy,ix, where z is the fastest
        // index).  These types will only be used by ranks < m_size/2
        int supersize[3] = { (int)src_local_n0, (int)gsize_0[1], 2*((int)nz_hat) };
        int subsize[3] = { (int)myN[0], (int)myN[1], (int)myN[2] };
        int start_lo[3] = { 0, 0, 0 };
        int start_hi[3] = { (int)src_local_n0/2, 0, 0 };

        MPI_Type_create_subarray(3, supersize, subsize, start_lo, MPI_ORDER_C, MPIREAL, &m_arylo_t);
        MPI_Type_create_subarray(3, supersize, subsize, start_hi, MPI_ORDER_C, MPIREAL, &m_aryhi_t);
        MPI_Type_commit(&m_arylo_t);
        MPI_Type_commit(&m_aryhi_t);

        _initialize_green(); // setup fft of Green's function (constant)
        m_buffer = new Real[myN[0] * myN[1] * myN[2]];
    }

    ~PoissonSolverScalarFFTW_MPI()
    {
        delete[] m_kernel;
        delete[] m_buffer;
        MPI_Type_free(&m_arylo_t);
        MPI_Type_free(&m_aryhi_t);
        this->dispose();
    }

    void solve()
    {
        // Note: _cub2fftw() is called from outside via public member call (for
        // efficiency)
        std::memset(data, 0, sizeof(Real)*2*alloc_local); // ensure zero's for multiple invocations of solve()
        _mpi_send_src();

#ifndef _FLOAT_PRECISION_
        fftw_execute(fwd);
#else
        fftwf_execute(fwd);
#endif

        _solve();

#ifndef _FLOAT_PRECISION_
        fftw_execute(bwd);
#else
        fftwf_execute(bwd);
#endif

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
