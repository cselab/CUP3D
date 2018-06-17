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

#include <vector>
#include <cassert>
#include <cstring>
#include <fftw3-mpi.h>

#include "Definitions.h"

#ifndef _FLOAT_PRECISION_
#define _FFTW_(s) fftw_##s
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#define MPIREAL MPI_DOUBLE
#else
#define _FFTW_(s) fftwf_##s
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#define MPIREAL MPI_FLOAT
#endif /* _FLOAT_PRECISION_ */

using namespace std;

#include <BlockInfo.h>

class My3DFFT_Infinite_MPI
{
public:
    My3DFFT_Infinite_MPI() : m_initialized(false) {}
    My3DFFT_Infinite_MPI(const My3DFFT_Infinite_MPI& c) = delete;

    My3DFFT_Infinite_MPI(const size_t N0, const size_t N1, const size_t N2, const Real h, const int desired_threads, const MPI_Comm comm)
    {
        _initialize(N0, N1, N2, h, desired_threads, comm);
    }

    ~My3DFFT_Infinite_MPI()
    {
        _FFTW_(free)(m_buf_tp);
        _FFTW_(free)(m_buf_full);
        _FFTW_(free)(m_kernel);
        _FFTW_(destroy_plan)(m_fwd_1D);
        _FFTW_(destroy_plan)(m_bwd_1D);
        _FFTW_(destroy_plan)(m_fwd_2D);
        _FFTW_(destroy_plan)(m_bwd_2D);
        _FFTW_(destroy_plan)(m_fwd_tp);
        _FFTW_(destroy_plan)(m_bwd_tp);
        _FFTW_(mpi_cleanup)();
    }

    // interface
    inline void setup(const size_t N0, const size_t N1, const size_t N2, const Real h, const int desired_threads, const MPI_Comm comm)
    {
        _initialize(N0, N1, N2, h, desired_threads, comm);
    }

    inline void put_data(const size_t ix, const size_t iy, const size_t iz, const Real v)
    {
        const size_t idx = this->get_data_linear_offset(ix,iy,iz);
        this->put_data(idx, v);
    }

    inline void put_data(const size_t idx, const Real v)
    {
        assert(idx < m_local_N0*m_NN1*2*m_Nzhat);
        m_buf_tp[idx] = v;
    }

    inline void get_data(const size_t ix, const size_t iy, const size_t iz, Real& v) const
    {
        const size_t idx = this->get_data_linear_offset(ix,iy,iz);
        this->get_data(idx, v);
    }

    void get_data(const size_t idx, Real& v) const
    {
        assert(idx < m_local_N0*m_NN1*2*m_Nzhat);
        v = m_buf_tp[idx];
    }

    inline size_t get_data_linear_offset(const size_t ix, const size_t iy, const size_t iz) const
    {
        return iz + 2*static_cast<size_t>(m_Nzhat)*(iy + static_cast<size_t>(m_NN1)*ix);
    }

    inline void clear()
    {
        std::memset(m_buf_tp, 0, 2*m_tp_size*sizeof(Real));
        std::memset(m_buf_full, 0, 2*m_full_size*sizeof(Real));
    }

    void transform_fwd()
    {
        assert(m_initialized);
        _FFTW_(execute)(m_fwd_2D);
        _FFTW_(execute)(m_fwd_tp);
        _copy_fwd_local();
        _FFTW_(execute)(m_fwd_1D);
    }

    void transform_bwd()
    {
        assert(m_initialized);
        _FFTW_(execute)(m_bwd_1D);
        _copy_bwd_local();
        _FFTW_(execute)(m_bwd_tp);
        _FFTW_(execute)(m_bwd_2D);
    }

    void convolve()
    {
        assert(m_initialized);

        mycomplex* const rho_hat = (mycomplex*)m_buf_full;
        const Real* const G_hat = m_kernel;
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < m_NN0t; ++i)
            for (ptrdiff_t j = 0; j < m_local_NN1; ++j)
                for (ptrdiff_t k = 0; k < m_Nzhat; ++k)
                {
                    const size_t linidx = k + m_Nzhat*(j + m_local_NN1*i);
                    rho_hat[linidx][0] *= G_hat[linidx] * m_norm_factor;
                    rho_hat[linidx][1] *= G_hat[linidx] * m_norm_factor;
                }
    }

    inline void sweep()
    {
        assert(m_initialized);
        this->transform_fwd();
        this->convolve();
        this->transform_bwd();
    }

private:
    bool m_initialized;

    // domain dimensions
    ptrdiff_t m_N0;  // Nx-points of original domain
    ptrdiff_t m_N1;  // Ny-points of original domain
    ptrdiff_t m_N2;  // Nz-points of original domain
    ptrdiff_t m_NN0; // Nx-points of padded domain
    ptrdiff_t m_NN1; // Ny-points of padded domain
    ptrdiff_t m_NN2; // Nz-points of padded domain
    ptrdiff_t m_NN0t; // Nx-points of padded domain (w/o periodic copies, actual transform size)
    ptrdiff_t m_NN1t; // Ny-points of padded domain (w/o periodic copies, actual transform size)
    ptrdiff_t m_NN2t; // Nz-points of padded domain (w/o periodic copies, actual transform size)
    Real m_h; // uniform grid spacing

    // MPI related
    MPI_Comm m_comm;
    int m_rank, m_size;
    ptrdiff_t m_local_N0, m_start_N0;
    ptrdiff_t m_local_NN0, m_start_NN0, m_local_NN1, m_start_NN1;

    // data buffers for input and transform.  Split into 2 buffers to exploit
    // smaller transpose matrix and fewer FFT's due to zero-padded domain.
    // This is at the cost of higher memory requirements.
    ptrdiff_t m_Nzhat; // for r2c transform
    ptrdiff_t m_tp_size;
    ptrdiff_t m_full_size;
    Real m_norm_factor; // FFT normalization factor
    Real* m_buf_tp;     // input, output, transpose and 2D FFTs (m_local_N0 x m_NN1 x 2m_Nzhat)
    Real* m_buf_full;   // full block of m_NN0t x m_local_NN1 x 2m_Nzhat for 1D FFTs
    Real* m_kernel;     // FFT of Green's function (real part, m_NN0t x m_local_NN1 x m_Nzhat)

    // FFTW plans
    myplan m_fwd_1D;
    myplan m_bwd_1D;
    myplan m_fwd_2D;
    myplan m_bwd_2D;
    myplan m_fwd_tp; // use FFTW's transpose facility
    myplan m_bwd_tp; // use FFTW's transpose facility

    // helpers
    void _initialize(const size_t N0, const size_t N1, const size_t N2, const Real h, const int desired_threads, const MPI_Comm comm)
    {
        // Setup MPI
        m_comm = comm;
        MPI_Comm_rank(m_comm, &m_rank);
        MPI_Comm_size(m_comm, &m_size);

        // initialize domain dimensions
        m_N0   = N0;       // Number of cells in simulation domain
        m_N1   = N1;
        m_N2   = N2;
        m_NN0  = 2*N0;     // Zero-padded domain for infinite domain
        m_NN1  = 2*N1;
        m_NN2  = 2*N2;
        m_NN0t = 2*N0 - 1; // Size of FFT (removing periodic copies)
        m_NN1t = 2*N1 - 1;
        m_NN2t = 2*N2 - 1;
        m_h    = h;        // uniform grid spacing

        // some checks
        _check_init(desired_threads);

        // x-decomposition of input data (number of 2D slices)
        m_local_N0 = m_N0 / m_size;
        m_start_N0 = m_local_N0 * m_rank;

        // x- and y-decompostion of extended domain (zero-padded) required for
        // transpose
        m_local_NN0 = m_NN0 / m_size;
        m_local_NN1 = m_NN1 / m_size;
        m_start_NN0 = m_local_NN0 * m_rank;
        m_start_NN1 = m_local_NN1 * m_rank;

        // buffer sizes
        m_Nzhat     = m_NN2t/2 + 1; // for symmetry in r2c transform
        m_tp_size   = m_local_N0  * m_NN1  * m_Nzhat;
        m_full_size = m_NN0t * m_local_NN1 * m_Nzhat;

        // FFT normalization factor
        m_norm_factor = 1.0 / (m_NN0t * m_NN1t * m_NN2t);

        // FFTW plans
        m_buf_tp   = _FFTW_(alloc_real)( 2*m_tp_size );
        m_buf_full = _FFTW_(alloc_real)( 2*m_full_size );

        // 1D plan
        {
            const int n[1] = {m_NN0t};
            const int howmany = m_local_NN1 * m_Nzhat;
            const int stride = m_local_NN1 * m_Nzhat;
            const int* embed = n;
            const int dist = 1;
            m_fwd_1D = _FFTW_(plan_many_dft)(1, n, howmany,
                    (mycomplex*)m_buf_full, embed, stride, dist,
                    (mycomplex*)m_buf_full, embed, stride, dist,
                    FFTW_FORWARD, FFTW_MEASURE);
            m_bwd_1D = _FFTW_(plan_many_dft)(1, n, howmany,
                    (mycomplex*)m_buf_full, embed, stride, dist,
                    (mycomplex*)m_buf_full, embed, stride, dist,
                    FFTW_BACKWARD, FFTW_MEASURE);
        }

        // 2D plan
        {
            const int n[2] = {m_NN1t, m_NN2t};
            const int howmany = m_local_N0;
            const int stride = 1;
            const int rembed[2] = {m_NN1, 2*m_Nzhat}; // unit: sizeof(Real)
            const int cembed[2] = {m_NN1, m_Nzhat};   // unit: sizeof(mycomplex)
            const int rdist = m_NN1 * 2*m_Nzhat;      // unit: sizeof(Real)
            const int cdist = m_NN1 * m_Nzhat;        // unit: sizeof(mycomplex)
            m_fwd_2D = _FFTW_(plan_many_dft_r2c)(2, n, howmany,
                    m_buf_tp, rembed, stride, rdist,
                    (mycomplex*)m_buf_tp, cembed, stride, cdist,
                    FFTW_MEASURE);
            m_bwd_2D = _FFTW_(plan_many_dft_c2r)(2, n, howmany,
                    (mycomplex*)m_buf_tp, cembed, stride, cdist,
                    m_buf_tp, rembed, stride, rdist,
                    FFTW_MEASURE);
        }

        // transpose plan
        m_fwd_tp = _FFTW_(mpi_plan_many_transpose)(m_N0, m_NN1, 2*m_Nzhat,
                m_local_N0, m_local_NN1,
                m_buf_tp, m_buf_tp,
                m_comm, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);

        m_bwd_tp = _FFTW_(mpi_plan_many_transpose)(m_NN1, m_N0, 2*m_Nzhat,
                m_local_NN1, m_local_N0,
                m_buf_tp, m_buf_tp,
                m_comm, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_IN);

        // compute Green's function
        _initialize_green();

        // to be sure
        this->clear();

        // cool
        m_initialized = true;
    }

    void _check_init(const int desired_threads)
    {
        if (m_N0 % m_size != 0 || m_NN1 % m_size != 0)
        {
            std::cout << "PoissonSolverScalarFFTW_cyclicConvolution.h: ERROR: Number of cells N0 and 2*N1 must be evenly divisible by the number of processes." << std::endl;
            abort();
        }

        int supported_threads;
        MPI_Query_thread(&supported_threads);
        if (supported_threads < MPI_THREAD_FUNNELED)
        {
            std::cout << "PoissonSolverScalarFFTW_cyclicConvolution.h: ERROR: MPI implementation does not support threads." << std::endl;
            abort();
        }

        const int retval = _FFTW_(init_threads)();
        if (retval == 0)
        {
            std::cout << "PoissonSolverScalarFFTW_cyclicConvolution.h: ERROR: Call to fftw_init_threads() returned zero." << std::endl;
            abort();
        }

        _FFTW_(plan_with_nthreads)(desired_threads);
        _FFTW_(mpi_init)();
    }

    void _copy_fwd_local()
    {
#pragma omp parallel for
        for (int i = 0; i < m_N0; ++i)
            for (int j = 0; j < m_local_NN1; ++j)
            {
                const Real* const src = m_buf_tp + 2*m_Nzhat*(j + m_local_NN1*i);
                Real* const dst = m_buf_full + 2*m_Nzhat*(j + m_local_NN1*i);
                std::memcpy(dst, src, 2*m_Nzhat*sizeof(Real));
            }
    }

    void _copy_bwd_local()
    {
#pragma omp parallel for
        for (int i = 0; i < m_N0; ++i)
            for (int j = 0; j < m_local_NN1; ++j)
            {
                const Real* const src = m_buf_full + 2*m_Nzhat*(j + m_local_NN1*i);
                Real* const dst = m_buf_tp + 2*m_Nzhat*(j + m_local_NN1*i);
                std::memcpy(dst, src, 2*m_Nzhat*sizeof(Real));
            }
    }

    void _initialize_green()
    {
        myplan green1D;
        myplan green2D;
        myplan greenTP;

        const ptrdiff_t tf_size = m_local_NN0 * m_NN1 * m_Nzhat;
        Real* tf_buf = _FFTW_(alloc_real)( 2*tf_size );

        // 1D plan
        {
            const int n[1] = {m_NN0t};
            const int howmany = m_local_NN1 * m_Nzhat;
            const int stride = m_local_NN1 * m_Nzhat;
            const int* embed = n;
            const int dist = 1;
            green1D = _FFTW_(plan_many_dft)(1, n, howmany,
                    (mycomplex*)tf_buf, embed, stride, dist,
                    (mycomplex*)tf_buf, embed, stride, dist,
                    FFTW_FORWARD, FFTW_MEASURE);
        }

        // 2D plan
        {
            const int n[2] = {m_NN1t, m_NN2t};
            const int howmany = m_local_NN0;
            const int stride = 1;
            const int rembed[2] = {m_NN1, 2*m_Nzhat}; // unit: sizeof(Real)
            const int cembed[2] = {m_NN1, m_Nzhat};   // unit: sizeof(mycomplex)
            const int rdist = m_NN1 * 2*m_Nzhat;      // unit: sizeof(Real)
            const int cdist = m_NN1 * m_Nzhat;        // unit: sizeof(mycomplex)
            green2D = _FFTW_(plan_many_dft_r2c)(2, n, howmany,
                    tf_buf, rembed, stride, rdist,
                    (mycomplex*)tf_buf, cembed, stride, cdist,
                    FFTW_MEASURE);
        }

        greenTP = _FFTW_(mpi_plan_many_transpose)(m_NN0, m_NN1, 2*m_Nzhat,
                m_local_NN0, m_local_NN1,
                tf_buf, tf_buf,
                m_comm, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);


        _set_kernel(tf_buf);

        _FFTW_(execute)(green2D);
        _FFTW_(execute)(greenTP);
        _FFTW_(execute)(green1D);

        const ptrdiff_t kern_size = m_NN0t * m_local_NN1 * m_Nzhat;
        m_kernel = _FFTW_(alloc_real)(kern_size); // FFT for this kernel is real
        std::memset(m_kernel, 0, kern_size*sizeof(Real));

        const mycomplex *const G_hat = (mycomplex *) tf_buf;
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < m_NN0t; ++i)
            for (ptrdiff_t j = 0; j < m_local_NN1; ++j)
                for (ptrdiff_t k = 0; k < m_Nzhat; ++k)
                {
                    const size_t linidx = k + m_Nzhat*(j + m_local_NN1*i);
                    m_kernel[linidx]  = G_hat[linidx][0]; // need real part only
                }

        _FFTW_(free)(tf_buf);
        _FFTW_(destroy_plan)(green1D);
        _FFTW_(destroy_plan)(green2D);
        _FFTW_(destroy_plan)(greenTP);
    }

    void _set_kernel(Real* const kern)
    {
        // This algorithm requires m_size >= 2 and m_N0 % m_size == 0

        // This factor is due to the discretization of the convolution
        // integtal.  It is composed of (h*h*h) * (-1/[4*pi*h]), where h is the
        // uniform grid spacing.  The first factor is the discrete volume
        // element of the convolution integral; the second factor belongs to
        // Green's function on a uniform mesh.
        const double fac = -m_h*m_h / (4.0*M_PI);

        // octant 000
        if (m_rank < m_size/2)
        {
            for (ptrdiff_t i = 0; i < m_local_NN0; ++i)
            {
                const double xi = m_start_NN0 + i;
                for (ptrdiff_t j = 0; j < m_N1; ++j)
                {
                    const double yi = j;
                    for (ptrdiff_t k = 0; k < m_N2; ++k)
                    {
                        const double zi = k;
                        const double r = std::sqrt(xi*xi + yi*yi + zi*zi);
                        const ptrdiff_t idx = k + 2*m_Nzhat*(j + m_NN1*i);
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
        const int vol = m_local_NN0 * m_N1 * m_N2;
        Real* const buf = new Real[vol];
        if (m_rank < m_size/2)
        {
            for (ptrdiff_t i = 0; i < m_local_NN0; ++i)
                for (ptrdiff_t j = 0; j < m_N1; ++j)
                    for (ptrdiff_t k = 0; k < m_N2; ++k)
                    {
                        const ptrdiff_t src = k + 2*m_Nzhat*(j + m_NN1*i);
                        const ptrdiff_t dst = k + m_N2*(j + m_N1*(m_local_NN0-1-i));
                        buf[dst] = kern[src]; // flip along x-dimension
                    }
            int dst_rank = m_size - m_rank - 1;
            MPI_Send(buf, vol, MPIREAL, dst_rank, m_rank, m_comm);
        }
        else
        {
            int src_rank = m_size - m_rank - 1;
            MPI_Recv(buf, vol, MPIREAL, src_rank, src_rank, m_comm, MPI_STATUSES_IGNORE);
            for (ptrdiff_t i = 0; i < m_local_NN0; ++i)
                for (ptrdiff_t j = 0; j < m_N1; ++j)
                    for (ptrdiff_t k = 0; k < m_N2; ++k)
                    {
                        const ptrdiff_t dst = k + 2*m_Nzhat*(j + m_NN1*i);
                        const ptrdiff_t src = k + m_N2*(j + m_N1*i);
                        kern[dst] = buf[src];
                    }
        }
        delete[] buf;

        // rotations
        // octant 010/110
        for (ptrdiff_t i = 0; i < m_local_NN0; ++i)
            for (ptrdiff_t j = 0; j < m_N1; ++j)
                for (ptrdiff_t k = 0; k < m_N2; ++k)
                {
                    const ptrdiff_t src = k + 2*m_Nzhat*((m_N1-1-j) + m_NN1*i);
                    const ptrdiff_t dst = k + 2*m_Nzhat*((m_N1+j) + m_NN1*i);
                    kern[dst] = kern[src];
                }

        // octant 011/111
        for (ptrdiff_t i = 0; i < m_local_NN0; ++i)
            for (ptrdiff_t j = 0; j < m_N1; ++j)
                for (ptrdiff_t k = 0; k < m_N2; ++k)
                {
                    const ptrdiff_t src = m_N2-1-k + 2*m_Nzhat*((m_N1-1-j) + m_NN1*i);
                    const ptrdiff_t dst = m_N2+k + 2*m_Nzhat*((m_N1+j) + m_NN1*i);
                    kern[dst] = kern[src];
                }

        // octant 001/101
        for (ptrdiff_t i = 0; i < m_local_NN0; ++i)
            for (ptrdiff_t j = 0; j < m_N1; ++j)
                for (ptrdiff_t k = 0; k < m_N2; ++k)
                {
                    const ptrdiff_t src = m_N2-1-k + 2*m_Nzhat*(j + m_NN1*i);
                    const ptrdiff_t dst = m_N2+k + 2*m_Nzhat*(j + m_NN1*i);
                    kern[dst] = kern[src];
                }
    }
};


template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_MPI
{
public:
    PoissonSolverScalarFFTW_MPI(const int desired_threads, TGrid& g) :
        m_grid(g),
        m_local_infos(g.getResidentBlocksInfo())
    {
        const ptrdiff_t N0  = m_grid.getBlocksPerDimension(0) * BlockType::sizeX;
        const ptrdiff_t N1  = m_grid.getBlocksPerDimension(1) * BlockType::sizeY;
        const ptrdiff_t N2  = m_grid.getBlocksPerDimension(2) * BlockType::sizeZ;
        const Real h        = m_grid.getBlocksInfo().front().h_gridpoint;
        const MPI_Comm comm = m_grid.getCartComm();
        m_fft.setup(N0, N1, N2, h, desired_threads, comm);
    }

    PoissonSolverScalarFFTW_MPI(const PoissonSolverScalarFFTW_MPI& c) = delete;
    ~PoissonSolverScalarFFTW_MPI() = default;

    void solve()
    {
        // Note: _cub2fftw() is called from outside via public member call (for
        // efficiency)

        m_fft.sweep();

        _fftw2cub();

        m_fft.clear(); // clear buffers for next invocation
    }

    inline size_t _offset(const int blockID) const
    {
        const BlockInfo &info = m_local_infos[blockID];
        return _offset(info);
    }

    inline size_t _offset_ext(const BlockInfo &info) const
    {
        assert(m_local_infos[info.blockID].blockID == info.blockID);
        return _offset(m_local_infos[info.blockID]);
    }

    inline size_t _offset(const BlockInfo &info) const
    {
        // info must be a local BlockInfo! (obtained from
        // grid->getResidentBlocksInfo())
        const int myIstart[3] = {
            info.index[0] * BlockType::sizeX,
            info.index[1] * BlockType::sizeY,
            info.index[2] * BlockType::sizeZ
        };
        return m_fft.get_data_linear_offset(myIstart[0], myIstart[1], myIstart[2]);
    }

    inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
    {
        return offset + m_fft.get_data_linear_offset(x, y, z);
    }

    inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs)
    {
        const size_t dest_index = _dest(offset, z, y, x);
        m_fft.put_data(dest_index, rhs);
    }

    void _cub2fftw()
    {
#pragma omp parallel for
        for(size_t i=0; i<m_local_infos.size(); ++i)
        {
            const BlockInfo info = m_local_infos[i];
            BlockType& b = *(BlockType*)info.ptrBlock;
            const size_t offset = _offset(info);

            for(int ix=0; ix<BlockType::sizeX; ix++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int iz=0; iz<BlockType::sizeZ; iz++)
                    {
                        const size_t dest_index = _dest(offset, iz, iy, ix);
                        m_fft.put_data(dest_index, b(ix,iy,iz).p);
                    }
        }
    }

private:
    typedef typename TGrid::BlockType BlockType;

    // computational grid
    TGrid& m_grid;

    // local grid blocks
    const std::vector<BlockInfo> m_local_infos;

    // FFT solver
    My3DFFT_Infinite_MPI m_fft;

protected:

    void _fftw2cub()
    {
#pragma omp parallel for
        for(size_t i=0; i<m_local_infos.size(); ++i) {
            const BlockInfo info = m_local_infos[i];
            BlockType& b = *(BlockType*)info.ptrBlock;
            const size_t offset = _offset(info);

            for(int ix=0; ix<BlockType::sizeX; ix++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int iz=0; iz<BlockType::sizeZ; iz++) {
                        const size_t src_index = _dest(offset, iz, iy, ix);
                        m_fft.get_data(src_index, b(ix,iy,iz).p);
                    }
        }
    }
};

#undef _FFTW_
#undef MPIREAL

#endif /* POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV */
