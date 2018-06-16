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

    void put_data();
    void get_data();

    void print_full(const std::string hint = "")
    {
        ostringstream fname;
        fname << "rank_" << m_rank << hint << ".dat";
        ofstream out(fname.str());
        for (int i = 0; i < m_NN0t; ++i)
        {
            for (int j = 0; j < m_local_NN1; ++j)
            {
                for (int k = 0; k < m_Nzhat; ++k)
                {
                    const int idx = k + 2*m_Nzhat*(j + m_local_NN1*i);
                    out << std::scientific << *(m_buf_full+idx) << ":" << std::scientific << *(m_buf_full+idx+1) << '\t';
                }
                out << std::endl;
            }
            out << std::endl;
            out << std::endl;
        }
    }

    void print_tp(const std::string hint = "", const double scale = 1.0)
    {
        ostringstream fname;
        fname << "rank_" << m_rank << hint << ".dat";
        ofstream out(fname.str());
        for (int i = 0; i < m_local_N0; ++i)
        {
            for (int j = 0; j < m_NN1; ++j)
            {
                for (int k = 0; k < 2*m_Nzhat; k+=2)
                {
                    const int idx = k + 2*m_Nzhat*(j + m_NN1*i);
                    out << std::scientific << *(m_buf_tp+idx) * scale << ":" << std::scientific << *(m_buf_tp+idx+1) * scale << '\t';
                }
                out << std::endl;
            }
            out << std::endl;
            out << std::endl;
        }
    }

    void initialize()
    {
        std::memset(m_buf_tp, 0, 2*m_tp_size*sizeof(Real));
        int d = m_local_N0*m_N1*m_N2*m_rank;
        for (int i = 0; i < m_local_N0; ++i)
        {
            for (int j = 0; j < m_N1; ++j)
            {
                for (int k = 0; k < m_N2; ++k)
                {
                    const int idx = k + 2*m_Nzhat*(j + m_NN1*i);
                    m_buf_tp[idx] = d++;
                }
            }
        }
    }

    void transform_fwd()
    {
        initialize();
        print_tp("initial");
        _FFTW_(execute)(m_fwd_2D);
        _FFTW_(execute)(m_fwd_tp);
        _copy_fwd_local();
        _FFTW_(execute)(m_fwd_1D);
        print_full("transformed");
    }

    void transform_bwd()
    {
        assert(m_initialized);
        _FFTW_(execute)(m_bwd_1D);
        _copy_bwd_local();
        _FFTW_(execute)(m_bwd_tp);
        _FFTW_(execute)(m_bwd_2D);
        // print_tp("final", 1.0/(m_NN0t*m_NN1t*m_NN2t));
        // print_tp("final");
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
        std::memset(m_buf_full, 0, 2*m_full_size*sizeof(Real));
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

        // // tmp
        // ostringstream fname;
        // fname << "rank_" << m_rank << "_Ghat.dat";
        // ofstream out(fname.str());
        // for (ptrdiff_t i = 0; i < m_NN0t; ++i)
        //     for (ptrdiff_t j = 0; j < m_local_NN1; ++j)
        //         for (ptrdiff_t k = 0; k < m_Nzhat; ++k)
        //         {
        //             const size_t linidx = k + m_Nzhat*(j + m_local_NN1*i);
        //             out << "(" << i << ", " << j+m_start_NN1 << ", " << k << "): " << std::scientific << G_hat[linidx][0] << '\t' << std::scientific << G_hat[linidx][1] << std::endl;
        //         }
        // out.close();

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
