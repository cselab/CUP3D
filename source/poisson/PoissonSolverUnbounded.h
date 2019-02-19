//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Fabian Wermelinger.
//
//  This algorithm uses the cyclic convolution method described in Eastwood and
//  Brownrigg (1979) for unbounded domains.
//  WARNING: This implementation only works with a 1D domain decomposition
//  along the x-coordinate.
#ifndef POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV
#define POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV

#include "PoissonSolver.h"

class PoissonSolverUnbounded : public PoissonSolver
{
  // domain dimensions
  const size_t m_N0 = gsize[0];  // Nx-points of original domain
  const size_t m_N1 = gsize[1];  // Ny-points of original domain
  const size_t m_N2 = gsize[2];  // Nz-points of original domain
  const size_t m_NN0 = 2*m_N0; // Nx-points of padded domain
  const size_t m_NN1 = 2*m_N1; // Ny-points of padded domain
  const size_t m_NN2 = 2*m_N2; // Nz-points of padded domain
  const size_t m_NN0t = 2*m_N0-1; //Nx of padded domain (w/o periodic copies, actual transform size)
  const size_t m_NN1t = 2*m_N1-1; //Ny of padded domain (w/o periodic copies, actual transform size)
  const size_t m_NN2t = 2*m_N2-1; //Nz of padded domain (w/o periodic copies, actual transform size)

  const size_t m_local_N0 = m_N0 / m_size;
  const size_t m_start_N0 = m_local_N0 * m_rank;
  const size_t m_local_NN0 = m_NN0/m_size, m_start_NN0 = m_local_NN0*m_rank;
  const size_t m_local_NN1 = m_NN1/m_size, m_start_NN1 = m_local_NN1*m_rank;

  // data buffers for input and transform.  Split into 2 buffers to exploit
  // smaller transpose matrix and fewer FFT's due to zero-padded domain.
  // This is at the cost of higher memory requirements.
  const size_t m_Nzhat = m_NN2t/2 + 1; // for symmetry in r2c transform
  const size_t m_tp_size   = m_local_N0  * m_NN1  * m_Nzhat;
  const size_t m_full_size = m_NN0t * m_local_NN1 * m_Nzhat;
  // FFT normalization factor
  const Real m_norm_factor = 1.0 / (m_NN0t*h * m_NN1t*h * m_NN2t*h);
  Real* data;   // input, output, transpose and 2D FFTs (m_local_N0 x m_NN1 x 2m_Nzhat)
  Real* m_buf_full; // full block of m_NN0t x m_local_NN1 x 2m_Nzhat for 1D FFTs
  Real* m_kernel;   // FFT of Green's function (real part, m_NN0t x m_local_NN1 x m_Nzhat)

  // FFTW plans
  myplan m_fwd_1D;
  myplan m_bwd_1D;
  myplan m_fwd_2D;
  myplan m_bwd_2D;
  myplan m_fwd_tp; // use FFTW's transpose facility
  myplan m_bwd_tp; // use FFTW's transpose facility

 public:
  PoissonSolverUnbounded(SimulationData&s) : PoissonSolver(s)
  {
    stridez = 2*m_Nzhat;
    stridey = m_NN1;

    if (m_N0 % m_size != 0 || m_NN1 % m_size != 0) {
      std::cout << "PoissonSolverUnbounded: ERROR: Number of cells N0 and 2*N1 must be evenly divisible by the number of processes." << std::endl;
      abort();
    }

    int supported_threads;
    MPI_Query_thread(&supported_threads);
    if (supported_threads < MPI_THREAD_FUNNELED) {
    std::cout << "PoissonSolverUnbounded: ERROR: MPI implementation does not support threads." << std::endl;
    abort();
    }

    const int retval = _FFTW_(init_threads)();
    if (retval == 0) {
      std::cout << "PoissonSolverUnbounded: ERROR: Call to fftw_init_threads() returned zero." << std::endl;
      abort();
    }
    const int desired_threads = omp_get_max_threads();
    _FFTW_(plan_with_nthreads)(desired_threads);
    _FFTW_(mpi_init)();

    // FFTW plans
    data   = _FFTW_(alloc_real)( 2*m_tp_size );
    m_buf_full = _FFTW_(alloc_real)( 2*m_full_size );

    // 1D plan
    {
      const int n[1] = {static_cast<int>(m_NN0t)};
      const int howmany = static_cast<int>( m_local_NN1 * m_Nzhat );
      const int stride  = static_cast<int>( m_local_NN1 * m_Nzhat );
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
      const int n[2] = {static_cast<int>(m_NN1t), static_cast<int>(m_NN2t)};
      const int howmany = static_cast<int>(m_local_N0);
      const int stride = 1;
      const int rembed[2] = {static_cast<int>(m_NN1), static_cast<int>(2*m_Nzhat)}; // unit: sizeof(Real)
      const int cembed[2] = {static_cast<int>(m_NN1), static_cast<int>(m_Nzhat)};   // unit: sizeof(mycomplex)
      const int rdist = static_cast<int>( m_NN1 * 2*m_Nzhat ); // unit: sizeof(Real)
      const int cdist = static_cast<int>( m_NN1 * m_Nzhat ); // unit: sizeof(mycomplex)
      m_fwd_2D = _FFTW_(plan_many_dft_r2c)(2, n, howmany,
              data, rembed, stride, rdist,
              (mycomplex*)data, cembed, stride, cdist,
              FFTW_MEASURE);
      m_bwd_2D = _FFTW_(plan_many_dft_c2r)(2, n, howmany,
              (mycomplex*)data, cembed, stride, cdist,
              data, rembed, stride, rdist,
              FFTW_MEASURE);
    }

    // transpose plan
    m_fwd_tp = _FFTW_(mpi_plan_many_transpose)(m_N0, m_NN1, 2*m_Nzhat,
            m_local_N0, m_local_NN1,
            data, data,
            m_comm, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);

    m_bwd_tp = _FFTW_(mpi_plan_many_transpose)(m_NN1, m_N0, 2*m_Nzhat,
            m_local_NN1, m_local_N0,
            data, data,
            m_comm, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_IN);

    _initialize_green();

    clear();
  }

  PoissonSolverUnbounded(const PoissonSolverUnbounded& c) = delete;
  ~PoissonSolverUnbounded()
  {
    _FFTW_(free)(data);
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

  void solve() override
  {
    // Note: _cub2fftw() is called from outside via public member call (for
    // efficiency)
    //_cub2fftw();

    _FFTW_(execute)(m_fwd_2D);
    _FFTW_(execute)(m_fwd_tp);
    _copy_fwd_local();
    _FFTW_(execute)(m_fwd_1D);

    {
      mycomplex* const rho_hat = (mycomplex*)m_buf_full;
      const Real* const G_hat = m_kernel;
      #pragma omp parallel for
      for (size_t i = 0; i < m_NN0t; ++i)
      for (size_t j = 0; j < m_local_NN1; ++j)
      for (size_t k = 0; k < m_Nzhat; ++k)
      {
        const size_t idx = k + m_Nzhat*(j + m_local_NN1*i);
        rho_hat[idx][0] *= G_hat[idx]; //normalization is carried on in G_hat
        rho_hat[idx][1] *= G_hat[idx]; //normalization is carried on in G_hat
      }
    }

    _FFTW_(execute)(m_bwd_1D);
    _copy_bwd_local();
    _FFTW_(execute)(m_bwd_tp);
    _FFTW_(execute)(m_bwd_2D);

    _fftw2cub();

    clear(); // clear buffers for next invocation
  }

 private:

  void _initialize_green()
  {
    myplan green1D;
    myplan green2D;
    myplan greenTP;

    const size_t tf_size = m_local_NN0 * m_NN1 * m_Nzhat;
    Real* tf_buf = _FFTW_(alloc_real)( 2*tf_size );

    // 1D plan
    {
    const int n[1] = {static_cast<int>(m_NN0t)};
    const int howmany = static_cast<int>( m_local_NN1 * m_Nzhat );
    const int stride  = static_cast<int>( m_local_NN1 * m_Nzhat );
    const int* embed = n;
    const int dist = 1;
    green1D = _FFTW_(plan_many_dft)(1, n, howmany,
            (mycomplex*)tf_buf, embed, stride, dist,
            (mycomplex*)tf_buf, embed, stride, dist,
            FFTW_FORWARD, FFTW_MEASURE);
    }

    // 2D plan
    {
    const int n[2] = {static_cast<int>(m_NN1t), static_cast<int>(m_NN2t)};
    const int howmany = static_cast<int>(m_local_NN0);
    const int stride = 1;
    const int rembed[2] = {static_cast<int>(m_NN1), static_cast<int>(2*m_Nzhat)}; // unit: sizeof(Real)
    const int cembed[2] = {static_cast<int>(m_NN1), static_cast<int>(m_Nzhat)};   // unit: sizeof(mycomplex)
    const int rdist = static_cast<int>( m_NN1 * 2*m_Nzhat );                      // unit: sizeof(Real)
    const int cdist = static_cast<int>( m_NN1 * m_Nzhat );                        // unit: sizeof(mycomplex)
    green2D = _FFTW_(plan_many_dft_r2c)(2, n, howmany,
            tf_buf, rembed, stride, rdist,
            (mycomplex*)tf_buf, cembed, stride, cdist,
            FFTW_MEASURE);
    }

    greenTP = _FFTW_(mpi_plan_many_transpose)(m_NN0, m_NN1, 2*m_Nzhat,
            m_local_NN0, m_local_NN1,
            tf_buf, tf_buf,
            m_comm, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);

    // This factor is due to the discretization of the convolution
    // integtal.  It is composed of (h*h*h) * (-1/[4*pi*h]), where h is the
    // uniform grid spacing.  The first factor is the discrete volume
    // element of the convolution integral; the second factor belongs to
    // Green's function on a uniform mesh.
    for (size_t i = 0; i < m_local_NN0; ++i)
    for (size_t j = 0; j < m_N1; ++j)
    for (size_t k = 0; k < m_N2; ++k)
    {
        const size_t I = m_start_NN0 + i;
        const Real xi = I>=m_N0? 2*m_N0-1 - I : I;
        const Real yi = j>=m_N1? 2*m_N1-1 - j : j;
        const Real zi = k>=m_N2? 2*m_N2-1 - k : k;
        const double r = std::sqrt(xi*xi + yi*yi + zi*zi);
        const size_t idx = k + 2*m_Nzhat*(j + m_NN1*i);
        if (r > 0) tf_buf[idx] = - h * h / (4*M_PI*r);
        else tf_buf[idx] = - Real(0.1924173658) * h * h;
    }

    _FFTW_(execute)(green2D);
    _FFTW_(execute)(greenTP);
    _FFTW_(execute)(green1D);

    const size_t kern_size = m_NN0t * m_local_NN1 * m_Nzhat;
    m_kernel = _FFTW_(alloc_real)(kern_size); // FFT for this kernel is real
    std::memset(m_kernel, 0, kern_size*sizeof(Real));

    const mycomplex *const G_hat = (mycomplex *) tf_buf;
    #pragma omp parallel for
    for (size_t i = 0; i < m_NN0t; ++i)
    for (size_t j = 0; j < m_local_NN1; ++j)
    for (size_t k = 0; k < m_Nzhat; ++k)
    {
      const size_t linidx = k + m_Nzhat*(j + m_local_NN1*i);
      m_kernel[linidx] = G_hat[linidx][0] * m_norm_factor;// need real part only
    }

    _FFTW_(free)(tf_buf);
    _FFTW_(destroy_plan)(green1D);
    _FFTW_(destroy_plan)(green2D);
    _FFTW_(destroy_plan)(greenTP);
  }

  void clear()
  {
    std::memset(data, 0, 2*m_tp_size*sizeof(Real));
    std::memset(m_buf_full, 0, 2*m_full_size*sizeof(Real));
  }

  void _copy_fwd_local()
  {
    #pragma omp parallel for
    for (size_t i = 0; i < m_N0; ++i)
    for (size_t j = 0; j < m_local_NN1; ++j) {
      const Real* const src = data + 2*m_Nzhat*(j + m_local_NN1*i);
      Real* const dst = m_buf_full + 2*m_Nzhat*(j + m_local_NN1*i);
      std::memcpy(dst, src, 2*m_Nzhat*sizeof(Real));
    }
  }

  void _copy_bwd_local()
  {
    #pragma omp parallel for
    for (size_t i = 0; i < m_N0; ++i)
    for (size_t j = 0; j < m_local_NN1; ++j) {
      const Real* const src = m_buf_full + 2*m_Nzhat*(j + m_local_NN1*i);
      Real* const dst = data + 2*m_Nzhat*(j + m_local_NN1*i);
      std::memcpy(dst, src, 2*m_Nzhat*sizeof(Real));
    }
  }
};

#undef MPIREAL

#endif /* POISSONSOLVERSCALARFFTW_CYCLICCONVOLUTION_H_NUOUYWFV */
