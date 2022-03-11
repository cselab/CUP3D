#include "Utils.h"
#include <mpi.h>
#include <stdexcept>

namespace cubismup3d {

void initMPI(int *argc, char ***argv)
{
  int provided;
  #ifdef CUP_ASYNC_DUMP
    const auto SECURITY = MPI_THREAD_MULTIPLE;
  #else
    const auto SECURITY = MPI_THREAD_FUNNELED;
  #endif
  MPI_Init_thread(argc, argv, SECURITY, &provided);
  if (provided < SECURITY) {
    fprintf(stderr, "ERROR: MPI implementation does not have required thread support\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void finalizeMPI(void)
{
  MPI_Finalize();
}

std::string computeNumBlocksArg(int cells, int blockSize)
{
  if (cells % blockSize != 0)
    throw std::runtime_error("num cells not divisible by block size");
  return std::to_string(cells / blockSize);
}

}  // cubismup3d
