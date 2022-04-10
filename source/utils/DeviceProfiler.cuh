#pragma once

#include <map>
#include <string>
#include <iostream>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <mpi.h>

#include <cuda_runtime.h>
#include "../include/helper_cuda.h"

// https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf?page=1&tab=scoredesc#tab-top
template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

class DeviceProfiler
{
public:
  DeviceProfiler(MPI_Comm m_comm) : m_comm_(m_comm)
  {
    MPI_Comm_rank(m_comm, &rank_);
  }

  ~DeviceProfiler() 
  {
    for (auto &[key, prof] : profs_)
    {
      checkCudaErrors(cudaEventDestroy(prof.start));
      checkCudaErrors(cudaEventDestroy(prof.stop));
    }
  }

  void startProfiler(std::string tag, cudaStream_t s)
  {
#ifdef BICGSTAB_PROFILER
    auto search_it = profs_.find(tag);
    if (search_it != profs_.end())
    {
      assert(!profs_[tag].started);
      profs_[tag].started = true;
      checkCudaErrors(cudaEventRecord(profs_[tag].start, s));
    }
    else 
    {
      ProfilerData dat = {true, 0.};
      profs_.emplace(tag, dat);
      checkCudaErrors(cudaEventCreate(&(profs_[tag].start)));
      checkCudaErrors(cudaEventCreate(&(profs_[tag].stop)));
      checkCudaErrors(cudaEventRecord(profs_[tag].start, s));
    }
#endif
  }

  void stopProfiler(std::string tag, cudaStream_t s)
  {
#ifdef BICGSTAB_PROFILER
    auto search_it = profs_.find(tag);
    assert(search_it != profs_.end());

    assert(profs_[tag].started);
    profs_[tag].started = false;
    checkCudaErrors(cudaEventRecord(profs_[tag].stop, s));
    checkCudaErrors(cudaEventSynchronize(profs_[tag].stop));

    float et = 0.;
    checkCudaErrors(cudaEventElapsedTime(&et, profs_[tag].start, profs_[tag].stop));
    profs_[tag].elapsed += et;
#endif
  }

  void print(std::string mt)
  {
#ifdef BICGSTAB_PROFILER
    std::string out;
		out += string_format("[DeviceProfiler] rank %i:\n", rank_);
    out += string_format("%10s:\t%.4e [ms]\n", mt.c_str(), profs_[mt].elapsed);
    for (auto &[key, prof] : profs_)
      if (key != mt)
        out += string_format("%10s:\t%.4e [ms]\t%6.2f%% of runtime\n", 
                             key.c_str(), prof.elapsed, profs_[key].elapsed / profs_[mt].elapsed * 100.);

    std::cout << out;
#endif
  }

protected:
  int rank_;
  MPI_Comm m_comm_;
  struct ProfilerData {
    bool started;
    float elapsed;
    cudaEvent_t start;
    cudaEvent_t stop;
  };
  std::map<std::string, ProfilerData> profs_;
};
