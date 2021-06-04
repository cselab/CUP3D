//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

#ifndef CUBISMUP_3D_LINEAR_INTERPOLATION_H
#define CUBISMUP_3D_LINEAR_INTERPOLATION_H

#include "Operator.h"
#include <cstdio>
#include <stdexcept>
#include <unordered_map>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

/* Distribute points over blocks. */
class InterpolationOperator : public Operator
{
  struct Particle {
    int id;
    double x[3];
  };
public:

  InterpolationOperator(SimulationData &s) :
    Operator(s),
    invCellSize_{
      s.grid->NX / s.grid->maxextent,
      s.grid->NY / s.grid->maxextent,
      s.grid->NZ / s.grid->maxextent,
    },
    maxLevel_{s.grid->levelMax}
  { }

  template <typename Array>
  void distributeParticlesToBlocks(Array &&points)
  {
    // Initialize mapping from (block index, level) to the block.
    const std::vector<cubism::BlockInfo> &vInfo = sim.vInfo();
    const int numBlocks = (int)vInfo.size();
    blocks_.reserve(numBlocks);
    for (int i = 0; i < numBlocks; ++i) {
      const auto &info = vInfo[i];
      const Key key{info.index[0], info.index[1], info.index[2], info.level};
      blocks_.insert(Map::value_type{key, i});
    }

    // Count how many particles to which block.
    blockOffsets_.resize(numBlocks + 1, 0);
    const int numParticles = (int)points.size();
#pragma omp parallel for
    for (int i = 0; i < numParticles; ++i) {
      const Particle p = _toParticle(points, i);
      int blockIdx = _getBlockIdx(p.x);
#pragma omp atomic
      ++blockOffsets_[blockIdx];
    }

    // Compute offsets from counts.
    for (int i = 0; i < numBlocks; ++i)
      blockOffsets_[i + 1] += blockOffsets_[i];
    assert(blockOffsets_.back() == numParticles);

    // Rerun to distribute particles.
    particles_.resize(numParticles);
#pragma omp parallel for
    for (int i = 0; i < numParticles; ++i) {
      const Particle p = _toParticle(points, i);
      int blockIdx = _getBlockIdx(p.x);

      int pos;
#pragma omp atomic capture
      pos = --blockOffsets_[blockIdx];

      particles_[pos] = p;
    }
  }

  using Operator::compute;

  std::pair<const Particle *, const Particle *> getBlockParticles(int blockIdx) const {
    assert(0 <= blockIdx && blockIdx < blockOffsets_.size() - 1);
    return {
      particles_.data() + blockOffsets_[blockIdx],
      particles_.data() + blockOffsets_[blockIdx + 1],
    };
  }

  inline int infoToBlockIdx(const cubism::BlockInfo &info) const noexcept {
    // compute() could send us this information.
    const Key key{info.index[0], info.index[1], info.index[2], info.level};
    auto it = blocks_.find(key);
    assert(it != blocks_.end());
    return it->second;
  }

  // Operator stuff we don't care about.
  void operator()(const double /* dt */) override { abort(); }
  std::string getName(void) { return "LinearCellCenteredInterpolation"; }

private:
  std::array<int, 3> _getGlobalCellIndex(const double x[3]) const
  {
    return {
      (int)(x[0] * invCellSize_[0]),
      (int)(x[1] * invCellSize_[1]),
      (int)(x[2] * invCellSize_[2]),
    };
  }

  /// Get index of the block that contains the given point.
  /// If there is no such (local) block, throw an error.
  int _getBlockIdx(const double x[3]) const
  {
    const auto cellIdx = _getGlobalCellIndex(x);
    const int bx = cellIdx[0] / (int)FluidBlock::sizeX;
    const int by = cellIdx[1] / (int)FluidBlock::sizeY;
    const int bz = cellIdx[2] / (int)FluidBlock::sizeZ;

    for (int level = 0; level <= maxLevel_; ++level) {
      const Key key{bx >> level, by >> level, bz >> level, level};
      const auto it = blocks_.find(key);
      if (it != blocks_.end())
        return it->second;
    }
    fprintf(stderr, "Point (%g %g %g) not in local blocks\n", x[0], x[1], x[2]);
    throw std::runtime_error("point not in local blocks");
  }

  template <typename Array>
  static Particle _toParticle(Array &&points, int i)
  {
    const auto out = points[i];
    Particle p;
    p.id = i;
    p.x[0] = out[0];
    p.x[1] = out[1];
    p.x[2] = out[2];
    return p;
  }

  const double invCellSize_[3];
  const int maxLevel_;

  struct Hash
  {
    inline size_t operator()(std::array<int, 4> a) const noexcept
    {
      // https://stackoverflow.com/questions/42701688/using-an-unordered-map-with-arrays-as-keys
      size_t h = 0;
      for (int i = 0; i < 4; ++i)
        h ^= std::hash<int>{}(a[i])  + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  std::vector<int> blockOffsets_;
  using Key = std::array<int, 4>;
  using Map = std::unordered_map<Key, int, Hash>;
  Map blocks_;
  std::vector<Particle> particles_;
};


/* Interpolate values for points contained in a single block. */
template <typename Getter, typename Setter>
class LinearCellCenteredInterpolationKernel
{
public:
  static constexpr std::array<int, 3> stencil_start{-1, -1, -1};
  static constexpr std::array<int, 3> stencil_end{2, 2, 2};
  const cubism::StencilInfo stencil;

  LinearCellCenteredInterpolationKernel(
      const InterpolationOperator &interpolator,
      Getter &getter,
      Setter &setter,
      std::vector<int> components) :
    stencil{-1, -1, -1, 2, 2, 2, true, std::move(components)},
    interpolator_{interpolator},
    getter_{getter},
    setter_{setter}
  { }

  template <typename Lab, typename BlockType>
  void operator()(Lab &lab,
                  const cubism::BlockInfo &info,
                  BlockType &o) const
  {
    typedef typename FluidGridMPI::BlockType Block;
    const int blockIdx = interpolator_.infoToBlockIdx(info);

    const auto span = interpolator_.getBlockParticles(blockIdx);
    const int level = info.level;
    const double invh = 1.0 / info.h_gridpoint;
    for (const auto *p = span.first; p != span.second; ++p) {

      // Position within a block, where the coordinate (0.0, 0.0) refers to the
      // center of the cell (0, 0).
      const double scaledPos[3] = {
        // -0.5 because we are interpolating cell-centered values.
        invh * (p->x[0] - info.origin[0]) - 0.5,
        invh * (p->x[1] - info.origin[1]) - 0.5,
        invh * (p->x[2] - info.origin[2]) - 0.5,
      };

      const int idx[3] = {
        // Due to rounding errors we have to make sure we don't go out of bounds.
        // (int)std::floor(scaledPos[0]),
        // (int)std::floor(scaledPos[1]),
        // (int)std::floor(scaledPos[2]),
        std::max(-1, std::min(FluidBlock::sizeX - 1, (int)std::floor(scaledPos[0]))),
        std::max(-1, std::min(FluidBlock::sizeY - 1, (int)std::floor(scaledPos[1]))),
        std::max(-1, std::min(FluidBlock::sizeZ - 1, (int)std::floor(scaledPos[2]))),
      };

      // Compute 1D weights.
      const double w[3] = {
        scaledPos[0] - idx[0],
        scaledPos[1] - idx[1],
        scaledPos[2] - idx[2],
      };

      // Do M2P interpolation.
      const double w000 = (1 - w[0]) * (1 - w[1]) * (1 - w[2]);
      const double w010 = (1 - w[0]) * (    w[1]) * (1 - w[2]);
      const double w100 = (    w[0]) * (1 - w[1]) * (1 - w[2]);
      const double w110 = (    w[0]) * (    w[1]) * (1 - w[2]);
      const double w001 = (1 - w[0]) * (1 - w[1]) * (    w[2]);
      const double w011 = (1 - w[0]) * (    w[1]) * (    w[2]);
      const double w101 = (    w[0]) * (1 - w[1]) * (    w[2]);
      const double w111 = (    w[0]) * (    w[1]) * (    w[2]);
      setter_(p->id,
              w000 * getter_(lab.read(idx[0]    , idx[1]    , idx[2]    ))
            + w010 * getter_(lab.read(idx[0]    , idx[1] + 1, idx[2]    ))
            + w100 * getter_(lab.read(idx[0] + 1, idx[1]    , idx[2]    ))
            + w110 * getter_(lab.read(idx[0] + 1, idx[1] + 1, idx[2]    ))
            + w001 * getter_(lab.read(idx[0]    , idx[1]    , idx[2] + 1))
            + w011 * getter_(lab.read(idx[0]    , idx[1] + 1, idx[2] + 1))
            + w101 * getter_(lab.read(idx[0] + 1, idx[1]    , idx[2] + 1))
            + w111 * getter_(lab.read(idx[0] + 1, idx[1] + 1, idx[2] + 1)));
    }
  }

private:
  const InterpolationOperator &interpolator_;
  Getter &getter_;
  Setter &setter_;
};



/*
 * Cell-centered mesh to particle linear interpolation.
 *
 * For each given point, interpolate the value of the field.
 *
 * Arguments:
 *   - points - Array of the points, where points support operator [] for
 *              accessing x, y, z coordinates.
 *   - getter - Lambda of a single argument (FluidElement), returning the value
 *              to be interpolated. The value should support addition and
 *              multiplication by a scalar. If you are interpolating more than
 *              one value simultaneously, check `utils/ScalarArray.h`.
 *   - setter - Lambda function of two arguments (point ID, interpolated value).
 *   - components - Stencil components.
 */
template <typename Array, typename Getter, typename Setter>
void linearCellCenteredInterpolation(
    SimulationData &sim,
    Array &&points,
    Getter &&getter,
    Setter &&setter,
    std::vector<int> components)
{
  InterpolationOperator interpolator{sim};
  interpolator.distributeParticlesToBlocks(std::forward<Array>(points));

  using Kernel = LinearCellCenteredInterpolationKernel<
      typename std::remove_reference<Getter>::type,
      typename std::remove_reference<Setter>::type>;
  Kernel kernel{interpolator, getter, setter, std::move(components)};
  interpolator.compute(kernel);
}


CubismUP_3D_NAMESPACE_END

#endif
