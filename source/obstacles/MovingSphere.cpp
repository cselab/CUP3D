#include "MovingSphere.h"
#include "ObstacleFromShape.h"

CubismUP_3D_NAMESPACE_BEGIN

class MovingSphereShape {
public:
  MovingSphereShape(MovingSphereFunc func, Real r) :
      posVelFunc_{std::move(func)}, r_{r}
  {
    setTime((Real)0.0);
  }

  bool isTouching(std::array<Real, 3> low, std::array<Real, 3> high) const {
    return !(low[0] > center_[0] + r_ || high[0] < center_[0] - r_
          || low[1] > center_[1] + r_ || high[1] < center_[1] - r_
          || low[2] > center_[2] + r_ || high[2] < center_[2] - r_);
  }

  Real signedDistance(std::array<Real, 3> position) const {
    const Real dx = position[0] - center_[0];
    const Real dy = position[1] - center_[1];
    const Real dz = position[2] - center_[2];
    return r_ - std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  std::array<Real, 3> comVelocity() const {
    return comVelocity_;
  }

  std::array<Real, 3> localRelativeVelocity(std::array<Real, 3> /* x */) const {
    return {0.0, 0.0, 0.0};  // The sphere is not rotating or anything.
  }

  Real lambdaFactor(Real /* time */) const noexcept {
    return (Real)1.0;
  }

  void setTime(Real t) {
    const auto out = posVelFunc_(t);
    center_[0] = out[0];
    center_[1] = out[1];
    center_[2] = out[2];
    comVelocity_[0] = out[3];
    comVelocity_[1] = out[4];
    comVelocity_[2] = out[5];
  }

private:
  MovingSphereFunc posVelFunc_;
  const Real r_;  // Radius.
  std::array<Real, 3> center_;
  std::array<Real, 3> comVelocity_;
};

class MovingSphere : public ObstacleFromShape<MovingSphereShape> {
public:
  using ObstacleFromShape<MovingSphereShape>::ObstacleFromShape;
};

std::shared_ptr<cubismup3d::Obstacle> createMovingSphere(
    SimulationData &s,
    const ObstacleArguments &args,
    MovingSphereFunc func,
    Real radius) {
  return std::make_shared<MovingSphere>(
        s, args, MovingSphereShape{std::move(func), radius});
}

CubismUP_3D_NAMESPACE_END
