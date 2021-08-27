#include <cmath>
#include <vector>

template <typename T>
struct Vector3 {
    T &operator[](int k) { return x_[k]; }
    const T &operator[](int k) const { return x_[k]; }

    friend Vector3 operator+(const Vector3 &a, const Vector3 &b) {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    }
    friend Vector3 operator-(const Vector3 &a, const Vector3 &b) {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }

    friend Vector3 cross(const Vector3 &a, const Vector3 &b) {
        return Vector3{
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        };
    }

    friend auto dot(Vector3 a, Vector3 b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    T x_[3];
};


inline Real pointTriangleSqrDistance(
        const Vector3<Vector3<Real>> &tri,
        const Vector3<Real> p) {
    const auto n = cross(tri[1] - tri[0], tri[2] - tri[0]);
    const auto p0 = p - tri[0];
    const auto p1 = p - tri[1];
    const auto p2 = p - tri[2];
    const auto n0 = cross(p1, p2);
    const auto n1 = cross(p2, p0);
    const auto n2 = cross(p0, p1);
    const bool inside0 = dot(n0, n) >= 0;
    const bool inside1 = dot(n1, n) >= 0;
    const bool inside2 = dot(n2, n) >= 0;
    if (inside0 && inside1 && inside2) {
        auto dot_ = dot(p - tri[0], n);
        return dot_ * dot_ / dot(n, n);
    }
    Vector3<Real> pu;
    Vector3<Real> vu;
    if (!inside0){
        vu = tri[2] - tri[1];
        pu = p1;
    } else if (!inside1) {
        vu = tri[0] - tri[2];
        pu = p2;
    } else {
        vu = tri[1] - tri[0];
        pu = p0;
    }
    const auto vuvu = dot(vu, vu);
    const auto vupu = dot(vu, pu);
    if (vupu <= 0) {
        return dot(pu, pu);
    } else if (vupu >= vuvu) {
        const auto pv = pu - vu;
        return dot(pv, pv);
    } else {
        const auto out = dot(pu, pu) - vupu * vupu / vuvu;
        return out >= 0 ? out : 0;
    }
}


class Mesh {
public:
    Mesh(std::vector<Vector3<Real>> x,
         std::vector<Vector3<int>> tri) :
        x_{std::move(x)}, tri_{std::move(tri)}
    { }

    Real nonConvexSDF(Vector3<Real> p) const {
        // Find the closest triangle and the distance to it.
        Vector3<Vector3<Real>> closest{};
        Real minSqrDist = 1e100;
        for (int i = 0; i < (int)tri_.size(); ++i) {
            Vector3<Vector3<Real>> t{
                x_[tri_[i][0]],
                x_[tri_[i][1]],
                x_[tri_[i][2]],
            };
            const Real sqrDist = pointTriangleSqrDistance(t, p);
            if (sqrDist < minSqrDist) {
                minSqrDist = sqrDist;
                closest = t;
            }
        }
        const auto dist = std::sqrt(minSqrDist);

        // Check on which side of the closest triangle we are.
        const auto n = cross(closest[1] - closest[0], closest[2] - closest[0]);
        const auto side = dot(n, p - closest[0]);
        return std::copysign(dist, side);
    }

private:
    std::vector<Vector3<Real>> x_;
    std::vector<Vector3<int>> tri_;
};
