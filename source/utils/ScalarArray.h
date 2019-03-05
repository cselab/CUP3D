//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland. All rights reserved.
//  Distributed under the terms of the “CC BY-NC 3.0” license.
//  No commercial use allowed without permission.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

#ifndef CubismUP_3D_utils_ScalarArray_h
#define CubismUP_3D_utils_ScalarArray_h

#include <array>

namespace cubismup3d {

/*
 * Extension of std::array that implements operators +, - between themselves
 * and operator * with a scalar.
 */
template <typename T, int N>
struct ScalarArray {
    T v[N];

    T &operator[](const size_t k) { return v[k]; }
    const T &operator[](const size_t k) const { return v[k]; }
    size_t size(void) const { return N; }
    T *data(void) { return v; }
    const T *data(void) const { return v; }

    // This is not the optimal implementation (values get constructed and then
    // replaced), but compiler should to the job of optimizing unnecessary
    // parts away.
    friend inline ScalarArray operator+(const ScalarArray &A,
                                        const ScalarArray &B) {
        ScalarArray<T, N> result;
        for (int i = 0; i < N; ++i)
            result[i] = A[i] + B[i];
        return result;
    }

    friend inline ScalarArray operator-(const ScalarArray &A,
                                        const ScalarArray &B) {
        ScalarArray<T, N> result;
        for (int i = 0; i < N; ++i)
            result[i] = A[i] - B[i];
        return result;
    }

    friend inline ScalarArray operator*(const T &A, const ScalarArray &B) {
        ScalarArray<T, N> result;
        for (int i = 0; i < N; ++i)
            result[i] = A * B[i];
        return result;
    }

    friend inline ScalarArray operator*(const ScalarArray &A, const T &B) {
        ScalarArray<T, N> result;
        for (int i = 0; i < N; ++i)
            result[i] = A[i] * B;
        return result;
    }
};

}  // cubismup3d

#endif
