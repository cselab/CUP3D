//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland. All rights reserved.
//  Distributed under the terms of the “CC BY-NC 3.0” license.
//  No commercial use allowed without permission.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

/*
 * Interface for using CubismUP_3D as a library.
 *
 * This code intentionally does not include any internal headers,
 * as they require predefined values as compilation flags.
 */

#ifndef CUBISMUP3D_H
#define CUBISMUP3D_H

#include <mpi.h>
#include <array>
#include <vector>

class Simulation;  // Internal class.

namespace cubismup3d {

/*
 * Callbacks and other information for `IF3D_ExternalObstacleOperator`.
 *
 * This structure enables the user to define a custom obstacle.
 */
struct ExternalObstacleSettings
{
    /*
     * Check if given box is touching (intersecting) the object.
     *
     * False positives are allowed.
     */
    bool (*is_touching_fn)(void *obj,
                           const double lo[3],
                           const double hi[3]) = nullptr;

    /*
     * Returns the signed distance to the object boundary.
     *
     * Positive values are to be returned for points inside the object,
     * negative for points outside of the object. Must be precise only close to
     * the obstacle surface.
     */
    double (*signed_distance_fn)(void *obj,
                                 double x,
                                 double y,
                                 double z) = nullptr;

    /*
     * Returns the local object velocity at the given location.
     */
    void (*velocity_fn)(void *obj,
                        double x,
                        double y,
                        double z,
                        double out[3]) = nullptr;

    /* Arbitrary pointer to be sent to the callbacks above. */
    void *obj = nullptr;

    /* Approx. length of the object. */
    double length = 0.0;

    /* Object center location. (center of mass? probably not so important) */
    std::array<double, 3> position{0.5, 0.5, 0.5};
};


/*
 * Parameters and settings of the simulation.
 */
struct SimulationConfig
{
    std::array<int, 3> nprocs{1, 1, 1};
    std::array<int, 3> bpd{1, 1, 1};
    bool dump2d = false;
    bool dump3d = false;
    double cfl = 0.1;
    double lambda = 1e5;
    std::array<double, 3> uinf{0.0, 0.0, 0.0};
    double nu;
    double tend;
    int visualization_freq = 0;       // After how many time steps to save?
    double visualization_time = 0.0;  // After how much time to save?
    std::vector<ExternalObstacleSettings> eos;
    std::string factory_content;      // Additional objects.
};


/*
 * Wrapper for the internal API.
 */
class SimulationWrapper
{
private:
    Simulation *S;
public:
    SimulationWrapper(MPI_Comm comm, const SimulationConfig &config);
    ~SimulationWrapper();

    /*
     * Get settings struct pointer attached to the `k`-th external obstacle.
     */
    ExternalObstacleSettings *get_eos(int k);

    /*
     * Get maximum time step `dt` ensuring stability.
     */
    double calc_max_timestep() const;

    /*
     * Perform one time step. It DOES NOT check if given `dt` is larger than the
     * current max allowed time step.
     *
     * Returns true if the simulation is finished.
     */
    bool timestep(double dt);

    /*
     * Get velocity values.
     *
     * Given points MUST belong to the current rank's subdomain.
     */
    void linear_interpolation(
            const double * __restrict__ x,
            const double * __restrict__ y,
            const double * __restrict__ z,
            double * __restrict__ vx,
            double * __restrict__ vy,
            double * __restrict__ vz,
            int N) const;
};

/*
 * Return block sizes (which is internally a compile-time constant).
 */
void get_block_sizes(int out[3]);


}  // cubismup3d

#endif
