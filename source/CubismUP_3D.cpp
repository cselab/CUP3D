//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

// DEPRECATED?

/*
 * Please check CubismUP_3D.h for information about the functions below.
 */

#include "CubismUP_3D.h"
#include "Simulation.h"
#include "obstacles/IF3D_ExternalObstacleOperator.h"
#include "obstacles/IF3D_LinearInterpolation.h"
#include "utils/ScalarArray.h"


namespace cubismup3d {

void get_block_sizes(int out[3])
{
    out[0] = FluidBlock::sizeX;
    out[1] = FluidBlock::sizeY;
    out[2] = FluidBlock::sizeZ;
}

SimulationWrapper::SimulationWrapper(MPI_Comm comm, const SimulationConfig &C)
{
    // In order not to change the internals, here we construct the factory and
    // other "command-line" arguments.
    std::string factory_content;
    for (const auto &eos : C.eos) {
        assert(eos.length > 0 && "For some reason, length > 0 is required.");
        factory_content += "IF3D_ExternalObstacleOperator"
                           " L=" + std::to_string(eos.length)
                         + " xpos=" + std::to_string(eos.position[0])
                         + " ypos=" + std::to_string(eos.position[1])
                         + " zpos=" + std::to_string(eos.position[2])
                         + '\n';
    }
    factory_content += C.factory_content;

    assert((!C.visualization_freq && !C.visualization_time) || (C.dump2d || C.dump3d));

    // Cannot really be const due to the annoying `char **` below.
    std::string _argv[] = {
        "-nprocsx", std::to_string(C.nprocs[0]),
        "-nprocsy", std::to_string(C.nprocs[1]),
        "-nprocsz", std::to_string(C.nprocs[2]),
        "-bpdx", std::to_string(C.bpd[0]),
        "-bpdy", std::to_string(C.bpd[1]),
        "-bpdz", std::to_string(C.bpd[2]),
        "-2Ddump", std::to_string((int)C.dump2d),
        "-3Ddump", std::to_string((int)C.dump3d),
        "-cfl", std::to_string(C.cfl),
        "-lambda", std::to_string(C.lambda),
        "-uinfx", std::to_string(C.uinf[0]),
        "-uinfy", std::to_string(C.uinf[1]),
        "-uinfz", std::to_string(C.uinf[2]),
        "-fdump", std::to_string(C.visualization_freq),
        "-tdump", std::to_string(C.visualization_time),
        "-nu", std::to_string(C.nu),
        "-tend", std::to_string(C.tend),
        "-factory-content", factory_content,
    };
    constexpr int argc = (int)(sizeof(_argv) / sizeof(_argv[0]));
    char *argv[argc];
    for (int i = 0; i < argc; ++i)
        argv[i] = const_cast<char *>(_argv[i].c_str());

    for (int i = 0; i < argc; i += 2)
        printf("==> %s = %s\n", _argv[i].c_str(), _argv[i + 1].c_str());

    ArgumentParser parser(argc, argv);
    S = new Simulation(comm, parser);
    S->rampup = false;

    const auto &obstacles = S->getObstacleVector();
    for (int i = 0; i < (int)C.eos.size(); ++i) {
        IF3D_ExternalObstacleOperator * const obstacle =
                static_cast<IF3D_ExternalObstacleOperator *>(obstacles[i]);
        obstacle->settings = C.eos[i];
    }
}

SimulationWrapper::~SimulationWrapper()
{
    delete S;
    S = nullptr;
}

ExternalObstacleSettings *SimulationWrapper::get_eos(int k) {
    const auto &obstacles = S->getObstacleVector();
    IF3D_ExternalObstacleOperator * const obstacle =
            static_cast<IF3D_ExternalObstacleOperator *>(obstacles[k]);
    return &obstacle->settings;
}

double SimulationWrapper::calc_max_timestep() const
{
    return S->calcMaxTimestep();
}

bool SimulationWrapper::timestep(const double dt)
{
    return S->timestep(dt);
}

void SimulationWrapper::linear_interpolation(
        const double * const __restrict__ x,
        const double * const __restrict__ y,
        const double * const __restrict__ z,
        double * const __restrict__ vx,
        double * const __restrict__ vy,
        double * const __restrict__ vz,
        const int N) const
{

    struct ArrayWrapper
    {
        const double * const __restrict__ x;
        const double * const __restrict__ y;
        const double * const __restrict__ z;
        const int N;

        std::array<double, 3> operator[](const int i) const
        {
            return {x[i], y[i], z[i]};
        }

        int size(void) const { return N; }
    };

    LinearInterpolation interpolation(S->grid);
    interpolation.interpolate(
            ArrayWrapper{x, y, z, N},
            // Getter (what to interpolate):
            [](const FluidElement &e) {
                return cubismup3d::ScalarArray<double, 3>{e.u, e.v, e.w};
            },
            // Setter (where to store the result):
            [vx, vy, vz](const int i,
                         const cubismup3d::ScalarArray<double, 3> &v) {
                vx[i] = v[0];
                vy[i] = v[1];
                vz[i] = v[2];
            },
            // FluidElement components for the stencil:
            {1, 2, 3});
}

}  // cubismup3d
