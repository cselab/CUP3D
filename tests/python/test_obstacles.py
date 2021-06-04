from utils import TestCaseEx
import sys
import libcubismup3d as cup
import unittest

class TestObstacles(TestCaseEx):
    def test_sphere(self):
        argv = [
            '-nprocsx', '1',
            '-nprocsy', '1',
            '-nprocsz', '1',
            '-bpdx', '8',
            '-bpdy', '4',
            '-bpdz', '4',
            '-extentx', '64.0',
            '-extenty', '-1',
            '-extentz', '-1',
            '-nsteps', '5',
            '-nu', '1.23',
            '-factory-content', 'Sphere L=10.0 xpos=10.0 xvel=100.0 yvel=0.0 zvel=0.0 bForcedInSimFrame=1 bFixFrameOfRef=1\n',
        ]

        s = cup.Simulation(argv)
        s.add_obstacle(cup.Sphere(position=[40.0, 10.0, 10.0], radius=8.0))

        # What do I test here?
        s.run()

    def test_moving_sphere(self):
        argv = [
            '-nprocsx', '1',
            '-nprocsy', '1',
            '-nprocsz', '1',
            '-bpdx', '8',
            '-bpdy', '4',
            '-bpdz', '4',
            '-extentx', '64.0',
            '-extenty', '-1',
            '-extentz', '-1',
            '-nsteps', '5',
            '-nu', '1.23',
        ]

        def pos_vel(t):
            # x, y, z, vx, vy, vz
            return [32.0 + 1.0 * t, 16.0, 16.0, 1.0, 0.0, 0.0]

        s = cup.Simulation(argv)
        s.add_obstacle(cup.create_moving_sphere(s, func=pos_vel, radius=10.0, length=20.0))
        s.run()
