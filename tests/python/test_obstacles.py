from utils import TestCaseEx
import sys
import libcubismup3d as cup3d
import unittest

class TestObstacles(TestCaseEx):
    def test_sphere(self):
        argv = [
            '-bpdx', '8',
            '-bpdy', '4',
            '-bpdz', '4',
            '-extentx', '64.0',
            '-extenty', '-1',
            '-extentz', '-1',
            '-nsteps', '3',
            '-nu', '1.23',
            '-factory-content', 'Sphere L=10.0 xpos=10.0 xvel=100.0 yvel=0.0 zvel=0.0 bForcedInSimFrame=1 bFixFrameOfRef=1\n',
            '-levelStart', '0',
            '-levelMax', '2',
        ]

        s = cup3d.Simulation(argv)
        s.add_obstacle(cup3d.Sphere(position=[40.0, 10.0, 10.0], radius=8.0))

        # For now test only that nothing crashes.
        s.run()
