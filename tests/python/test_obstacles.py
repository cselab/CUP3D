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
            '-bpdx', '64',
            '-bpdy', '32',
            '-bpdz', '32',
            '-extentx', '64.0',
            '-extenty', '-1',
            '-extentz', '-1',
            '-nsteps', '5',
            '-nu', '1.23',
            '-factory-content', 'Sphere L=10.0 xpos=10.0 xvel=100.0 yvel=0.0 zvel=0.0 bForcedInSimFrame=1 bFixFrameOfRef=1\n',
        ]

        S = cup.Simulation(argv)
        S.add_obstacle(cup.Sphere(position=[40.0, 10.0, 10.0], radius=8.0))

        # What do I test here?
        S.run()
