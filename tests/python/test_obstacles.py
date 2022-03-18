from utils import TestCaseEx
import libcubismup3d as cup3d

import numpy as np

class TestObstacles(TestCaseEx):
    def test_sphere_chi_moments(self):
        """Test sphere chi field. Test loading the chi field."""
        extent = 10.0
        h = extent / 128
        r = 1.5
        xpos, ypos, zpos = (3.125, 2.5, 2.5)  # Pick multiples of h.
        argv = [
            '-bpdx', '8',
            '-bpdy', '4',
            '-bpdz', '4',
            '-extentx', str(extent),
            '-extenty', '-1',
            '-extentz', '-1',
            '-nsteps', '1',
            '-nu', '1.23',
            '-factory-content', f'Sphere L={2 * r} xpos={xpos} ypos={ypos} zpos={zpos} xvel=10.0 yvel=0.0 zvel=0.0 bForcedInSimFrame=1 bFixFrameOfRef=1\n',
            '-levelStart', '0',
            '-levelMax', '2',
        ]

        s = cup3d.Simulation(argv)
        s.run()

        with self.assertRaises(TypeError):
            out = np.empty((64, 64, 128), dtype=np.float32)
            s.fields.chi.to_uniform(out)  # wrong dtype

        with self.assertRaises(TypeError):
            out = np.empty((10, 12, 14), dtype=np.float64)
            s.fields.chi.to_uniform(out)  # wrong size

        with self.assertRaises(TypeError):
            out = np.empty((65, 65, 129), dtype=np.float64)[1:, 1:, 1:]
            s.fields.chi.to_uniform(out)  # not contiguous

        chi = np.empty((64, 64, 128), dtype=np.float64)
        chi2 = s.fields.chi.to_uniform(chi)
        self.assertEqual(chi.__array_interface__, chi2.__array_interface__)

        x = (0.5 + np.arange(128)) * h - xpos
        y = (0.5 + np.arange(64)) * h - ypos
        z = (0.5 + np.arange(64)) * h - zpos
        x = x[None, None, :]
        y = y[None, :, None]
        z = z[:, None, None]

        # Test moments.
        h3 = h ** 3
        M = chi.sum() * h3
        Mx = (x * chi).sum() * h3
        My = (y * chi).sum() * h3
        Mz = (z * chi).sum() * h3
        Mxx = (x * x * chi).sum() * h3
        Myy = (y * y * chi).sum() * h3
        Mzz = (z * z * chi).sum() * h3
        Mxy = (x * y * chi).sum() * h3
        Myz = (y * z * chi).sum() * h3
        Mzx = (z * x * chi).sum() * h3
        self.assertAlmostEqual(M, 4 / 3 * np.pi * r ** 3, places=3)
        self.assertLess(abs(Mx), 1e-9)
        self.assertLess(abs(My), 1e-9)
        self.assertLess(abs(Mz), 1e-9)
        self.assertAlmostEqual(Mxx, 4 / 15 * np.pi * r ** 5, places=1)
        self.assertAlmostEqual(Myy, 4 / 15 * np.pi * r ** 5, places=1)
        self.assertAlmostEqual(Mzz, 4 / 15 * np.pi * r ** 5, places=1)
        self.assertLess(abs(Mxy), 1e-9)
        self.assertLess(abs(Myz), 1e-9)
        self.assertLess(abs(Mzx), 1e-9)
