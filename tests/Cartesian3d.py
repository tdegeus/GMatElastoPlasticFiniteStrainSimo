import unittest
import numpy as np
import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat

class Test_main(unittest.TestCase):

    def test_Elastic(self):

        K = 12.3
        G = 45.6

        gamma = 0.02

        F = np.array([
                [1.0 + gamma, 0.0, 0.0],
                [0.0, 1.0 / (1.0 + gamma), 0.0],
                [0.0, 0.0, 1.0]])

        Sig = np.array([
            [G * 2.0 * np.log(1.0 + gamma), 0.0, 0.0],
            [0.0, - G * 2.0 * np.log(1.0 + gamma), 0.0],
            [0.0, 0.0, 0.0]])

        mat = GMat.Elastic(K, G)
        mat.setDefGrad(F)
        Sig = mat.Stress()

        self.assertTrue(np.allclose(mat.Stress(), Sig))

    def test_Array2d(self):

        K = 12.3
        G = 45.6

        gamma = 0.02

        F = np.array([
                [1.0 + gamma, 0.0, 0.0],
                [0.0, 1.0 / (1.0 + gamma), 0.0],
                [0.0, 0.0, 1.0]])

        Sig = np.array([
            [G * 2.0 * np.log(1.0 + gamma), 0.0, 0.0],
            [0.0, - G * 2.0 * np.log(1.0 + gamma), 0.0],
            [0.0, 0.0, 0.0]])

        nelem = 3
        nip = 2
        mat = GMat.Array2d([nelem, nip])
        ndim = 3

        I = np.ones([nelem, nip], dtype='int')
        mat.setElastic(I, K, G)

        f = np.zeros((nelem, nip, ndim, ndim))
        sig = np.zeros((nelem, nip, ndim, ndim))

        for e in range(nelem):
            for q in range(nip):
                f[e, q, :, :] = F
                sig[e, q, :, :] = Sig

        mat.setDefGrad(f)

        self.assertTrue(np.allclose(mat.Stress(), sig))

if __name__ == '__main__':

    unittest.main()
