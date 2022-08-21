import unittest

import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import numpy as np


class Test_main(unittest.TestCase):
    """ """

    def test_Epseq_Sigeq(self):

        A = np.zeros((2, 3, 3, 3))
        A[..., 0, 1] = 1
        A[..., 1, 0] = 1

        self.assertTrue(np.allclose(GMat.Epseq(A), 2 / np.sqrt(3) * np.ones(A.shape[:-2])))
        self.assertTrue(np.allclose(GMat.Sigeq(A), np.sqrt(3.0) * np.ones(A.shape[:-2])))

    def test_Strain(self):

        shape = [2, 3]
        gamma = np.random.random(shape)
        F = tensor.Array2d(shape).I2
        F[..., 0, 0] = 1 + gamma
        F[..., 1, 1] = 1 / (1 + gamma)

        Eps = np.zeros_like(F)
        Eps[..., 0, 0] = np.log(1 + gamma)
        Eps[..., 1, 1] = -np.log(1 + gamma)

        self.assertTrue(np.allclose(GMat.Strain(F), Eps))

    def test_Elastic(self):

        shape = [2, 3]
        mat = GMat.Elastic2d(
            K=np.random.random(shape),
            G=np.random.random(shape),
        )

        gamma = np.random.random(shape)

        mat.F[..., 0, 0] = 1 + gamma
        mat.F[..., 1, 1] = 1 / (1 + gamma)
        mat.refresh()

        Sig = np.zeros_like(mat.F)
        Sig[..., 0, 0] = 2 * mat.G * np.log(1 + gamma)
        Sig[..., 1, 1] = -2 * mat.G * np.log(1 + gamma)

        self.assertTrue(np.allclose(mat.Sig, Sig))


if __name__ == "__main__":

    unittest.main()
