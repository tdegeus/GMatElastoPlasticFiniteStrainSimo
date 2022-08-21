import unittest

import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.hdf5") as data:

            K = data["K"][...]
            G = data["G"][...]
            tauy0 = data["tauy0"][...]
            H = data["H"][...]

            plastic = GMat.Matrix(K.shape[0], K.shape[1])
            elastic = GMat.Matrix(K.shape[0], K.shape[1])

            for i in range(K.shape[0]):
                for j in range(K.shape[1]):
                    iden = np.zeros(K.shape, dtype=bool)
                    iden[i, j] = True
                    plastic.setLinearHardening(iden, K[i, j], G[i, j], tauy0[i, j], H[i, j])
                    elastic.setElastic(iden, K[i, j], G[i, j])

            for i in range(20):

                plastic.increment()
                F = data[f"/data/{i:d}/F"][...]

                self.assertTrue(
                    np.allclose(plastic.Stress(F), data[f"/data/{i:d}/plastic/Stress"][...], 1e-3)
                )
                self.assertTrue(
                    np.allclose(
                        plastic.Tangent(F)[1], data[f"/data/{i:d}/plastic/Tangent"][...], 1e-3
                    )
                )
                self.assertTrue(
                    np.allclose(plastic.Epsp(), data[f"/data/{i:d}/plastic/epsp"][...], 1e-3)
                )

                self.assertTrue(
                    np.allclose(elastic.Stress(F), data[f"/data/{i:d}/elastic/Stress"][...], 1e-3)
                )
                self.assertTrue(
                    np.allclose(
                        elastic.Tangent(F)[1], data[f"/data/{i:d}/elastic/Tangent"][...], 1e-3
                    )
                )


if __name__ == "__main__":

    unittest.main()
