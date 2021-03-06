import h5py
import numpy as np
import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import unittest

class Test(unittest.TestCase):

    def test_main(self):

        with h5py.File('Cartesian3d_random.hdf5', 'r') as data:

            mat = GMat.Array2d(data['/shape'][...])

            I = data['/LinearHardening/I'][...]
            idx = data['/LinearHardening/idx'][...]
            K = data['/LinearHardening/K'][...]
            G = data['/LinearHardening/G'][...]
            tauy0 = data['/LinearHardening/tauy0'][...]
            H = data['/LinearHardening/H'][...]

            mat.setLinearHardening(I, idx, K, G, tauy0, H)

            I = data['/Elastic/I'][...]
            idx = data['/Elastic/idx'][...]
            K = data['/Elastic/K'][...]
            G = data['/Elastic/G'][...]

            mat.setElastic(I, idx, K, G)

            for i in range(20):

                mat.increment()

                F = data['/random/{0:d}/F'.format(i)][...]

                mat.setDefGrad(F)

                self.assertTrue(np.allclose(mat.Stress(), data['/random/{0:d}/Stress'.format(i)][...]))
                self.assertTrue(np.allclose(mat.Tangent(), data['/random/{0:d}/Tangent'.format(i)][...]))
                self.assertTrue(np.allclose(mat.Epsp(), data['/random/{0:d}/Epsp'.format(i)][...]))

if __name__ == '__main__':

    unittest.main()
