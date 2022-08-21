import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import h5py
import numpy as np

with h5py.File("Cartesian3d_random.hdf5", "w") as data:

    shape = [1000, 4]

    plastic = GMat.LinearHardening2d(
        K=np.random.random(shape),
        G=np.random.random(shape),
        tauy0=np.random.random(shape),
        H=np.random.random(shape),
    )

    elastic = GMat.Elastic2d(K=plastic.K, G=plastic.G)

    data["K"] = plastic.K
    data["G"] = plastic.G
    data["tauy0"] = plastic.tauy0
    data["H"] = plastic.H

    I2 = tensor.Array2d(shape).I2

    for i in range(20):

        plastic.increment()
        plastic.F = I2 + np.random.random(shape + [3, 3])
        elastic.F = np.copy(plastic.F)

        data[f"/data/{i:d}/F"] = plastic.F
        data[f"/data/{i:d}/plastic/Stress"] = plastic.Sig
        data[f"/data/{i:d}/plastic/Tangent"] = plastic.C
        data[f"/data/{i:d}/plastic/epsp"] = plastic.epsp
        data[f"/data/{i:d}/elastic/Stress"] = elastic.Sig
        data[f"/data/{i:d}/elastic/Tangent"] = elastic.C
