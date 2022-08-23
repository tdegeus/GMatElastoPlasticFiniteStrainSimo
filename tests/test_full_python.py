import itertools
import unittest

import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import numpy as np

# turn of warning for zero division
# (which occurs in the linearization of the logarithmic strain)
np.seterr(divide="ignore", invalid="ignore")

# ---------------------------------------------- GRID ----------------------------------------------

Nx = 4  # number of voxels in x-direction
Ny = 3  # number of voxels in y-direction
Nz = 2  # number of voxels in z-direction
shape = [Nx, Ny, Nz]  # number of voxels as list: [Nx, Ny, Nz]


def tensor2cpp(A):
    ret = np.empty([Nx, Ny, Nz, 3, 3])
    for i in range(3):
        for j in range(3):
            ret[..., i, j] = A[i, j, ...]
    return ret


def tensor4cpp(A):
    ret = np.empty([Nx, Ny, Nz, 3, 3, 3, 3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for m in range(3):
                    ret[..., i, j, k, m] = A[i, j, k, m, ...]
    return ret


# --------------------------------------- TENSOR OPERATIONS ----------------------------------------

# tensor operations / products: np.einsum enables index notation, avoiding loops
# e.g. ddot42 performs $C_ij = A_ijkl B_lk$ for the entire grid
def trans2(A2):
    return np.einsum("ijxyz->jixyz", A2)


def ddot22(A2, B2):
    return np.einsum("ijxyz,jixyz->xyz", A2, B2)


def ddot42(A4, B2):
    return np.einsum("ijklxyz,lkxyz->ijxyz", A4, B2)


def ddot44(A4, B4):
    return np.einsum("ijklxyz,lkmnxyz->ijmnxyz", A4, B4)


def dot11(A1, B1):
    return np.einsum("ixyz,ixyz->xyz", A1, B1)


def dot22(A2, B2):
    return np.einsum("ijxyz,jkxyz->ikxyz", A2, B2)


def dot24(A2, B4):
    return np.einsum("ijxyz,jkmnxyz->ikmnxyz", A2, B4)


def dot42(A4, B2):
    return np.einsum("ijklxyz,lmxyz  ->ijkmxyz", A4, B2)


def dyad22(A2, B2):
    return np.einsum("ijxyz,klxyz->ijklxyz", A2, B2)


def dyad11(A1, B1):
    return np.einsum("ixyz,jxyz->ijxyz", A1, B1)


# eigenvalue decomposition of 2nd-order tensor: return in convention i,j,x,y,z
# NB requires to swap default order of NumPy (in in/output)
def eig2(A2):
    def swap1i(A1):
        return np.einsum("xyzi->ixyz", A1)

    def swap2(A2):
        return np.einsum("ijxyz->xyzij", A2)

    def swap2i(A2):
        return np.einsum("xyzij->ijxyz", A2)

    vals, vecs = np.linalg.eig(swap2(A2))
    vals = swap1i(vals)
    vecs = swap2i(vecs)
    return vals, vecs


# logarithm of grid of 2nd-order tensors
def ln2(A2):
    vals, vecs = eig2(A2)
    return sum(np.log(vals[i]) * dyad11(vecs[:, i], vecs[:, i]) for i in range(3))


# exponent of grid of 2nd-order tensors
def exp2(A2):
    vals, vecs = eig2(A2)
    return sum(np.exp(vals[i]) * dyad11(vecs[:, i], vecs[:, i]) for i in range(3))


# determinant of grid of 2nd-order tensors
def det2(A2):
    return (
        A2[0, 0] * A2[1, 1] * A2[2, 2]
        + A2[0, 1] * A2[1, 2] * A2[2, 0]
        + A2[0, 2] * A2[1, 0] * A2[2, 1]
    ) - (
        A2[0, 2] * A2[1, 1] * A2[2, 0]
        + A2[0, 1] * A2[1, 0] * A2[2, 2]
        + A2[0, 0] * A2[1, 2] * A2[2, 1]
    )


# inverse of grid of 2nd-order tensors
def inv2(A2):
    A2det = det2(A2)
    A2inv = np.empty([3, 3, Nx, Ny, Nz])
    A2inv[0, 0] = (A2[1, 1] * A2[2, 2] - A2[1, 2] * A2[2, 1]) / A2det
    A2inv[0, 1] = (A2[0, 2] * A2[2, 1] - A2[0, 1] * A2[2, 2]) / A2det
    A2inv[0, 2] = (A2[0, 1] * A2[1, 2] - A2[0, 2] * A2[1, 1]) / A2det
    A2inv[1, 0] = (A2[1, 2] * A2[2, 0] - A2[1, 0] * A2[2, 2]) / A2det
    A2inv[1, 1] = (A2[0, 0] * A2[2, 2] - A2[0, 2] * A2[2, 0]) / A2det
    A2inv[1, 2] = (A2[0, 2] * A2[1, 0] - A2[0, 0] * A2[1, 2]) / A2det
    A2inv[2, 0] = (A2[1, 0] * A2[2, 1] - A2[1, 1] * A2[2, 0]) / A2det
    A2inv[2, 1] = (A2[0, 1] * A2[2, 0] - A2[0, 0] * A2[2, 1]) / A2det
    A2inv[2, 2] = (A2[0, 0] * A2[1, 1] - A2[0, 1] * A2[1, 0]) / A2det
    return A2inv


# ------------------------ INITIATE (IDENTITY) TENSORS ------------------------

# identity tensor (single tensor)
i = np.eye(3)
# identity tensors (grid)
I2 = np.einsum("ij,xyz", i, np.ones([Nx, Ny, Nz]))
I4 = np.einsum("ijkl,xyz->ijklxyz", np.einsum("il,jk", i, i), np.ones([Nx, Ny, Nz]))
I4rt = np.einsum("ijkl,xyz->ijklxyz", np.einsum("ik,jl", i, i), np.ones([Nx, Ny, Nz]))
I4s = (I4 + I4rt) / 2.0
II = dyad22(I2, I2)

# ------------------------------------- CONSTITUTIVE RESPONSE --------------------------------------


# constitutive response to a certain loading and history
# NB: completely uncoupled from the FFT-solver, but implemented as a regular
#     grid of quadrature points, to have an efficient code;
#     each point is completely independent, just evaluated at the same time
def constitutive(F, F_t, be_t, ep_t, K, mu, tauy0, H):

    # function to compute linearization of the logarithmic Finger tensor
    def dln2_d2(A2):
        vals, vecs = eig2(A2)
        K4 = np.zeros([3, 3, 3, 3, Nx, Ny, Nz])
        for m, n in itertools.product(range(3), repeat=2):
            gc = (np.log(vals[n]) - np.log(vals[m])) / (vals[n] - vals[m])
            gc[vals[n] == vals[m]] = (1.0 / vals[m])[vals[n] == vals[m]]
            K4 += gc * dyad22(dyad11(vecs[:, m], vecs[:, n]), dyad11(vecs[:, m], vecs[:, n]))
        return K4

    # elastic stiffness tensor
    C4e = K * II + 2.0 * mu * (I4s - 1.0 / 3.0 * II)

    # trial state
    Fdelta = dot22(F, inv2(F_t))
    be_s = dot22(Fdelta, dot22(be_t, trans2(Fdelta)))
    lnbe_s = ln2(be_s)
    tau_s = ddot42(C4e, lnbe_s) / 2.0
    taum_s = ddot22(tau_s, I2) / 3.0
    taud_s = tau_s - taum_s * I2
    taueq_s = np.sqrt(3.0 / 2.0 * ddot22(taud_s, taud_s))
    N_s = 3.0 / 2.0 * taud_s / taueq_s
    phi_s = taueq_s - (tauy0 + H * ep_t)
    phi_s = 1.0 / 2.0 * (phi_s + np.abs(phi_s))

    # return map
    dgamma = phi_s / (H + 3.0 * mu)
    ep = ep_t + dgamma
    tau = tau_s - 2.0 * dgamma * N_s * mu
    lnbe = lnbe_s - 2.0 * dgamma * N_s
    be = exp2(lnbe)
    J = det2(F)
    sig = tau / J

    # consistent tangent operator
    a0 = dgamma * mu / taueq_s
    a1 = mu / (H + 3.0 * mu)
    C4ep = (
        ((K - 2.0 / 3.0 * mu) / 2.0 + a0 * mu) * II
        + (1.0 - 3.0 * a0) * mu * I4s
        + 2.0 * mu * (a0 - a1) * dyad22(N_s, N_s)
    )
    dlnbe4_s = dln2_d2(be_s)
    dbe4_s = 2.0 * dot42(I4s, be_s)
    K4 = (C4e / 2.0) * (phi_s <= 0.0).astype(np.float64) + C4ep * (phi_s > 0.0).astype(np.float64)
    K4 = ddot44(K4, ddot44(dlnbe4_s, dbe4_s))
    K4 = dot42(-I4rt, tau) + K4

    return sig, K4 / J, be, ep


class Test_main(unittest.TestCase):
    """ """

    def test_Elastic(self):

        # material parameters
        K = np.random.random(shape)
        mu = np.random.random(shape)
        tauy0 = 100000 * np.random.random(shape)
        H = np.random.random(shape)
        mat = GMat.Elastic3d(K, mu)

        # stress, deformation gradient, plastic strain, elastic Finger tensor
        # NB "_t" signifies that it concerns the value at the previous increment
        ep_t = np.zeros([Nx, Ny, Nz])
        F = np.array(I2, copy=True)
        F_t = np.array(I2, copy=True)
        be_t = np.array(I2, copy=True)

        # initialize macroscopic incremental loading
        ninc = 100

        # incremental deformation
        for inc in range(1, ninc):

            # deformation gradient
            F = np.copy(I2) + 0.1 * inc * np.random.random([3, 3] + [Nx, Ny, Nz])

            # compute constitutive response
            # - Python
            Sig, C, be, ep = constitutive(F, F_t, be_t, ep_t, K, mu, tauy0, H)
            # - C++
            mat.F = tensor2cpp(F)
            # - compare
            self.assertTrue(np.allclose(mat.Sig, tensor2cpp(Sig), atol=1e-4, rtol=1e-2))
            self.assertTrue(np.allclose(mat.C, tensor4cpp(C), atol=1e-4, rtol=1e-2))

            # end-of-increment: update history
            # - Python
            F_t = np.array(F, copy=True)
            be_t = np.array(be, copy=True)
            ep_t = np.array(ep, copy=True)

    def test_LinearHardening(self):

        # material parameters
        K = np.random.random(shape)
        mu = np.random.random(shape)
        tauy0 = np.random.random(shape)
        H = np.random.random(shape)
        mat = GMat.LinearHardening3d(K, mu, tauy0, H)

        # stress, deformation gradient, plastic strain, elastic Finger tensor
        # NB "_t" signifies that it concerns the value at the previous increment
        ep_t = np.zeros([Nx, Ny, Nz])
        F = np.array(I2, copy=True)
        F_t = np.array(I2, copy=True)
        be_t = np.array(I2, copy=True)

        # initialize macroscopic incremental loading
        ninc = 100

        # incremental deformation
        for inc in range(1, ninc):

            # deformation gradient
            F = np.copy(I2) + 0.1 * inc * np.random.random([3, 3] + [Nx, Ny, Nz])

            # compute constitutive response
            # - Python
            Sig, C, be, ep = constitutive(F, F_t, be_t, ep_t, K, mu, tauy0, H)
            # - C++
            mat.F = tensor2cpp(F)
            # - compare
            self.assertTrue(np.allclose(mat.epsp, ep, atol=1e-4, rtol=1e-2))
            self.assertTrue(np.allclose(mat.Sig, tensor2cpp(Sig), atol=1e-4, rtol=1e-2))
            self.assertTrue(np.allclose(mat.C, tensor4cpp(C), atol=1e-4, rtol=1e-2))

            # end-of-increment: update history
            # - Python
            F_t = np.array(F, copy=True)
            be_t = np.array(be, copy=True)
            ep_t = np.array(ep, copy=True)
            # - C++
            mat.increment()


if __name__ == "__main__":

    unittest.main()
