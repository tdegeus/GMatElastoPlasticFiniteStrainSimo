import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.style.use(["goose", "goose-latex"])
except FileNotFoundError:
    pass


def consistency(mat, plastic):

    # tensor products

    def A2_dot_B2(A2, B2):
        return np.einsum("ij, jk -> ik", A2, B2)

    def A4_dot_B2(A4, B2):
        return np.einsum("ijkl, lm -> ijkm", A4, B2)

    def A2_dot_B4(A2, B4):
        return np.einsum("ij, jkmn -> ikmn", A2, B4)

    def A4_ddot_B2(A4, B2):
        return np.einsum("ijkl, lk -> ij", A4, B2)

    def norm(A2):
        return np.abs(np.einsum("ij, ji", A2, A2))

    # pre-loading

    ninc = 301

    for ilam, lam in enumerate(np.linspace(0.0, 1.0, ninc)):

        if plastic:
            mat.increment()

        F0 = np.array(
            [
                [1.0 + lam, 0.0, 0.0],
                [0.0, 1.0 / (1.0 + lam), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

        mat.F = F0
        Sig0 = np.copy(mat.Sig)
        C = np.copy(mat.C)

        Tau0 = Sig0 * np.linalg.det(F0)

        P0 = A2_dot_B2(Tau0, np.linalg.inv(F0).T)
        K = A4_dot_B2(A2_dot_B4(np.linalg.inv(F0), C), np.linalg.inv(F0).T)

    # consistency check

    x = np.logspace(-16, 0, 100)
    y = np.zeros(x.shape)

    for i in range(len(x)):

        dF = np.random.random((3, 3)) * x[i]
        F = F0 + dF

        mat.F = F
        Sig = np.copy(mat.Sig)

        Tau = Sig * np.linalg.det(F)

        P = A2_dot_B2(Tau, np.linalg.inv(F).T)

        dP = P - P0

        y[i] = norm(dP - (A4_ddot_B2(K, dF.T)).T) / norm(dP)

    return (x, y)


# Plot

fig, ax = plt.subplots()

ax.plot(*consistency(GMat.Elastic0d(10.0, 1.0), False), c="b", label=r"Elastic")

ax.plot(
    *consistency(GMat.LinearHardening0d(10.0, 1.0, 1.0, 1.0), True), c="r", label=r"LinearHardening"
)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlim([1e-16, 1e0])
ax.set_ylim([1e-16, 1e0])

ax.set_xlabel(r"$|| \delta \bm{\varepsilon} ||$")
ax.set_ylabel(r"$\eta$")

gplt.plot_powerlaw(
    -2,
    0.0,
    1.0,
    0.5,
    axis=ax,
    units="relative",
    color="k",
    linewidth=1,
    label=r"rounding error $\sim || \delta \bm{\varepsilon} ||^{-2}$",
)

gplt.plot_powerlaw(
    +2,
    0.5,
    0.0,
    0.5,
    axis=ax,
    units="relative",
    color="k",
    linewidth=1,
    linestyle="--",
    label=r"linearisation error $\sim || \delta \bm{\varepsilon} ||^{+2}$",
)

ax.legend()

fig.savefig("consistency.pdf")
plt.close(fig)
