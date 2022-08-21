import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.style.use(["goose", "goose-latex"])
except FileNotFoundError:
    pass


def StressStrain(mat, plastic):

    ninc = 301

    epseq = np.zeros(ninc)
    sigeq = np.zeros(ninc)

    for ilam, lam in enumerate(np.linspace(0.0, 1.0, ninc)):

        if plastic:
            mat.increment()

        mat.F = np.array(
            [
                [1.0 + lam, 0.0, 0.0],
                [0.0, 1.0 / (1.0 + lam), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

        epseq[ilam] = GMat.Epseq(GMat.Strain(mat.F))
        sigeq[ilam] = GMat.Sigeq(mat.Sig)

    return (epseq, sigeq)


# Plot

fig, ax = plt.subplots()

elas = GMat.Elastic0d(10.0, 1.0)
plas = GMat.LinearHardening0d(10.0, 1.0, 1.0, 1.0)

ax.plot(*StressStrain(elas, False), c="b", label=r"Elastic")
ax.plot(*StressStrain(plas, True), c="r", label=r"LinearHardening")

ax.plot(ax.get_xlim(), plas.tauy0 * np.ones(2), c="r", ls="--", lw=1.0)
ax.plot(plas.tauy0 / (3.0 * plas.G) * np.ones(2), ax.get_ylim(), c="r", ls="--", lw=1.0)

ax.set_xlabel(r"$\varepsilon_\mathrm{eq}$")
ax.set_ylabel(r"$\sigma_\mathrm{eq}$")

ax.legend()

fig.savefig("stress-strain.pdf")
plt.show()
