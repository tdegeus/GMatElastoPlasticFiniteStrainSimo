
import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import numpy as np


def EQ(a,b):
    assert np.abs(a-b) < 1.e-12

def ALLEQ(a, b):
    assert np.allclose(a, b)


K = 12.3
G = 45.6

gamma = 0.02

F = np.array([
        [1.0 + gamma, 0.0, 0.0],
        [0.0, 1.0 / (1.0 + gamma), 0.0],
        [0.0, 0.0, 1.0]])


# Elastic


mat = GMat.Elastic(K, G)

Sig = mat.Stress(F)

EQ(Sig[0,0], G * +2.0 * np.log(1.0 + gamma))
EQ(Sig[1,1], G * -2.0 * np.log(1.0 + gamma))
EQ(Sig[2,2], 0.0)
EQ(Sig[0,1], 0.0)
EQ(Sig[0,2], 0.0)
EQ(Sig[1,0], 0.0)
EQ(Sig[1,2], 0.0)
EQ(Sig[2,0], 0.0)
EQ(Sig[2,1], 0.0)


# Matrix


nelem = 2
nip = 2

mat = GMat.Matrix(nelem, nip)

# all rows: elastic
I = np.ones([nelem, nip], dtype='int')
mat.setElastic(I, K, G)

f = np.zeros((nelem, nip, 3, 3))
for i in range(3):
    for j in range(3):
        f[:, :, i, j] = F[i, j]

Sig = mat.Stress(f)

for e in range(nelem):
    for q in range(nip):
        EQ(Sig[e,q,0,0], G * +2.0 * np.log(1.0 + gamma))
        EQ(Sig[e,q,1,1], G * -2.0 * np.log(1.0 + gamma))

ALLEQ(Sig[:,:,0,1], 0.0)
ALLEQ(Sig[:,:,1,0], 0.0)
ALLEQ(Sig[:,:,:,2], 0.0)
ALLEQ(Sig[:,:,2,:], 0.0)


print('All checks passed')
