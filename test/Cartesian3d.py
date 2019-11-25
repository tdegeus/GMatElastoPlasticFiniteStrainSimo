
import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as GMat
import numpy as np

# ==================================================================================================

def EQ(a,b):
  assert np.abs(a-b) < 1.e-12

def ALLEQ(a, b):
  assert np.allclose(a, b)

# ==================================================================================================

# material model
# - parameters
K = 12.3
G = 45.6
# - model
mat = GMat.Elastic(K,G)

# simple shear + volumetric deformation
# - parameters
gamma = 0.02
# - strain
F = [[1.0 + gamma, 0.0                , 0.0],
     [0.0        , 1.0 / (1.0 + gamma), 0.0],
     [0.0        , 0.0                , 1.0]]
# - stress
Sig = mat.Stress(F)
# - analytical solution
EQ(Sig[0,0], G * +2.0 * np.log(1.0 + gamma))
EQ(Sig[1,1], G * -2.0 * np.log(1.0 + gamma))
EQ(Sig[2,2], 0)
EQ(Sig[0,1], 0)
EQ(Sig[0,2], 0)
EQ(Sig[1,0], 0)
EQ(Sig[1,2], 0)
EQ(Sig[2,0], 0)
EQ(Sig[2,1], 0)

# ==================================================================================================

# parameters
K = 12.3
G = 45.6

# allocate matrix
nelem = 2
nip = 2
mat = GMat.Matrix(nelem, nip)

# all rows: elastic
I = np.ones([nelem, nip], dtype='int')
mat.setElastic(I,K,G)

# simple shear + volumetric deformation
# - parameters
gamma = 0.02;
# - strain
F = np.zeros((nelem, nip, 3, 3))
F[:,:,0,0] = 1.0 + gamma
F[:,:,1,1] = 1.0 / (1.0 + gamma)
F[:,:,2,2] = 1.0
# - stress
Sig = mat.Stress(F)

# - analytical solution
EQ(Sig[0,0,0,0], G * +2.0 * np.log(1.0 + gamma)); EQ(sig[0,1,0,0], G * +2.0 * np.log(1.0 + gamma))
EQ(Sig[0,0,1,1], G * -2.0 * np.log(1.0 + gamma)); EQ(sig[0,1,1,1], G * -2.0 * np.log(1.0 + gamma))
ALLEQ(Sig[:,:,0,1], 0)
ALLEQ(Sig[:,:,1,1], 0)
ALLEQ(Sig[:,:,:,2], 0)
ALLEQ(Sig[:,:,2,:], 0)

# ==================================================================================================

print('All checks passed')
