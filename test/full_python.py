import GMatElastoPlasticFiniteStrainSimo.Cartesian3d as gmat
import numpy as np
import scipy.sparse.linalg as sp
import itertools

# turn of warning for zero division
# (which occurs in the linearization of the logarithmic strain)
np.seterr(divide='ignore', invalid='ignore')

# ---------------------------------------------- GRID ----------------------------------------------

Nx     = 2           # number of voxels in x-direction
Ny     = 1           # number of voxels in y-direction
Nz     = 1           # number of voxels in z-direction
shape  = [Nx,Ny,Nz]  # number of voxels as list: [Nx,Ny,Nz]

# --------------------------------------- TENSOR OPERATIONS ----------------------------------------

# tensor operations / products: np.einsum enables index notation, avoiding loops
# e.g. ddot42 performs $C_ij = A_ijkl B_lk$ for the entire grid
trans2 = lambda A2   : np.einsum('ijxyz          ->jixyz  ',A2   )
ddot22 = lambda A2,B2: np.einsum('ijxyz  ,jixyz  ->xyz    ',A2,B2)
ddot42 = lambda A4,B2: np.einsum('ijklxyz,lkxyz  ->ijxyz  ',A4,B2)
ddot44 = lambda A4,B4: np.einsum('ijklxyz,lkmnxyz->ijmnxyz',A4,B4)
dot11  = lambda A1,B1: np.einsum('ixyz   ,ixyz   ->xyz    ',A1,B1)
dot22  = lambda A2,B2: np.einsum('ijxyz  ,jkxyz  ->ikxyz  ',A2,B2)
dot24  = lambda A2,B4: np.einsum('ijxyz  ,jkmnxyz->ikmnxyz',A2,B4)
dot42  = lambda A4,B2: np.einsum('ijklxyz,lmxyz  ->ijkmxyz',A4,B2)
dyad22 = lambda A2,B2: np.einsum('ijxyz  ,klxyz  ->ijklxyz',A2,B2)
dyad11 = lambda A1,B1: np.einsum('ixyz   ,jxyz   ->ijxyz  ',A1,B1)

# eigenvalue decomposition of 2nd-order tensor: return in convention i,j,x,y,z
# NB requires to swap default order of NumPy (in in/output)
def eig2(A2):
    swap1i    = lambda A1: np.einsum('xyzi ->ixyz ',A1)
    swap2     = lambda A2: np.einsum('ijxyz->xyzij',A2)
    swap2i    = lambda A2: np.einsum('xyzij->ijxyz',A2)
    vals,vecs = np.linalg.eig(swap2(A2))
    vals      = swap1i(vals)
    vecs      = swap2i(vecs)
    return vals,vecs

# logarithm of grid of 2nd-order tensors
def ln2(A2):
    vals,vecs = eig2(A2)
    return sum([np.log(vals[i])*dyad11(vecs[:,i],vecs[:,i]) for i in range(3)])

# exponent of grid of 2nd-order tensors
def exp2(A2):
    vals,vecs = eig2(A2)
    return sum([np.exp(vals[i])*dyad11(vecs[:,i],vecs[:,i]) for i in range(3)])

# determinant of grid of 2nd-order tensors
def det2(A2):
    return (A2[0,0]*A2[1,1]*A2[2,2]+A2[0,1]*A2[1,2]*A2[2,0]+A2[0,2]*A2[1,0]*A2[2,1])-\
           (A2[0,2]*A2[1,1]*A2[2,0]+A2[0,1]*A2[1,0]*A2[2,2]+A2[0,0]*A2[1,2]*A2[2,1])

# inverse of grid of 2nd-order tensors
def inv2(A2):
    A2det = det2(A2)
    A2inv = np.empty([3,3,Nx,Ny,Nz])
    A2inv[0,0] = (A2[1,1]*A2[2,2]-A2[1,2]*A2[2,1])/A2det
    A2inv[0,1] = (A2[0,2]*A2[2,1]-A2[0,1]*A2[2,2])/A2det
    A2inv[0,2] = (A2[0,1]*A2[1,2]-A2[0,2]*A2[1,1])/A2det
    A2inv[1,0] = (A2[1,2]*A2[2,0]-A2[1,0]*A2[2,2])/A2det
    A2inv[1,1] = (A2[0,0]*A2[2,2]-A2[0,2]*A2[2,0])/A2det
    A2inv[1,2] = (A2[0,2]*A2[1,0]-A2[0,0]*A2[1,2])/A2det
    A2inv[2,0] = (A2[1,0]*A2[2,1]-A2[1,1]*A2[2,0])/A2det
    A2inv[2,1] = (A2[0,1]*A2[2,0]-A2[0,0]*A2[2,1])/A2det
    A2inv[2,2] = (A2[0,0]*A2[1,1]-A2[0,1]*A2[1,0])/A2det
    return A2inv

# ------------------------ INITIATE (IDENTITY) TENSORS ------------------------

# identity tensor (single tensor)
i    = np.eye(3)
# identity tensors (grid)
I    = np.einsum('ij,xyz'           ,                  i   ,np.ones([Nx,Ny,Nz]))
I4   = np.einsum('ijkl,xyz->ijklxyz',np.einsum('il,jk',i,i),np.ones([Nx,Ny,Nz]))
I4rt = np.einsum('ijkl,xyz->ijklxyz',np.einsum('ik,jl',i,i),np.ones([Nx,Ny,Nz]))
I4s  = (I4+I4rt)/2.
II   = dyad22(I,I)

# ------------------------------------- CONSTITUTIVE RESPONSE --------------------------------------

# constitutive response to a certain loading and history
# NB: completely uncoupled from the FFT-solver, but implemented as a regular
#     grid of quadrature points, to have an efficient code;
#     each point is completely independent, just evaluated at the same time
def constitutive(F,F_t,be_t,ep_t):

    # function to compute linearization of the logarithmic Finger tensor
    def dln2_d2(A2):
        vals,vecs = eig2(A2)
        K4        = np.zeros([3,3,3,3,Nx,Ny,Nz])
        for m, n in itertools.product(range(3),repeat=2):
            gc  = (np.log(vals[n])-np.log(vals[m]))/(vals[n]-vals[m])
            gc[vals[n]==vals[m]] = (1.0/vals[m])[vals[n]==vals[m]]
            K4 += gc*dyad22(dyad11(vecs[:,m],vecs[:,n]),dyad11(vecs[:,m],vecs[:,n]))
        return K4

    # elastic stiffness tensor
    C4e      = K*II+2.*mu*(I4s-1./3.*II)

    # trial state
    Fdelta   = dot22(F,inv2(F_t))
    be_s     = dot22(Fdelta,dot22(be_t,trans2(Fdelta)))
    lnbe_s   = ln2(be_s)
    tau_s    = ddot42(C4e,lnbe_s)/2.
    taum_s   = ddot22(tau_s,I)/3.
    taud_s   = tau_s-taum_s*I
    taueq_s  = np.sqrt(3./2.*ddot22(taud_s,taud_s))
    N_s      = 3./2.*taud_s/taueq_s
    phi_s    = taueq_s-(tauy0+H*ep_t)
    phi_s    = 1./2.*(phi_s+np.abs(phi_s))

    # return map
    dgamma   = phi_s/(H+3.*mu)
    ep       = ep_t  +   dgamma
    tau      = tau_s -2.*dgamma*N_s*mu
    lnbe     = lnbe_s-2.*dgamma*N_s
    be       = exp2(lnbe)
    J        = det2(F)
    sig      = tau / J

    # consistent tangent operator
    a0       = dgamma*mu/taueq_s
    a1       = mu/(H+3.*mu)
    C4ep     = ((K-2./3.*mu)/2.+a0*mu)*II+(1.-3.*a0)*mu*I4s+2.*mu*(a0-a1)*dyad22(N_s,N_s)
    dlnbe4_s = dln2_d2(be_s)
    dbe4_s   = 2.*dot42(I4s,be_s)
    K4       = (C4e/2.)*(phi_s<=0.).astype(np.float)+C4ep*(phi_s>0.).astype(np.float)
    K4       = ddot44(K4,ddot44(dlnbe4_s,dbe4_s))
    K4       = dot42(-I4rt,tau)+K4

    return sig, K4, be, ep

# phase indicator:
phase  = np.zeros([Nx,Ny,Nz]); phase[1,0,0] = 1.
# function to convert material parameters to grid of scalars
param  = lambda M0,M1: M0*np.ones([Nx,Ny,Nz])*(1.-phase)+\
                       M1*np.ones([Nx,Ny,Nz])*    phase
# material parameters
K      = param(0.833,0.833)  # bulk      modulus
mu     = param(0.386,0.386)  # shear     modulus
H      = param(0.004,0.008)  # hardening modulus
tauy0  = param(1000.,0.003)  # initial yield stress

# --------------------------------------- C++ IMPLEMENTATION ---------------------------------------

Isoft = (1. - phase).astype(np.uint).reshape(2,1)
Ihard = (     phase).astype(np.uint).reshape(2,1)

mat = gmat.Matrix(2,1)
mat.setElastic        (Isoft, K[0,0,0], mu[0,0,0])
mat.setLinearHardening(Ihard, K[1,0,0], mu[1,0,0], tauy0[1,0,0], H[1,0,0])
mat.check()

# -------------------------------------------- LOADING ---------------------------------------------

# stress, deformation gradient, plastic strain, elastic Finger tensor
# NB "_t" signifies that it concerns the value at the previous increment
ep_t = np.zeros([    Nx,Ny,Nz])
P    = np.zeros([3,3,Nx,Ny,Nz])
F    = np.array(I,copy=True)
F_t  = np.array(I,copy=True)
be_t = np.array(I,copy=True)

# initialize macroscopic incremental loading
ninc = 50
lam  = 0.0

# incremental deformation
for inc in range(1,ninc):

    # uniformly set deformation gradient (pure-shear)
    lam   += 0.2/float(ninc)
    F      = np.array(I,copy=True)
    F[0,0] =    (1.+lam)
    F[1,1] = 1./(1.+lam)

    # compute constitutive response
    # - Python
    Sig, C, be, ep = constitutive(F, F_t, be_t, ep_t)
    # - reorder
    C_F = np.zeros((2,1,3,3))
    C_F[0,0,:,:] = F[:,:,0,0,0]
    C_F[1,0,:,:] = F[:,:,1,0,0]
    # - C++
    C_Sig, C_C = mat.Tangent(C_F)
    C_ep = mat.Epsp()
    # - compare
    assert np.abs(ep[0,0,0] - C_ep[0,0]) < 1.e-12
    assert np.abs(ep[1,0,0] - C_ep[1,0]) < 1.e-12
    assert np.linalg.norm(Sig[:,:,0,0,0] - C_Sig[0,0,:,:]) < 1.e-12
    assert np.linalg.norm(Sig[:,:,1,0,0] - C_Sig[1,0,:,:]) < 1.e-12
    assert np.linalg.norm(C[:,:,:,:,0,0,0] - C_C[0,0,:,:,:,:]) < 1.e-12
    assert np.linalg.norm(C[:,:,:,:,1,0,0] - C_C[1,0,:,:,:,:]) < 1.e-12

    # end-of-increment: update history
    # - Python
    F_t  = np.array(F ,copy=True)
    be_t = np.array(be,copy=True)
    ep_t = np.array(ep,copy=True)
    # - C++
    mat.increment()
