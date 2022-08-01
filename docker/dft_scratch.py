#http://forum.psicode.org/t/dft-from-scratch-doesnt-converge-to-the-correct-result/2423/3

import os,sys
import numpy as np
import math
sys.path.insert(1, os.path.abspath('/usr/local/psi4/lib/'))
import psi4

functional  = "SVWN"
basis       = "sto-3g"
E_conv      = 1e-6
D_conv      = 1e-6
SCF_maxiter = 100

n_orbs = 8
n_elec = n_orbs
n_occ  = n_elec//2
R      = 2.5

psi4.core.clean()
# Define the geometry
psi4.core.set_output_file("H{}_R{}_{}_{}_Psi4.dat".format(n_orbs,R,basis,functional),True)
geometry = "0 1\n"
for d in range(n_orbs//2):
    geometry += "H 0. 0. {}\n".format(- (R/2. + d*R))
    geometry += "H 0. 0. {}\n".format(+ (R/2. + d*R))
geometry += "symmetry c1\n"
geometry += "nocom\n"
geometry += "noreorient\n"

psi4.geometry(geometry)

psi4.set_options({'basis': basis,'save_jk':True, 'debug':1, 'print':5})
dft_e, dft_wfn = psi4.energy(functional, return_wfn=True)

Hcore = dft_wfn.H().clone()
E_nuc = dft_wfn.get_energies('Nuclear')
mints = psi4.core.MintsHelper(dft_wfn.basisset())
I = np.asarray(mints.ao_eri())

# Overlap matrix in the AO basis:
S_AO = dft_wfn.S().np
# KS-MO coefficient matrix in the AO basis:
C_AO = dft_wfn.Ca().np
# Compute the inverse square root of the overlap matrix S
S_eigval, S_eigvec = np.linalg.eigh(S_AO)
S_sqrt_inv = S_eigvec @ np.diag((S_eigval)**(-1./2.)) @ S_eigvec.T
C_transformation = np.linalg.inv(S_sqrt_inv)

# Construct the SAD Guess for the initial density D_AO
psi4.core.prepare_options_for_module("SCF")
sad_basis_list = psi4.core.BasisSet.build(dft_wfn.molecule(), "ORBITAL",
    psi4.core.get_global_option("BASIS"), puream=dft_wfn.basisset().has_puream(),
                                     return_atomlist=True)
sad_fitting_list = psi4.core.BasisSet.build(dft_wfn.molecule(), "DF_BASIS_SAD",
    psi4.core.get_option("SCF", "DF_BASIS_SAD"), puream=dft_wfn.basisset().has_puream(),
                                       return_atomlist=True)

# Use Psi4 SADGuess object to build the SAD Guess
SAD = psi4.core.SADGuess.build_SAD(dft_wfn.basisset(), sad_basis_list)
SAD.set_atomic_fit_bases(sad_fitting_list)
SAD.compute_guess()
D_AO = SAD.Da()

# Initialize the potential object
V_xc = dft_wfn.Da().clone()
Vpotential = dft_wfn.V_potential()

Etot        = 0
Delta_E     = 1e8
convergence = True
dRMS        = 1e8
Fock_list   = []
DIIS_error  = []

for SCF_ITER in range(1, SCF_maxiter+1):

    # Compute the Coulomb potential
    J_coulomb = np.einsum('pqrs,rs->pq', I, D_AO.np)
    # Compute the XC potential with the new density matrix in the AO basis
    Vpotential.set_D([D_AO])
    Vpotential.compute_V([V_xc])

    # Compute the Fock matrix in the AO basis:
    F_AO = Hcore.np + 2*J_coulomb + V_xc.np

    # DIIS
    diis_e = np.einsum('ij,jk,kl->il', F_AO, D_AO.np, S_AO) - np.einsum('ij,jk,kl->il', S_AO, D_AO.np, F_AO)
    diis_e = S_sqrt_inv @ diis_e @ S_sqrt_inv
    Fock_list.append(F_AO)
    DIIS_error.append(diis_e)
    dRMS = np.mean(diis_e**2)**0.5

    if SCF_ITER >= 2:

        # Limit size of DIIS vector
        diis_count = len(Fock_list)
        if diis_count > 6:
            # Remove oldest vector
            del Fock_list[0]
            del DIIS_error[0]
            diis_count -= 1

        # Build error matrix B, [Pulay:1980:393], Eqn. 6, LHS
        B = np.empty((diis_count + 1, diis_count + 1))
        B[-1, :] = -1
        B[:, -1] = -1
        B[-1, -1] = 0
        for num1, e1 in enumerate(DIIS_error):
            for num2, e2 in enumerate(DIIS_error):
                if num2 > num1: continue
                val = np.einsum('ij,ij->', e1, e2)
                B[num1, num2] = val
                B[num2, num1] = val

        # normalize
        B[:-1, :-1] /= np.abs(B[:-1, :-1]).max()

        # Build residual vector, [Pulay:1980:393], Eqn. 6, RHS
        resid = np.zeros(diis_count + 1)
        resid[-1] = -1

        # Solve Pulay equations, [Pulay:1980:393], Eqn. 6
        ci = np.linalg.solve(B, resid)

        # Calculate new fock matrix as linear
        # combination of previous fock matrices
        F_AO = np.zeros_like(F_AO)
        for num, c in enumerate(ci[:-1]):
            F_AO += c * Fock_list[num]

    # Build the Fock matrix in the OAO basis:
    F_OAO = S_sqrt_inv @ F_AO @ S_sqrt_inv.T
    eigvals,eigvecs = np.linalg.eigh(F_OAO)

    # Transform back to the AO basis:
    C_occ = S_sqrt_inv @ eigvecs[:,:n_occ] # AO --> OAO = C_transformation = S^{-1/2}, OAO --> AO = S^{1/2}

    # Compute the density matrix
    D_AO.np[:] = np.einsum('pi,qi->pq', C_occ, C_occ, optimize=True)

    # Compute the Hxc energy with the new density:
    Hxcpot_energy = 2*np.einsum('pq,pq->', (J_coulomb + V_xc.np), D_AO.np) # factors 2 because D_AO is only alpha-D_AO.
    Vpotential.set_D([D_AO])
    Vpotential.compute_V([V_xc]) # otherwise it doesn't change the EHxc energy...
    EHxc = Vpotential.quadrature_values()["FUNCTIONAL"]
    # Compute the Hxc potential contribution (old potential * new density !)
    Etot_new = 2*np.sum(eigvals[:n_occ]) + EHxc - Hxcpot_energy + E_nuc
    Delta_E  = abs(Etot_new - Etot)
    Etot     = Etot_new

    # Print results
    print("*"*10 + " ITERATION {:3d} ".format(SCF_ITER) + "*"*10)
    print("J matrix         : {}".format(J_coulomb))
    #print("F matrix         : {}".format(F_AO))
    #print("V(xc) matrix     : {}".format(V_xc.np))
    #print("C_occ matrix     : {}".format(C_occ))
    #print("D_AO matrix      : {}".format(D_AO.np))
    #print("KS energies      : {}".format(eigvals))
    #print("KS energy orbs   : {:16.8f}".format(2*np.sum(eigvals[:n_occ])))
    #print("Hxcpot energy    : {:16.8f}".format(Hxcpot_energy))
    print("Functional energy: {}".format(EHxc))
    print("Energy (hartree) : {:16.8f}".format(Etot))
    print("Delta Energy     : {:16.8f}".format(Delta_E))
    print("dRMS             : {:16.8f}\n".format(dRMS))

    if (Delta_E < E_conv) and (dRMS < D_conv):
        print("*"*10 + " SUCCESS " + "*"*10)
        print("R = ",R)
        print("Iteration        : {:16d}".format(SCF_ITER))
        print("DFT energy       : {:16.8f}".format(Etot))
        print("DFT energy psi4  : {:16.8f}".format(dft_e))
        break

    if SCF_ITER == SCF_maxiter:
        psi4.core.clean()
        raise Exception("Maximum number of SCF cycles exceeded.")
