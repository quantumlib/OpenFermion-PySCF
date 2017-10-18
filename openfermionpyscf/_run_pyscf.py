#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""Class and functions to store quantum chemistry data."""
from __future__ import absolute_import

from functools import reduce

import numpy
import pyscf
from pyscf import ci, cc, fci, mp

from openfermion.config import *


def prepare_pyscf_molecule(molecule):
    """
    This function creates and saves a pyscf input file.

    Args:
        molecule: An instance of the MolecularData class.

    Returns:
        pyscf_molecule: A pyscf molecule instance.
    """
    pyscf_molecule = pyscf.gto.Mole()
    pyscf_molecule.atom = molecule.geometry
    pyscf_molecule.basis = molecule.basis
    pyscf_molecule.spin = molecule.multiplicity - 1
    pyscf_molecule.charge = molecule.charge
    pyscf_molecule.symmetry = False
    pyscf_molecule.build()

    return pyscf_molecule


def compute_scf(pyscf_molecule):
    """
    Perform a Hartree-Fock calculation.

    Args:
        pyscf_molecule: A pyscf molecule instance.

    Returns:
        pyscf_scf: A PySCF "SCF" calculation object.
    """
    if pyscf_molecule.spin:
        pyscf_scf = pyscf.scf.ROHF(pyscf_molecule)
    else:
        pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    return pyscf_scf


def compute_integrals(pyscf_molecule, pyscf_scf):
    """
    Compute the 1-electron and 2-electron integrals.

    Args:
        pyscf_molecule: A pyscf molecule instance.
        pyscf_scf: A PySCF "SCF" calculation object.

    Returns:
        one_electron_integrals: An N by N array storing h_{pq}
        two_electron_integrals: An N by N by N by N array storing h_{pqrs}.
    """
    # Get one electrons integrals.
    n_orbitals = pyscf_scf.mo_coeff.shape[0]
    one_electron_compressed = reduce(numpy.dot, (pyscf_scf.mo_coeff.T,
                                                 pyscf_scf.get_hcore(),
                                                 pyscf_scf.mo_coeff))
    one_electron_integrals = one_electron_compressed.reshape(
        n_orbitals, n_orbitals).astype(float)

    # Get two electron integrals in compressed format.
    two_electron_compressed = pyscf.ao2mo.kernel(pyscf_molecule,
                                                 pyscf_scf.mo_coeff)
    two_electron_integrals = numpy.empty((n_orbitals, n_orbitals,
                                          n_orbitals, n_orbitals))

    # Unpack symmetry.
    n_pairs = n_orbitals * (n_orbitals + 1) // 2
    if two_electron_compressed.ndim == 2:

        # Case of 4-fold symmetry.
        assert(two_electron_compressed.size == n_pairs ** 2)
        pq = 0
        for p in range(n_orbitals):
            for q in range(p + 1):
                rs = 0
                for r in range(n_orbitals):
                    for s in range(r + 1):
                        pqrs_value = two_electron_compressed[pq, rs]
                        two_electron_integrals[p, s, r, q] = float(pqrs_value)
                        two_electron_integrals[q, s, r, p] = float(pqrs_value)
                        two_electron_integrals[p, r, s, q] = float(pqrs_value)
                        two_electron_integrals[q, r, s, p] = float(pqrs_value)
                        rs += 1
                pq += 1
    else:

        # Case of 8-fold symmetry.
        assert(two_electron_compressed.size == n_pairs * (n_pairs + 1) // 2)
        pq = 0
        pqrs = 0
        for p in range(n_orbitals):
            for q in range(p + 1):
                rs = 0
                for r in range(p + 1):
                    for s in range(r + 1):
                        if pq >= rs:
                            pqrs_value = two_electron_compressed[pqrs]
                            two_electron_integrals[p, s,
                                                   r, q] = float(pqrs_value)
                            two_electron_integrals[q, s,
                                                   r, p] = float(pqrs_value)
                            two_electron_integrals[p, r,
                                                   s, q] = float(pqrs_value)
                            two_electron_integrals[q, r,
                                                   s, p] = float(pqrs_value)
                            two_electron_integrals[s, p,
                                                   q, r] = float(pqrs_value)
                            two_electron_integrals[s, q,
                                                   p, r] = float(pqrs_value)
                            two_electron_integrals[r, p,
                                                   q, s] = float(pqrs_value)
                            two_electron_integrals[r, q,
                                                   p, s] = float(pqrs_value)
                            pqrs += 1
                        rs += 1
                pq += 1

    # Return.
    return one_electron_integrals, two_electron_integrals


def run_pyscf(molecule,
              run_scf=True,
              run_mp2=False,
              run_cisd=False,
              run_ccsd=False,
              run_fci=False,
              verbose=False):
    """
    This function runs a pyscf calculation.

    Args:
        molecule: An instance of the MolecularData class.
        run_scf: Optional boolean to run SCF calculation.
        run_mp2: Optional boolean to run MP2 calculation.
        run_cisd: Optional boolean to run CISD calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        run_fci: Optional boolean to FCI calculation.
        verbose: Boolean whether to print calculation results to screen.

    Returns:
        molecule: The updated MolecularData object.
    """
    # Prepare pyscf molecule.
    pyscf_molecule = prepare_pyscf_molecule(molecule)
    molecule.n_orbitals = int(pyscf_molecule.nbas)
    molecule.n_qubits = 2 * molecule.n_orbitals
    molecule.nuclear_repulsion = float(pyscf.gto.energy_nuc(pyscf_molecule))

    # Run SCF.
    pyscf_scf = compute_scf(pyscf_molecule)
    pyscf_scf.verbose = 0
    molecule.hf_energy = float(pyscf_scf.kernel())
    if verbose:
        print('Hartree-Fock energy for {} ({} electrons) is {}.'.format(
            molecule.name, molecule.n_electrons, molecule.hf_energy))

    # Populate fields.
    molecule.canonical_orbitals = pyscf_scf.mo_coeff.astype(float)
    molecule.orbital_energies = pyscf_scf.mo_energy.astype(float)

    # Get integrals.
    one_body_integrals, two_body_integrals = compute_integrals(
        pyscf_molecule, pyscf_scf)
    molecule.one_body_integrals = one_body_integrals
    molecule.two_body_integrals = two_body_integrals
    molecule.overlap_integrals = pyscf_scf.get_ovlp()

    # Run MP2.
    if run_mp2:
        pyscf_mp2 = pyscf.mp.MP2(pyscf_scf)
        pyscf_mp2.verbose = 0
        molecule.mp2_energy = molecule.hf_energy + pyscf_mp2.kernel()[0]
        if verbose:
            print('MP2 energy for {} ({} electrons) is {}.'.format(
                molecule.name, molecule.n_electrons, molecule.mp2_energy))

    # Disable non-singlets for now for CI and CC due to not being implemented
    if (run_cisd or run_ccsd) and molecule.multiplicity != 1:
        print("WARNING: CISD and CCSD not implemented for "
              "non-singlet states in PySCF.")
        run_cisd = run_ccsd = False

    # Run CISD.
    if run_cisd:
        pyscf_cisd = pyscf.ci.CISD(pyscf_scf)
        pyscf_cisd.verbose = 0
        pyscf_cisd.kernel()
        molecule.cisd_energy = molecule.hf_energy + pyscf_cisd.e_corr
        if verbose:
            print('CISD energy for {} ({} electrons) is {}.'.format(
                molecule.name, molecule.n_electrons, molecule.cisd_energy))

    # Run CCSD.
    if run_ccsd:
        pyscf_ccsd = pyscf.cc.CCSD(pyscf_scf)
        pyscf_ccsd.verbose = 0
        pyscf_ccsd.kernel()
        molecule.ccsd_energy = molecule.hf_energy + pyscf_ccsd.e_corr
        if verbose:
            print('CCSD energy for {} ({} electrons) is {}.'.format(
                molecule.name, molecule.n_electrons, molecule.ccsd_energy))

    # Run FCI.
    if run_fci:
        pyscf_fci = pyscf.fci.FCI(pyscf_molecule, pyscf_scf.mo_coeff)
        pyscf_fci.verbose = 0
        molecule.fci_energy = pyscf_fci.kernel()[0]
        if verbose:
            print('FCI energy for {} ({} electrons) is {}.'.format(
                molecule.name, molecule.n_electrons, molecule.fci_energy))

    # Return updated molecule instance.
    molecule.save()
    return molecule
