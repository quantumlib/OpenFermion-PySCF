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

"""Class to use pyscf program to access quantum chemistry data."""

import numpy
from functools import reduce
from pyscf import ao2mo
from pyscf import scf
from pyscf.cc.addons import spatial2spin
from openfermion.hamiltonians import MolecularData


class PyscfMolecularData(MolecularData):

    """A derived class from openfermion.hamiltonians.MolecularData. This class
    is created to store the PySCF method objects as well as molecule data from
    a fixed basis set at a fixed geometry that is obtained from PySCF
    electronic structure packages. This class provides an interface to access
    the PySCF Hartree-Fock, MP, CI, Coupled-Cluster methods and their energies,
    density matrices and wavefunctions.

    Attributes:
        _pyscf_data(dict): To store PySCF method objects temporarily.
    """
    def __init__(self, geometry=None, basis=None, multiplicity=None,
                 charge=0, description="", filename="", data_directory=None):
        MolecularData.__init__(self, geometry, basis, multiplicity,
                               charge, description, filename, data_directory)
        self._pyscf_data = {}

    @property
    def canonical_orbitals(self):
        """Hartree-Fock canonical orbital coefficients (represented on AO
        basis).
        """
        if self._canonical_orbitals is None:
            scf = self._pyscf_data.get('scf', None)
            self._canonical_orbitals = scf.mo_coeff
        return self._canonical_orbitals

    @property
    def overlap_integrals(self):
        """Overlap integrals for AO basis."""
        if self._overlap_integrals is None:
            scf = self._pyscf_data.get('scf', None)
            self._overlap_integrals = scf.get_ovlp()
        return self._overlap_integrals

    @property
    def one_body_integrals(self):
        """A 2D array for one-body Hamiltonian (H_core) in the MO
        representation."""
        if self._one_body_integrals is None:
            scf = self._pyscf_data.get('scf', None)
            mo = self.canonical_orbitals
            h_core = scf.get_hcore()
            self._one_body_integrals = reduce(numpy.dot, (mo.T, h_core, mo))
        return self._one_body_integrals

    @property
    def two_body_integrals(self):
        """A 4-dimension array for electron repulsion integrals in the MO
        representation.  The integrals are computed as
        h[p,q,r,s]=\int \phi_p(x)* \phi_q(y)* V_{elec-elec} \phi_r(y) \phi_s(x) dxdy
        """
        if self._two_body_integrals is None:
            mol = self._pyscf_data.get('mol', None)
            mo = self.canonical_orbitals
            n_orbitals = mo.shape[1]

            eri = ao2mo.kernel(mol, mo)
            eri = ao2mo.restore(1, # no permutation symmetry
                                eri, n_orbitals)
            # See PQRS convention in OpenFermion.hamiltonians.molecular_data
            # h[p,q,r,s] = (ps|qr) = pyscf_eri[p,s,q,r]
            self._two_body_integrals = numpy.asarray(
                eri.transpose(0, 2, 3, 1), order='C')
        return self._two_body_integrals

    @property
    def cisd_one_rdm(self):
        r"""A 2-dimension array for CISD one-particle density matrix in the MO
        representation.  d[p,q] = < a^\dagger_p a_q >
        """
        if self._cisd_one_rdm is None:
            cisd = self._pyscf_data.get('cisd', None)
            if cisd is None:
                return None

            mf = self._pyscf_data.get('scf', None)
            if isinstance(mf, scf.uhf.UHF):
                raise ValueError('Spin trace for UCISD density matrix.')

            rdm1 = cisd.make_rdm1()
            if isinstance(mf, scf.rohf.ROHF):
                rdm1 = rdm1[0] + rdm1[1]

# pyscf one_rdm is computed as dm1[p,q] = <a^\dagger_q a_p>
            self._cisd_one_rdm = rdm1.T
        return self._cisd_one_rdm

    @property
    def cisd_two_rdm(self):
        r"""A 4-dimension array for CISD two-particle density matrix in the MO
        representation.  D[p,q,r,s] = < a^\dagger_p a^\dagger_q a_r a_s >
        """
        if self._cisd_two_rdm is None:
            cisd = self._pyscf_data.get('cisd', None)
            if cisd is None:
                return None

            mf = self._pyscf_data.get('scf', None)
            if isinstance(mf, scf.uhf.UHF):
                raise ValueError('Spin trace for UCISD density matrix.')

            rdm2 = cisd.make_rdm2()
            if isinstance(mf, scf.rohf.ROHF):
                aa, ab, bb = rdm2
                rdm2 = aa + bb + ab + ab.transpose(2, 3, 0, 1)

# pyscf.ci.cisd.make_rdm2 convention
#       dm2[p,s,q,r] = <a^\dagger_p a^\dagger_q a_r a_s>.
# the two_body_tensor in openfermion.ops._interaction_rdm.InteractionRDM
#       tbt[p,q,r,s] = <a^\dagger_p a^\dagger_q a_r a_s>.
            self._cisd_two_rdm = rdm2.transpose(0, 2, 3, 1)
        return self._cisd_two_rdm

    @property
    def ccsd_one_rdm(self):
        r"""A 2-dimension array for CCSD one-particle density matrix in the MO
        representation.  d[p,q] = < a^\dagger_p a_q >
        """
        ccsd = self._pyscf_data.get('ccsd', None)
        if ccsd is None:
            return None

        mf = self._pyscf_data.get('scf', None)
        if isinstance(mf, scf.uhf.UHF):
            raise ValueError('Spin trace for UCCSD density matrix.')

        rdm1 = ccsd.make_rdm1()
        if isinstance(mf, scf.rohf.ROHF):
            rdm1 = rdm1[0] + rdm1[1]
        return rdm1.T

    @property
    def ccsd_two_rdm(self):
        r"""A 4-dimension array for CCSD two-particle density matrix in the MO
        representation.  D[p,q,r,s] = < a^\dagger_p a^\dagger_q a_r a_s >
        """
        ccsd = self._pyscf_data.get('ccsd', None)
        if ccsd is None:
            return None

        mf = self._pyscf_data.get('scf', None)
        if isinstance(mf, scf.uhf.UHF):
            raise ValueError('Spin trace for UCCSD density matrix.')

        rdm2 = ccsd.make_rdm2()
        if isinstance(mf, scf.rohf.ROHF):
            aa, ab, bb = rdm2
            rdm2 = aa + bb + ab + ab.transpose(2, 3, 0, 1)
        return rdm2.transpose(0, 2, 3, 1)

    @property
    def mp2_one_rdm(self):
        r"""A 2-dimension array for MP2 one-particle density matrix in the MO
        representation.  d[p,q] = < a^\dagger_p a_q >
        """
        mp2 = self._pyscf_data.get('mp2', None)
        if mp2 is None:
            return None

        mf = self._pyscf_data.get('scf', None)
        if isinstance(mf, scf.uhf.UHF):
            raise ValueError('Spin trace for UMP2 density matrix.')

        rdm1 = mp2.make_rdm1()
        if isinstance(mf, scf.rohf.ROHF):
            rdm1 = rdm1[0] + rdm1[1]
        return rdm1.T

    @property
    def mp2_two_rdm(self):
        r"""A 4-dimension array for MP2 two-particle density matrix in the MO
        representation.  D[p,q,r,s] = < a^\dagger_p a^\dagger_q a_r a_s >
        """
        mp2 = self._pyscf_data.get('mp2', None)
        if mp2 is None:
            return None

        mf = self._pyscf_data.get('scf', None)
        if isinstance(mf, scf.uhf.UHF):
            raise ValueError('Spin trace for UMP2 density matrix.')

        rdm2 = mp2.make_rdm2()
        if isinstance(mf, scf.rohf.ROHF):
            aa, ab, bb = rdm2
            rdm2 = aa + bb + ab + ab.transpose(2, 3, 0, 1)
        return rdm2.transpose(0, 2, 3, 1)

    @property
    def fci_one_rdm(self):
        r"""A 2-dimension array for FCI one-particle density matrix in the MO
        representation.  d[p,q] = < a^\dagger_p a_q >
        """
        if self._fci_one_rdm is None:
            fci = self._pyscf_data.get('fci', None)
            if fci is None:
                return None

            mf = self._pyscf_data.get('scf', None)
            if isinstance(mf, scf.uhf.UHF):
                raise ValueError('Spin trace for UHF-FCI density matrices.')

            norb = self.canonical_orbitals.shape[1]
            nelec = self.n_electrons
            self._fci_one_rdm = fci.make_rdm1(fci.ci, norb, nelec).T
        return self._fci_one_rdm

    @property
    def fci_two_rdm(self):
        r"""A 4-dimension array for FCI two-particle density matrix in the MO
        representation.  D[p,q,r,s] = < a^\dagger_p a^\dagger_q a_r a_s >
        """
        if self._fci_two_rdm is None:
            fci = self._pyscf_data.get('fci', None)
            if fci is None:
                return None

            mf = self._pyscf_data.get('scf', None)
            if isinstance(mf, scf.uhf.UHF):
                raise ValueError('Spin trace for UHF-FCI density matrix.')

            norb = self.canonical_orbitals.shape[1]
            nelec = self.n_electrons
            fci_rdm2 = fci.make_rdm2(fci.ci, norb, nelec)
            self._fci_two_rdm = fci_rdm2.transpose(0, 2, 3, 1)
        return self._fci_two_rdm

    @property
    def ccsd_single_amps(self):
        r"""A 2-dimension array t[a,i] for CCSD single excitation amplitudes
        where a is virtual index and i is occupied index.
        """
        if self._ccsd_single_amps is None:
            ccsd = self._pyscf_data.get('ccsd', None)
            if ccsd is None:
                return None

            t1 = spatial2spin(ccsd.t1)
            no, nv = t1.shape
            nmo = no + nv
            self._ccsd_single_amps = numpy.zeros((nmo, nmo))
            self._ccsd_single_amps[no:,:no] = t1.T

        return self._ccsd_single_amps

    @property
    def ccsd_double_amps(self):
        r"""A 4-dimension array t[a,i,b,j] for CCSD double excitation amplitudes
        where a, b are virtual indices and i, j are occupied indices.
        """
        if self._ccsd_double_amps is None:
            ccsd = self._pyscf_data.get('ccsd', None)
            if ccsd is None:
                return None

            t2 = spatial2spin(ccsd.t2)
            no, nv = t2.shape[1:3]
            nmo = no + nv
            self._ccsd_double_amps = numpy.zeros((nmo, nmo, nmo, nmo))
            self._ccsd_double_amps[no:,:no,no:,:no] = .5 * t2.transpose(2,0,3,1)

        return self._ccsd_double_amps
