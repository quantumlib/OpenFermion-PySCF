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

from functools import reduce

import numpy
from openfermion.hamiltonians import MolecularData
from pyscf import ao2mo


class MolecularData(MolecularData):
    def __init__(self, geometry=None, basis=None, multiplicity=None,
                 charge=0, description="", filename="", data_directory=None):
        super(MolecularData, self).__init__(
            geometry, basis, multiplicity,
            charge, description, filename, data_directory)
        self._pyscf_data = {}

    @property
    def canonical_orbitals(self):
        if self._canonical_orbitals is None:
            scf = self._pyscf_data.get('scf', None)
            self._canonical_orbitals = scf.mo_coeff
        return self._canonical_orbitals

    @property
    def overlap_integrals(self):
        if self._overlap_integrals is None:
            scf = self._pyscf_data.get('scf', None)
            self._overlap_integrals = scf.get_ovlp()
        return self._overlap_integrals

    @property
    def one_body_integrals(self):
        if self._one_body_integrals is None:
            scf = self._pyscf_data.get('scf', None)
            mo = self.canonical_orbitals
            h_core = scf.get_hcore()
            self._one_body_integrals = reduce(numpy.dot, (mo.T, h_core, mo))
        return self._one_body_integrals

    @property
    def two_body_integrals(self):
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
        if self._cisd_one_rdm is None:
            cisd = self._pyscf_data.get('cisd', None)
# pyscf one_rdm is computed as dm1[p,q] = <a^\dagger_q a_p>
            self._cisd_one_rdm = cisd.make_rdm1().T
        return self._cisd_one_rdm

    @property
    def cisd_two_rdm(self):
        if self._cisd_two_rdm is None:
            cisd = self._pyscf_data.get('cisd', None)
# pyscf.ci.cisd.make_rdm2 convention
#       dm2[p,s,q,r] = <a^\dagger_p a^\dagger_q a_r a_s>.
# the two_body_tensor in openfermion.ops._interaction_rdm.InteractionRDM
#       tbt[p,q,r,s] = <a^\dagger_p a^\dagger_q a_r a_s>.
            self._cisd_two_rdm = cisd.make_rdm2().transpose(0, 2, 3, 1)
        return self._cisd_two_rdm

    @property
    def ccsd_one_rdm(self):
        ccsd = self._pyscf_data.get('ccsd', None)
        return ccsd.make_rdm1().T

    @property
    def ccsd_two_rdm(self):
        ccsd = self._pyscf_data.get('ccsd', None)
        return ccsd.make_rdm2().transpose(0, 2, 3, 1)

    @property
    def fci_one_rdm(self):
        if self._fci_one_rdm is None:
            fci = self._pyscf_data.get('fci', None)
            norb = self.canonical_orbitals.shape[1]
            nelec = self.n_electrons
            self._fci_one_rdm = fci.make_rdm1(fci.ci, norb, nelec).T
        return self._fci_one_rdm

    @property
    def fci_two_rdm(self):
        if self._fci_two_rdm is None:
            fci = self._pyscf_data.get('fci', None)
            norb = self.canonical_orbitals.shape[1]
            nelec = self.n_electrons
            fci_rdm2 = fci.make_rdm2(fci.ci, norb, nelec)
            self._fci_two_rdm = fci_rdm2.transpose(0, 2, 3, 1)
        return self._fci_two_rdm

    @property
    def ccsd_single_amps(self):
        if self._ccsd_single_amps is None:
            ccsd = self._pyscf_data.get('ccsd', None)
            self._ccsd_single_amps = ccsd.t1
        return self._ccsd_single_amps.T

    @property
    def ccsd_double_amps(self):
        if self._ccsd_double_amps is None:
            ccsd = self._pyscf_data.get('ccsd', None)
            self._ccsd_double_amps = ccsd.t2
        return self._ccsd_double_amps.transpose(2, 0, 3, 1)
